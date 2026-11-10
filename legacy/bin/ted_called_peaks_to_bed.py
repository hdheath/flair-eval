#!/usr/bin/env python
"""
ted_called_peaks_to_bed.py — extract TED's GLOBAL TSS/TTS peak CALLS from a
production ted_log.tsv into IGV-loadable BED9 tracks (the ACTUAL peaks TED
called, not a re-derivation).

WHICH ROWS ARE THE CALLED PEAKS (confirmed in
flair-ted-52826/src/flair/ted.py::ted_global_collect_approved_peaks):
  Every detected global pile-up peak is logged once at ted.py:1209-1222 with
    stage  = 'global_peak'                                  (ted.py:1212)
    status = 'pass'   if reality >= threshold  (APPROVED -> became an iso end)
           = 'reject' if reality <  threshold  (no-filter diagnostic; kept)   (ted.py:1213)
    junc_id    = f'{chrom}:None-None:{strand}:0j'           (ted.py:1210)
    cluster_id = f'global_{end_type}_{peak_pos}'            (ted.py:1211)
    tss_pos = peak_pos, tts_pos = -1   for a TSS peak       (ted.py:1217-1218)
    tts_pos = peak_pos, tss_pos = -1   for a TTS peak
  Per-peak reality lives in TED_tss_reality (TSS rows) or TED_tts_reality
  (TTS rows); the other end's reality column is -1 (ted.py:1204-1205).

  The OTHER two ted_log writers are NOT scored peak calls and are skipped:
    - stage='global_assign' status='kept' (ted.py:1532): self-position read
      diagnostic (tss_pos/tts_pos default -1) — read-assignment placeholder.
    - stage='global_assign' status='pass' (ted.py:1607): per-isoform snapped
      (tss,tts) read-assignment rows (real chrom:start-end:strand:Nj junc_id,
      BOTH ends set) — the reads assigned to peaks, repeated per isoform.

The dRNA TTS force-pass is implemented as a lowered threshold (θ_tts=0.50 vs
0.84 for cDNA), so the in-row `status` already encodes it correctly per chem;
no special-casing here — we faithfully emit what production computed.

DEDUP: peaks at the same (chrom, pos, strand) within an end_type are deduped
keeping the max-reality instance (which is the approved one when any passed,
since pass requires reality >= θ).

Output (BED9, IGV-loadable): <outdir>/<lib>.tss_peaks.bed and
<lib>.tts_peaks.bed.  score = round(reality*1000) clamped 0..1000.
itemRgb: approved=0,170,0 (green); rejected=200,0,0 (red).
name = "<APPROVED|REJECTED>|r=<reality>|d=<depth>|m=<model>|a=<annot>|n=<n_reads>",
where reality = w_d*d + w_m*m + w_a*a (the three RAW component scores, each in
[0,1]) so the name shows WHY a peak passed/failed, not just the blended reality:
  d = TED_tss_depth / TED_tts_depth  (pile-up depth percentile score, S_d)
  m = TED_tss_model / TED_tts_model  (sequence/XGBoost artifact score, S_m)
  a = TED_tss_annot / TED_tts_annot  (annotation-proximity score, S_a; 0 for novel)
"""
import argparse
import os
import sys

# ted_log column indices (0-based), per _TED_LOG_HEADER (ted.py:404-413)
C_JUNC = 0
C_CLUSTER = 1
C_STAGE = 2
C_STATUS = 3
C_NREADS = 4
C_TSS_POS = 6
C_TTS_POS = 7
C_TSS_REALITY = 11
C_TTS_REALITY = 12
C_TSS_DEPTH = 13
C_TTS_DEPTH = 14
C_TSS_MODEL = 16
C_TTS_MODEL = 17
C_TSS_ANNOT = 18
C_TTS_ANNOT = 19
C_THR_TSS = 22
C_THR_TTS = 23
NCOL = 26

RGB_APPROVED = "0,170,0"
RGB_REJECTED = "200,0,0"


def parse_junc(junc_id):
    """'chrom:None-None:strand:0j' -> (chrom, strand). None on malformed."""
    parts = junc_id.split(":")
    if len(parts) < 4:
        return None, None
    chrom = parts[0]
    strand = parts[2]
    if strand not in ("+", "-"):
        return None, None
    return chrom, strand


def extract(log_path, end_type):
    """Stream the ted_log; collect deduped called peaks for one end_type.

    Returns dict keyed by (chrom, pos, strand) ->
    (reality, status, n_reads, depth, model, annot), keeping the max-reality
    instance per key (approved beats rejected because pass requires
    reality >= threshold).
    """
    pos_col = C_TSS_POS if end_type == "tss" else C_TTS_POS
    real_col = C_TSS_REALITY if end_type == "tss" else C_TTS_REALITY
    depth_col = C_TSS_DEPTH if end_type == "tss" else C_TTS_DEPTH
    model_col = C_TSS_MODEL if end_type == "tss" else C_TTS_MODEL
    annot_col = C_TSS_ANNOT if end_type == "tss" else C_TTS_ANNOT
    best = {}
    n_rows = n_pass = n_reject = n_skip_malformed = 0
    thresholds = set()
    with open(log_path) as fh:
        for line in fh:
            if not line or line[0] == "#":
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < NCOL:
                continue
            if f[C_STAGE] != "global_peak":
                # also silently drops repeated header rows (f[C_STAGE]=='stage')
                continue
            # this end_type's position must be set (the other is -1)
            try:
                pos = int(f[pos_col])
            except ValueError:
                continue
            if pos < 0:
                continue
            chrom, strand = parse_junc(f[C_JUNC])
            if chrom is None:
                n_skip_malformed += 1
                continue
            status = f[C_STATUS]
            try:
                reality = float(f[real_col])
            except ValueError:
                reality = -1.0
            try:
                n_reads = int(f[C_NREADS])
            except ValueError:
                n_reads = 0
            def _fcol(idx):
                try:
                    return float(f[idx])
                except (ValueError, IndexError):
                    return -1.0
            depth = _fcol(depth_col)
            model = _fcol(model_col)
            annot = _fcol(annot_col)
            thr = f[C_THR_TSS] if end_type == "tss" else f[C_THR_TTS]
            thresholds.add(thr)
            n_rows += 1
            if status == "pass":
                n_pass += 1
            elif status == "reject":
                n_reject += 1
            key = (chrom, pos, strand)
            prev = best.get(key)
            if prev is None or reality > prev[0]:
                best[key] = (reality, status, n_reads, depth, model, annot)
    return best, dict(
        n_rows=n_rows, n_pass=n_pass, n_reject=n_reject,
        n_skip_malformed=n_skip_malformed, thresholds=sorted(thresholds),
    )


def write_bed(best, end_type, lib, out_path):
    """Write deduped peaks to BED9. Returns (n_approved, n_rejected)."""
    label = end_type.upper()
    rows = []
    n_app = n_rej = 0
    for (chrom, pos, strand), (reality, status, n_reads, depth, model, annot) in best.items():
        approved = (status == "pass")
        if approved:
            n_app += 1
        else:
            n_rej += 1
        # IGV score 0..1000 from reality (reality is typically 0..1; -1 sentinel
        # would only occur if a column was malformed -> clamp to 0).
        sc = int(round(max(0.0, min(1.0, reality)) * 1000))
        rgb = RGB_APPROVED if approved else RGB_REJECTED
        # name carries reality AND its three raw components (d/m/a) so IGV hover
        # shows the full decomposition: reality = w_d*d + w_m*m + w_a*a.
        name = "%s|%s|r=%.4f|d=%.3f|m=%.3f|a=%.3f|n=%d" % (
            label, "APPROVED" if approved else "REJECTED",
            reality, depth, model, annot, n_reads)
        start = pos          # BED is 0-based; tss_pos/tts_pos are already 0-based genomic
        end = pos + 1
        rows.append((chrom, start, end, name, sc, strand, start, end, rgb))
    rows.sort(key=lambda r: (r[0], r[1], r[5]))
    with open(out_path, "w") as out:
        out.write('track name="%s.%s_peaks" itemRgb="On"\n' % (lib, end_type))
        for r in rows:
            out.write("%s\t%d\t%d\t%s\t%d\t%s\t%d\t%d\t%s\n" % r)
    return n_app, n_rej


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--ted-log", required=True)
    ap.add_argument("--lib", required=True, help="library label for filenames/track names")
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    summary = {}
    for end_type in ("tss", "tts"):
        best, stats = extract(args.ted_log, end_type)
        out_path = os.path.join(args.outdir, "%s.%s_peaks.bed" % (args.lib, end_type))
        n_app, n_rej = write_bed(best, end_type, args.lib, out_path)
        summary[end_type] = dict(
            path=out_path, raw_rows=stats["n_rows"],
            raw_pass=stats["n_pass"], raw_reject=stats["n_reject"],
            dedup_approved=n_app, dedup_rejected=n_rej,
            dedup_total=n_app + n_rej, thresholds=stats["thresholds"],
            malformed=stats["n_skip_malformed"],
        )

    print("=== %s ===" % args.lib)
    for end_type in ("tss", "tts"):
        s = summary[end_type]
        print("[%s] %s" % (end_type.upper(), s["path"]))
        print("    raw global_peak rows: %d (pass=%d reject=%d) | thresholds=%s | malformed_junc=%d"
              % (s["raw_rows"], s["raw_pass"], s["raw_reject"], s["thresholds"], s["malformed"]))
        print("    deduped peaks: %d  (APPROVED=%d  REJECTED=%d)"
              % (s["dedup_total"], s["dedup_approved"], s["dedup_rejected"]))
    sys.stdout.flush()


if __name__ == "__main__":
    main()

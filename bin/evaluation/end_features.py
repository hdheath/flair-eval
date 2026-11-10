#!/usr/bin/env python3
"""Post-assembly end features for TED-tagged isoform BEDs — the inputs for end filtering
that do NOT need the assembler re-run.

Three signals, all computable from (isoforms.bed + read map + the aligned BAM + genome):

  3' DOWNSTREAM A-RICHNESS (internal-priming proxy). Fraction of A in the N bp of GENOMIC
  sequence immediately downstream of the called 3' end (strand-aware: revcomp on minus), plus
  the longest A homopolymer in that window, plus the best canonical PAS hexamer found in
  [-50, 0] upstream. An oligo-dT read that mispriIsmed inside an A-rich stretch shows high
  a_frac with no PAS; a real cleavage site shows the opposite.

  5' UNTEMPLATED G (PacBio cap signature). In a TSO prep, RT adds untemplated C's at the cap,
  so the molecule reads [TSO]GGG[transcript 5' ...]. The G run sits in the soft clip IMMEDIATELY
  adjacent to the aligned block — for a forward alignment the TRAILING bases of the leading clip,
  for a reverse alignment the LEADING bases of the trailing clip read as C (pysam reports
  query_sequence in reference orientation). Reported per end as the percentage of that end's
  reads with EXACTLY one untemplated G (capg1) and with >=1 (capg_any).

  5' CLIP LENGTH (ONT TSO remnant). ONT PCS114 reads keep 20-30 bp of adapter/SSP remnant at the
  5' clip while PacBio reads are trimmed to 1-4 nt, so the clip LENGTH distribution is itself the
  ONT strand-switch signal. Reported per end as median/mean clip length and the fraction of reads
  whose clip is >= --tso-min-len.

Conventions follow ~/.claude-science-scratch/{capg,tso}.py so numbers are comparable to the
13_08_26 end-support analyses. Reads are attributed to an end when their own 5' (or 3') terminus
is within --window bp of it, so a read only votes on the end it actually reached.

Two tables are written: one per isoform, one per distinct END (isoforms sharing a TSS/TTS are
pooled — "percent of the reads in a peak" is an end-level quantity, not an isoform-level one).
"""
import argparse, collections, os, sys

COMP = {"A": "T", "C": "G", "G": "C", "T": "A", "N": "N"}
PAS = ("AATAAA", "ATTAAA", "AGTAAA", "TATAAA", "CATAAA", "GATAAA", "AATATA",
       "AATACA", "AATAGA", "ACTAAA", "AACAAA", "AAGAAA")


def rc(s):
    return "".join(COMP.get(c, "N") for c in reversed(s))


def clip5(aln):
    """The 5' soft clip in READ orientation ('' when the end is flush, None when unusable)."""
    q, ct = aln.query_sequence, aln.cigartuples
    if q is None or not ct:
        return None
    if not aln.is_reverse:
        return q[:ct[0][1]].upper() if ct[0][0] == 4 else ""
    return rc(q[len(q) - ct[-1][1]:].upper()) if ct[-1][0] == 4 else ""


def g_run(clip):
    """Untemplated G run abutting the alignment."""
    g = 0
    for ch in reversed(clip):
        if ch == "G":
            g += 1
        else:
            break
    return g


def parse_bed(path):
    """name -> (chrom, strand, tss, tts). Strand-aware; BED is 0-based half-open."""
    iso = {}
    with open(path) as fh:
        for ln in fh:
            f = ln.rstrip("\n").split("\t")
            if len(f) < 12:
                continue
            c, s, e, name, strand = f[0], int(f[1]), int(f[2]), f[3], f[5]
            tss, tts = (s, e) if strand == "+" else (e, s)
            iso[name] = (c, strand, tss, tts)
    return iso


def parse_map(path):
    """read name -> isoform name (FLAIR writes 'isoform<TAB>read1,read2,...')."""
    r2i = {}
    with open(path) as fh:
        for ln in fh:
            p = ln.rstrip("\n").split("\t")
            if len(p) < 2 or not p[1]:
                continue
            for r in p[1].split(","):
                if r:
                    r2i[r] = p[0]
    return r2i


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--bed", required=True, help="TED .isoforms.bed (tagged or filtered)")
    ap.add_argument("--read-map", required=True, help=".isoform.read.map.txt from --generate_map")
    ap.add_argument("--bam", required=True, help="the aligned BAM the assembly was built from")
    ap.add_argument("--genome", required=True, help="genome FASTA (indexed)")
    ap.add_argument("--out-prefix", required=True)
    ap.add_argument("--a-window", type=int, default=20,
                    help="bp downstream of the 3' end scored for A-richness (default 20)")
    ap.add_argument("--window", type=int, default=50,
                    help="bp: a read votes on an end when its own terminus is this close (default 50)")
    ap.add_argument("--tso-min-len", type=int, default=8,
                    help="5' clip length counted as a TSO remnant for the ONT fraction (default 8)")
    ap.add_argument("--chroms", nargs="*", default=None, help="limit to these chromosomes")
    args = ap.parse_args()

    import pysam

    iso = parse_bed(args.bed)
    if args.chroms:
        keep = set(args.chroms)
        iso = {k: v for k, v in iso.items() if v[0] in keep}
    if not iso:
        sys.exit("no isoforms parsed from --bed (check --chroms)")
    r2i = parse_map(args.read_map)
    sys.stderr.write(f"[end_features] {len(iso):,} isoforms, {len(r2i):,} mapped reads\n")

    # ---- 3' genomic features (genome only) -------------------------------
    fa = pysam.FastaFile(args.genome)
    lens = dict(zip(fa.references, fa.lengths))
    three = {}
    for name, (c, strand, tss, tts) in iso.items():
        n = args.a_window
        if c not in lens:
            continue
        if strand == "+":
            s, e = tts, min(tts + n, lens[c])
            down = fa.fetch(c, s, e).upper()
            up = fa.fetch(c, max(0, tts - 50), tts).upper()
        else:
            s, e = max(0, tts - n), tts
            down = rc(fa.fetch(c, s, e).upper())
            up = rc(fa.fetch(c, tts, min(tts + 50, lens[c])).upper())
        if not down:
            continue
        run = best = 0
        for ch in down:
            run = run + 1 if ch == "A" else 0
            best = max(best, run)
        motif = ""
        for m in PAS:
            if m in up:
                motif = m
                break
        three[name] = (down.count("A") / len(down), best, motif)

    # ---- 5' read features (one pass over the BAM) ------------------------
    agg = collections.defaultdict(lambda: [0, 0, 0, [], 0])  # n, capg_any, capg1, cliplens, n_3p
    bam = pysam.AlignmentFile(args.bam, "rb")
    chroms = args.chroms or sorted({v[0] for v in iso.values()})
    for c in chroms:
        try:
            it = bam.fetch(c)
        except (ValueError, KeyError):
            continue
        for aln in it:
            if aln.is_secondary or aln.is_supplementary or aln.is_unmapped:
                continue
            name = r2i.get(aln.query_name)
            if name is None:
                continue
            rec = iso.get(name)
            if rec is None or rec[0] != c:
                continue
            _, strand, tss, tts = rec
            r5 = aln.reference_start if not aln.is_reverse else aln.reference_end
            r3 = aln.reference_end if not aln.is_reverse else aln.reference_start
            a = agg[name]
            if abs(r5 - tss) <= args.window:
                clip = clip5(aln)
                if clip is not None:
                    a[0] += 1
                    g = g_run(clip)
                    if g >= 1:
                        a[1] += 1
                    if g == 1:
                        a[2] += 1
                    a[3].append(len(clip))
            if abs(r3 - tts) <= args.window:
                a[4] += 1

    def med(v):
        if not v:
            return ""
        v = sorted(v)
        n = len(v)
        return v[n // 2] if n % 2 else (v[n // 2 - 1] + v[n // 2]) / 2

    iso_out = args.out_prefix + ".isoform_features.tsv"
    end_out = args.out_prefix + ".end_features.tsv"
    tss_pool = collections.defaultdict(lambda: [0, 0, 0, []])
    tts_pool = collections.defaultdict(lambda: [0, [], [], []])

    with open(iso_out, "w") as fh:
        fh.write("isoform\tchrom\tstrand\ttss\ttts\tn_reads_5p\tcapg_any_pct\tcapg1_pct\t"
                 "clip5_median\tclip5_tso_pct\tn_reads_3p\ta%d_frac\ta%d_maxrun\tpas_motif\n"
                 % (args.a_window, args.a_window))
        for name, (c, strand, tss, tts) in sorted(iso.items()):
            n, ga, g1, clips, n3 = agg.get(name, [0, 0, 0, [], 0])
            af, mr, mo = three.get(name, ("", "", ""))
            tso_pct = (100.0 * sum(1 for x in clips if x >= args.tso_min_len) / n) if n else ""
            fh.write("\t".join(str(x) for x in [
                name, c, strand, tss, tts, n,
                f"{100.0*ga/n:.1f}" if n else "", f"{100.0*g1/n:.1f}" if n else "",
                med(clips), f"{tso_pct:.1f}" if n else "", n3,
                f"{af:.3f}" if af != "" else "", mr, mo]) + "\n")
            k5 = (c, strand, tss)
            p = tss_pool[k5]
            p[0] += n; p[1] += ga; p[2] += g1; p[3].extend(clips)
            k3 = (c, strand, tts)
            q = tts_pool[k3]
            q[0] += n3
            if af != "":
                q[1].append(af); q[2].append(mr); q[3].append(mo)

    with open(end_out, "w") as fh:
        fh.write("axis\tchrom\tstrand\tpos\tn_reads\tcapg_any_pct\tcapg1_pct\tclip5_median\t"
                 "clip5_tso_pct\ta%d_frac\ta%d_maxrun\tpas_motif\n" % (args.a_window, args.a_window))
        for (c, st, pos), (n, ga, g1, clips) in sorted(tss_pool.items()):
            tso_pct = (100.0 * sum(1 for x in clips if x >= args.tso_min_len) / n) if n else ""
            fh.write("\t".join(str(x) for x in [
                "5p", c, st, pos, n,
                f"{100.0*ga/n:.1f}" if n else "", f"{100.0*g1/n:.1f}" if n else "",
                med(clips), f"{tso_pct:.1f}" if n else "", "", "", ""]) + "\n")
        for (c, st, pos), (n, afs, runs, motifs) in sorted(tts_pool.items()):
            fh.write("\t".join(str(x) for x in [
                "3p", c, st, pos, n, "", "", "", "",
                f"{max(afs):.3f}" if afs else "", max(runs) if runs else "",
                next((m for m in motifs if m), "")]) + "\n")
    sys.stderr.write(f"[end_features] wrote {iso_out} and {end_out}\n")


if __name__ == "__main__":
    main()

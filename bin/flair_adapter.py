#!/usr/bin/env python3
"""flair_adapter.py: run any FLAIR version through one interface (docs/CONTRACTS.md section 2).

It runs *inside* the FLAIR environment and is the only thing in flair-eval that calls `flair`.
It is stdlib-only, because that environment can be any FLAIR version on any Python >= 3.8.

Usage (contract):
  flair_adapter.py align        --reads R --genome G [--gtf A] --threads N --out DIR
  flair_adapter.py assemble     --bam B --reads R --genome G --gtf A --args "..." --threads N --out DIR [--regions BED]
  flair_adapter.py quantify     --isoforms-fa F --reads R --sample S --threads N --out DIR
  flair_adapter.py productivity --isoforms-bed B --genome G --gtf A --out DIR
  flair_adapter.py detect       [--json]          (print the detected CLI profile; not a pipeline step)

Optional extras (not needed by core; see docs/flair_versions.md):
  quantify:     --isoforms-bed B --bam B --genome G   (FLAIR dev quantify counts from a genome BAM)
  every step:   --flair-bin PATH   --keep-work

The FLAIR version is detected by probing `flair --version` and the `--help` texts, so a branch or
commit is handled by what its CLI offers, not by its version string. FLAIR >= 3 runs
`flair transcriptome`; 2.x runs `correct` + `collapse`. A step the version cannot run writes
`skipped.json` and exits 0. FLAIR is never patched; version quirks are listed in
docs/flair_versions.md.

Every step writes `adapter_log.json` in DIR (exact commands, FLAIR version, exit codes, wall time).
"""
import argparse
import datetime
import gzip
import json
import os
import re
import shlex
import shutil
import subprocess
import sys
import time

ADAPTER_VERSION = "1.0.0"
LOG_SCHEMA = "flair-eval/adapter_log/1"
SKIP_SCHEMA = "flair-eval/skipped/1"
PRODUCTIVITY_CLASSES = ("PRO", "PTC", "NST", "NGO")
# Fixed internal sample/condition/batch for FLAIR quantify manifests. FLAIR 3.0.1 rejects '_'
# in these, and the normalised counts.tsv does not carry them, so the real sample name is only
# recorded in adapter_log.json.
QUANT_SID, QUANT_COND, QUANT_BATCH = "s1", "c1", "b1"
# The flair4_ted / CLAUDE.md rule: without these, TED byregion deadlocks.
THREAD_ENV = {"OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1", "MKL_NUM_THREADS": "1",
              "XGBOOST_NTHREAD": "1"}


class AdapterError(Exception):
    """Bad input or an unusable FLAIR install (exit code 2)."""


class StepFailed(Exception):
    """A FLAIR command failed (exit code = FLAIR's, or 1)."""

    def __init__(self, msg, code=1):
        super().__init__(msg)
        self.code = code or 1


class Skip(Exception):
    """This FLAIR version cannot run this step: write skipped.json and exit 0."""


def now_iso():
    return datetime.datetime.now(datetime.timezone.utc).replace(microsecond=0).isoformat()


def open_text(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path)


def child_env(extra=None):
    env = dict(os.environ)
    env.update(THREAD_ENV)
    if extra:
        env.update(extra)
    return env


# --------------------------------------------------------------------------------------------
# Run log
# --------------------------------------------------------------------------------------------
class RunLog:
    """Collects everything adapter_log.json records and runs commands."""

    def __init__(self, step, out_dir, args_ns):
        self.step = step
        self.out_dir = out_dir
        self.logs_dir = os.path.join(out_dir, "logs")
        self.t0 = time.time()
        self.data = {
            "schema": LOG_SCHEMA,
            "adapter_version": ADAPTER_VERSION,
            "step": step,
            "status": "running",
            "started": now_iso(),
            "finished": None,
            "wall_seconds": None,
            "argv": sys.argv[:],
            "inputs": {k: v for k, v in vars(args_ns).items() if k not in ("func",)},
            "flair_version": None,
            "flair_bin": None,
            "profile": None,
            "commands": [],
            "warnings": [],
            "outputs": {},
            "notes": [],
        }

    def warn(self, msg):
        self.data["warnings"].append(msg)
        sys.stderr.write("flair_adapter.py: warning: %s\n" % msg)

    def note(self, msg):
        self.data["notes"].append(msg)

    def run(self, argv, cwd, label, env=None, stdout_path=None, check=True):
        """Run one command, log it, raise StepFailed on a non-zero exit when check is set."""
        os.makedirs(self.logs_dir, exist_ok=True)
        n = len(self.data["commands"]) + 1
        log_rel = os.path.join("logs", "%02d_%s.log" % (n, re.sub(r"[^A-Za-z0-9_.-]", "_", label)))
        log_path = os.path.join(self.out_dir, log_rel)
        entry = {"label": label, "argv": [str(a) for a in argv],
                 "cmd": " ".join(shlex.quote(str(a)) for a in argv), "cwd": cwd,
                 "log": log_rel, "stdout": None, "exit_code": None, "wall_seconds": None}
        self.data["commands"].append(entry)
        t = time.time()
        with open(log_path, "w") as logf:
            logf.write("# %s\n# cwd: %s\n" % (entry["cmd"], cwd))
            logf.flush()
            out_fh = open(stdout_path, "w") if stdout_path else logf
            try:
                proc = subprocess.run([str(a) for a in argv], cwd=cwd, env=env or child_env(),
                                      stdout=out_fh, stderr=logf)
                code = proc.returncode
            except OSError as e:
                logf.write("\n# could not start: %s\n" % e)
                code = 127
            finally:
                if stdout_path:
                    out_fh.close()
        if stdout_path:
            entry["stdout"] = os.path.relpath(stdout_path, self.out_dir) \
                if stdout_path.startswith(self.out_dir) else stdout_path
        entry["exit_code"] = code
        entry["wall_seconds"] = round(time.time() - t, 3)
        if check and code != 0:
            tail = ""
            try:
                with open(log_path) as fh:
                    tail = "".join(fh.readlines()[-15:])
            except OSError:
                pass
            raise StepFailed("%s exited %d (log %s)\n%s" % (label, code, log_path, tail), code)
        return code

    def finish(self, status):
        self.data["status"] = status
        self.data["finished"] = now_iso()
        self.data["wall_seconds"] = round(time.time() - self.t0, 3)
        os.makedirs(self.out_dir, exist_ok=True)
        with open(os.path.join(self.out_dir, "adapter_log.json"), "w") as fh:
            json.dump(self.data, fh, indent=2, default=str)
            fh.write("\n")


# --------------------------------------------------------------------------------------------
# FLAIR CLI detection
# --------------------------------------------------------------------------------------------
LONG_FLAG = re.compile(r"(?<![\w-])(--[A-Za-z][\w-]*)")


def help_flags(text):
    return sorted(set(LONG_FLAG.findall(text or "")))


def usage_block(text):
    """The argparse usage block: from 'usage:' to the first blank line."""
    m = re.search(r"usage:(.*?)(?:\n\s*\n|\Z)", text or "", re.S)
    return m.group(1) if m else ""


def strip_brackets(s):
    out, depth = [], 0
    for ch in s:
        if ch == "[":
            depth += 1
        elif ch == "]":
            depth = max(0, depth - 1)
        elif depth == 0:
            out.append(ch)
    return "".join(out)


def required_long_flags(text):
    """Long names of required options, also when the usage line shows only the short form."""
    req = set(LONG_FLAG.findall(strip_brackets(usage_block(text))))
    shorts = set(re.findall(r"(?<![\w-])(-[A-Za-z])\b", strip_brackets(usage_block(text))))
    for short in shorts:
        # help lines look like "  -b GENOME_ALIGNED_BAM, --genome_aligned_bam ..." or "-b, --x"
        m = re.search(r"^\s*%s\b[^\n]*?(--[A-Za-z][\w-]*)" % re.escape(short), text or "", re.M)
        if m:
            req.add(m.group(1))
    # argparse groups titled "required ..." (TED quantify lists --isoform_bed there while its
    # usage line shows it in brackets)
    for m in re.finditer(r"^required[^\n]*:\n((?:[ \t]+[^\n]*\n?)+)", text or "", re.M):
        for line in m.group(1).splitlines():
            mm = re.match(r"^\s{1,4}(?:-[A-Za-z][^,\n]*,\s*)?(--[A-Za-z][\w-]*)", line)
            if mm:
                req.add(mm.group(1))
    return sorted(req)


def parse_subcommands(text):
    text = text or ""
    m = re.search(r"\{([a-zA-Z_,]+)\}", text)              # 3.0.x argparse choices
    if m:
        return sorted(set(m.group(1).split(",")))
    m = re.search(r"module:\s*([a-zA-Z_, ]+)", text)       # 2.x hand-written usage
    if m:
        return sorted({s.strip() for s in m.group(1).split(",") if s.strip()})
    subs = set()                                           # dev: indented "name   description"
    in_block = False
    for line in text.splitlines():
        if re.match(r"^subcommands:", line):
            in_block = True
            continue
        if in_block:
            mm = re.match(r"^\s{4}([a-z][a-zA-Z_]*)\s{2,}\S", line)
            if mm:
                subs.add(mm.group(1))
            elif line.strip() and not line.startswith(" "):
                in_block = False
    return sorted(subs)


def parse_version(text):
    for line in (text or "").splitlines()[::-1]:
        m = re.search(r"FLAIR\s+v?(\S+)", line)
        if m:
            return m.group(1)
    m = re.search(r"(\d+\.\d+\S*)", text or "")
    return m.group(1) if m else None


def version_major(v):
    m = re.match(r"(\d+)", v or "")
    return int(m.group(1)) if m else None


def capture(argv, env=None):
    try:
        p = subprocess.run(argv, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                           env=env or child_env(), universal_newlines=True, timeout=600)
        return p.returncode, p.stdout, p.stderr
    except (OSError, subprocess.TimeoutExpired) as e:
        return 127, "", str(e)


def find_flair(flair_bin=None):
    cand = flair_bin or os.environ.get("FLAIR_EVAL_FLAIR") or shutil.which("flair")
    if not cand or not os.path.exists(cand):
        raise AdapterError("no `flair` executable found (looked at --flair-bin, $FLAIR_EVAL_FLAIR, PATH)")
    return os.path.abspath(cand)


def env_tool(flair, name):
    """A tool from the FLAIR env (next to `flair`), else from PATH."""
    side = os.path.join(os.path.dirname(flair), name)
    if os.path.exists(side) and os.access(side, os.X_OK):
        return side
    return shutil.which(name)


def detect_profile(flair):
    """Probe the FLAIR CLI and return what the adapter needs to know about it."""
    code, out, err = capture([flair, "--version"])
    version = parse_version(out) or parse_version(err)
    if version is None:
        raise AdapterError("`%s --version` gave no version (exit %d): %s" % (flair, code, (out + err)[-500:]))
    _, out, err = capture([flair, "--help"])
    subs = parse_subcommands(out + "\n" + err)
    prof = {"flair_version": version, "flair_major": version_major(version), "subcommands": subs}

    def sub_help(name):
        if name not in subs:
            return None
        c, o, e = capture([flair, name, "--help"])
        return o + "\n" + e

    tx = sub_help("transcriptome")
    major = prof["flair_major"]
    if tx is not None and (major is None or major >= 3):
        prof["assemble_cli"] = "transcriptome"
        flags = help_flags(tx)
        req = required_long_flags(tx)
        prof["transcriptome"] = {
            "flags": flags,
            "required": req,
            "bam_flag": "-b",
            "needs_sample_name": "--sample_name" in req,
            "parallel_flag": "--parallel_mode" if "--parallel_mode" in flags else (
                "--parallelmode" if "--parallelmode" in flags else None),
            "temp_dir": "--temp_dir" in flags,
            "generate_map": "--generate_map" in flags,
        }
    elif "collapse" in subs and "correct" in subs:
        prof["assemble_cli"] = "correct+collapse"
        ch, co = sub_help("correct"), sub_help("collapse")
        prof["correct"] = {"flags": help_flags(ch), "required": required_long_flags(ch)}
        cflags = help_flags(co)
        prof["collapse"] = {"flags": cflags, "required": required_long_flags(co),
                            "temp_dir": "--temp_dir" in cflags,
                            "generate_map": "--generate_map" in cflags}
    else:
        prof["assemble_cli"] = None

    prof["align"] = "align" in subs

    qh = sub_help("quantify")
    if qh is None:
        prof["quantify_mode"] = None
    else:
        qflags = help_flags(qh)
        qreq = required_long_flags(qh)
        mode = "bam_manifest" if "--manifest" in qflags else "reads_manifest"
        prof["quantify_mode"] = mode
        prof["quantify"] = {"flags": qflags, "required": qreq,
                            "needs_isoform_bed": "--isoform_bed" in qreq or mode == "bam_manifest",
                            "temp_dir": "--temp_dir" in qflags}

    pp = env_tool(flair, "predictProductivity")
    if pp:
        _, o, e = capture([pp, "--help"])
        pflags = help_flags(o + "\n" + e)
        prof["predict_productivity"] = {"bin": pp, "flags": pflags,
                                        "append_column": "--append_column" in pflags}
    else:
        prof["predict_productivity"] = None
    return prof


# --------------------------------------------------------------------------------------------
# File helpers
# --------------------------------------------------------------------------------------------
def require_file(path, what):
    if not path:
        raise AdapterError("%s is required" % what)
    if not os.path.exists(path):
        raise AdapterError("%s not found: %s" % (what, path))
    return os.path.abspath(path)


def read_bed(path, min_cols=12):
    rows = []
    with open_text(path) as fh:
        for ln, line in enumerate(fh, 1):
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < min_cols:
                raise AdapterError("%s line %d: %d columns, need >= %d" % (path, ln, len(f), min_cols))
            rows.append(f)
    return rows


def bed12_blocks(f):
    start = int(f[1])
    sizes = [int(x) for x in f[10].rstrip(",").split(",") if x != ""]
    starts = [int(x) for x in f[11].rstrip(",").split(",") if x != ""]
    if len(sizes) != int(f[9]) or len(starts) != int(f[9]):
        raise AdapterError("BED row %s: blockCount %s does not match blocks" % (f[3], f[9]))
    return [(start + s, start + s + z) for s, z in zip(starts, sizes)]


def norm_bed12(f):
    """First 12 columns, integers checked, score clamped to 0..1000, no trailing commas."""
    out = f[:12]
    try:
        score = int(float(out[4])) if out[4] not in ("", ".") else 0
    except ValueError:
        score = 0
    out[4] = str(max(0, min(1000, score)))
    out[10] = out[10].rstrip(",")
    out[11] = out[11].rstrip(",")
    for i in (1, 2, 6, 7, 9):
        int(out[i])
    if out[5] not in ("+", "-", "."):
        raise AdapterError("BED row %s: bad strand %r" % (out[3], out[5]))
    if out[8] in ("", "."):
        out[8] = "0"
    bed12_blocks(out)
    return out


def gtf_attr(attrs, key):
    m = re.search(r'(?:^|;)\s*%s\s+"?([^";]*)"?' % re.escape(key), attrs)
    return m.group(1).strip() if m else None


def gtf_transcript_genes(path):
    """{transcript_id: gene_id} from a GTF (transcript or exon lines)."""
    out = {}
    if not path or not os.path.exists(path):
        return out
    with open_text(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] not in ("transcript", "exon"):
                continue
            tid, gid = gtf_attr(f[8], "transcript_id"), gtf_attr(f[8], "gene_id")
            if tid and gid and tid not in out:
                out[tid] = gid
    return out


MISSING_GENE = {"", ".", "NA", "None", "none", "nan"}


def assign_genes(names, native_rows, native_gtf):
    """gene_id per isoform name. Order: native BED col 13 (dev/TED) -> native GTF (matched as
    '<transcript_id>_<gene_id>' == BED name, the 2.x/3.0.1 naming, or transcript_id == name) ->
    the last '_' token of the name -> novel_gene_<n>."""
    col13 = {}
    for f in native_rows:
        if len(f) >= 13 and f[12].strip() not in MISSING_GENE:
            col13[f[3]] = f[12].strip()
    t2g = gtf_transcript_genes(native_gtf)
    by_pair = {"%s_%s" % (t, g): g for t, g in t2g.items()}
    genes, sources, n_novel = {}, {"bed_col13": 0, "gtf": 0, "name_suffix": 0, "novel": 0}, 0
    for name in names:
        g, src = col13.get(name), "bed_col13"
        if g is None:
            g, src = by_pair.get(name), "gtf"
        if g is None and name in t2g:
            g, src = t2g[name], "gtf"
        if g is None and "_" in name:
            tail = name.rsplit("_", 1)[1]
            if tail not in MISSING_GENE:
                g, src = tail, "name_suffix"
        if g is None or g in MISSING_GENE:
            n_novel += 1
            g, src = "novel_gene_%d" % n_novel, "novel"
        genes[name] = g
        sources[src] += 1
    return genes, sources


def write_gtf(rows, genes, path, source="FLAIR"):
    with open(path, "w") as out:
        for f in rows:
            name, chrom, strand = f[3], f[0], f[5]
            gid = genes[name]
            blocks = bed12_blocks(f)
            attrs = 'gene_id "%s"; transcript_id "%s";' % (gid, name)
            out.write("\t".join([chrom, source, "transcript", str(int(f[1]) + 1), f[2], ".",
                                 strand, ".", attrs]) + "\n")
            order = list(range(len(blocks)))
            if strand == "-":
                order = order[::-1]
            for k, i in enumerate(order, 1):
                s, e = blocks[i]
                out.write("\t".join([chrom, source, "exon", str(s + 1), str(e), ".", strand, ".",
                                     attrs + ' exon_number "%d";' % k]) + "\n")


def read_fasta(path):
    seqs, name, buf = {}, None, []
    with open_text(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(buf)
                name, buf = line[1:].split()[0] if line[1:].split() else "", []
            elif name is not None:
                buf.append(line.strip())
    if name is not None:
        seqs[name] = "".join(buf)
    return seqs


def fasta_names(path):
    names = []
    with open_text(path) as fh:
        for line in fh:
            if line.startswith(">"):
                parts = line[1:].split()
                names.append(parts[0] if parts else "")
    return names


def first_existing(paths):
    for p in paths:
        if p and os.path.exists(p):
            return p
    return None


def native_sibling(paths, fname="isoforms.flair.bed"):
    """Find `fname` next to any of `paths` (also behind symlinks, as Nextflow stages inputs)."""
    for p in paths:
        if not p:
            continue
        for d in (os.path.dirname(os.path.abspath(p)), os.path.dirname(os.path.realpath(p))):
            c = os.path.join(d, fname)
            if os.path.exists(c):
                return c
    return None


def split_user_args(arg_str, known_flags, log, what):
    toks = shlex.split(arg_str or "")
    unknown = [t.split("=", 1)[0] for t in toks if t.startswith("--")
               and t.split("=", 1)[0] not in known_flags]
    if unknown:
        log.warn("%s: flags not in this FLAIR's `%s --help`: %s (FLAIR 3.0.1 ignores unknown "
                 "flags silently; stricter versions fail)" % (what, what, " ".join(unknown)))
    return toks


def has_flag(toks, *names):
    for t in toks:
        if t.split("=", 1)[0] in names:
            return True
    return False


def safe_name(s):
    return re.sub(r"[^A-Za-z0-9.-]", "-", s or "") or "sample"


# --------------------------------------------------------------------------------------------
# Steps
# --------------------------------------------------------------------------------------------
def step_align(a, log, prof, flair, work):
    if not prof.get("align"):
        raise Skip("this FLAIR (%s) has no `flair align` subcommand" % prof["flair_version"])
    reads = require_file(a.reads, "--reads")
    genome = require_file(a.genome, "--genome")
    if a.gtf:
        log.note("--gtf is accepted for interface symmetry; FLAIR align defaults do not use it")
    prefix = os.path.join(work, "aligned")
    log.run([flair, "align", "-g", genome, "-r", reads, "-t", str(a.threads), "-o", prefix],
            cwd=work, label="flair_align")
    bam = prefix + ".bam"
    if not os.path.exists(bam):
        raise StepFailed("flair align wrote no %s" % bam)
    shutil.move(bam, os.path.join(a.out, "aligned.bam"))
    bai_src = first_existing([bam + ".bai", prefix + ".bai"])
    if bai_src:
        shutil.move(bai_src, os.path.join(a.out, "aligned.bam.bai"))
    else:
        samtools = env_tool(flair, "samtools")
        if not samtools:
            raise StepFailed("flair align wrote no index and samtools is not available")
        log.run([samtools, "index", os.path.join(a.out, "aligned.bam")], cwd=work, label="samtools_index")
    log.data["outputs"] = {"aligned.bam": os.path.getsize(os.path.join(a.out, "aligned.bam"))}


def subset_bam(bam, regions, work, log, flair):
    samtools = env_tool(flair, "samtools")
    if not samtools:
        raise AdapterError("--regions needs samtools in the FLAIR env")
    sub = os.path.join(work, "regions.bam")
    log.run([samtools, "view", "-b", "-L", regions, "-o", sub, bam], cwd=work, label="samtools_view_regions")
    log.run([samtools, "index", sub], cwd=work, label="samtools_index_regions")
    return sub


def bam_read_names(bam, work, log, flair):
    samtools = env_tool(flair, "samtools")
    out = os.path.join(work, "regions.read_names.txt")
    log.run([samtools, "view", bam], cwd=work, label="samtools_view_names", stdout_path=out)
    names = set()
    with open(out) as fh:
        for line in fh:
            if line and not line.startswith("@"):
                names.add(line.split("\t", 1)[0])
    return names


def subset_reads(reads, keep, out_path):
    """Write the reads in `keep` (FASTA or FASTQ, gz or not) to out_path; returns count kept."""
    kept = 0
    with open_text(reads) as fh, open(out_path, "w") as out:
        first = fh.readline()
        if not first:
            return 0
        fastq = first.startswith("@")
        line = first
        if fastq:
            while line:
                rec = [line, fh.readline(), fh.readline(), fh.readline()]
                if rec[0][1:].split()[0] in keep:
                    out.writelines(rec)
                    kept += 1
                line = fh.readline()
        else:
            write = False
            while line:
                if line.startswith(">"):
                    write = line[1:].split()[0] in keep
                    kept += write
                if write:
                    out.write(line)
                line = fh.readline()
    return kept


def find_native_outputs(prefix):
    return {
        "bed": first_existing([prefix + ".isoforms.bed"]),
        "gtf": first_existing([prefix + ".isoforms.gtf"]),
        "fa": first_existing([prefix + ".isoforms.fa"]),
        "map": first_existing([prefix + ".isoform.read.map.txt", prefix + ".isoform.read.map.tsv",
                               prefix + ".isoforms.read.map.txt"]),
        "counts": first_existing([prefix + ".isoform.counts.txt", prefix + ".isoform.counts.tsv"]),
    }


def run_transcriptome(a, log, prof, flair, work, bam, genome, gtf, user_toks):
    tp = prof["transcriptome"]
    prefix = os.path.join(work, "flair")
    cmd = [flair, "transcriptome", "-b", bam, "-g", genome, "-f", gtf, "-t", str(a.threads), "-o", prefix]
    if tp["needs_sample_name"] and not has_flag(user_toks, "--sample_name"):
        cmd += ["--sample_name", safe_name(a.sample)]
    if tp["parallel_flag"] and not has_flag(user_toks, "--parallel_mode", "--parallelmode"):
        cmd += [tp["parallel_flag"], "byregion"]      # CLAUDE.md: byregion, never bychrom
    if tp["generate_map"] and not has_flag(user_toks, "--generate_map"):
        cmd += ["--generate_map"]
    if tp["temp_dir"] and not has_flag(user_toks, "--temp_dir"):
        tmp = os.path.join(work, "tmp")
        os.makedirs(tmp, exist_ok=True)
        cmd += ["--temp_dir", tmp]
    cmd += user_toks
    log.run(cmd, cwd=work, label="flair_transcriptome")
    return prefix


def run_correct_collapse(a, log, prof, flair, work, bam, reads, genome, gtf, user_toks):
    # 2.x: reads BED12 from the shared BAM, with FLAIR's own converter (keeps flag 0/16 only,
    # as 2.x `flair align` does), then correct -> collapse.
    py = env_tool(flair, "python") or env_tool(flair, "python3") or sys.executable
    reads_bed = os.path.join(work, "reads.bed")
    log.run([py, "-m", "flair.bam2Bed12", "-i", bam], cwd=work, label="bam2Bed12", stdout_path=reads_bed)
    corr = os.path.join(work, "correct")
    cmd = [flair, "correct", "-q", reads_bed, "-f", gtf, "-g", genome, "-t", str(a.threads), "-o", corr]
    log.run(cmd, cwd=work, label="flair_correct")
    corrected = first_existing([corr + "_all_corrected.bed"])
    if not corrected:
        raise StepFailed("flair correct wrote no %s_all_corrected.bed" % corr)
    cp = prof["collapse"]
    prefix = os.path.join(work, "flair")
    cmd = [flair, "collapse", "-g", genome, "-r", reads, "-q", corrected, "-f", gtf,
           "-t", str(a.threads), "-o", prefix]
    if cp["generate_map"] and not has_flag(user_toks, "--generate_map"):
        cmd += ["--generate_map"]
    if cp["temp_dir"] and not has_flag(user_toks, "--temp_dir"):
        tmp = os.path.join(work, "tmp")
        os.makedirs(tmp, exist_ok=True)
        cmd += ["--temp_dir", tmp + os.sep]
    cmd += user_toks
    log.run(cmd, cwd=work, label="flair_collapse")
    return prefix


def step_assemble(a, log, prof, flair, work):
    if prof.get("assemble_cli") is None:
        raise Skip("this FLAIR (%s) has neither `transcriptome` nor `correct`+`collapse`" % prof["flair_version"])
    bam = require_file(a.bam, "--bam")
    genome = require_file(a.genome, "--genome")
    gtf = require_file(a.gtf, "--gtf")
    reads = require_file(a.reads, "--reads") if a.reads else None
    if not os.path.exists(bam + ".bai") and not os.path.exists(re.sub(r"\.bam$", ".bai", bam)):
        log.warn("no .bai next to %s; FLAIR needs an indexed BAM" % bam)
    if a.regions:
        regions = require_file(a.regions, "--regions")
        bam = subset_bam(bam, regions, work, log, flair)
        log.note("BAM restricted to --regions with samtools view -L (reads overlapping the regions)")
    if prof["assemble_cli"] == "transcriptome":
        user_toks = split_user_args(a.args, prof["transcriptome"]["flags"], log, "transcriptome")
        prefix = run_transcriptome(a, log, prof, flair, work, bam, genome, gtf, user_toks)
    else:
        if not reads:
            raise AdapterError("FLAIR %s assembles with collapse, which needs --reads" % prof["flair_version"])
        user_toks = split_user_args(a.args, prof["collapse"]["flags"], log, "collapse")
        if a.regions:
            keep = bam_read_names(bam, work, log, flair)
            sub_reads = os.path.join(work, "regions.reads" + (".fq" if fastx_is_fastq(reads) else ".fa"))
            n = subset_reads(reads, keep, sub_reads)
            log.note("reads restricted to the %d reads in the region BAM for collapse" % n)
            reads = sub_reads
        prefix = run_correct_collapse(a, log, prof, flair, work, bam, reads, genome, gtf, user_toks)
    normalise_assembly(a, log, find_native_outputs(prefix))


def fastx_is_fastq(path):
    with open_text(path) as fh:
        return fh.read(1) == "@"


def normalise_assembly(a, log, nat):
    if not nat["bed"]:
        raise StepFailed("FLAIR wrote no isoforms.bed")
    out = a.out
    native_dir = os.path.join(out, "native")
    os.makedirs(native_dir, exist_ok=True)
    native_rows = read_bed(nat["bed"])
    ncols = sorted({len(f) for f in native_rows})
    log.data["outputs"]["native_bed_columns"] = ncols
    # native BED kept verbatim: dev quantify and inline productivity need FLAIR's own columns.
    shutil.copyfile(nat["bed"], os.path.join(out, "isoforms.flair.bed"))
    if nat["gtf"]:
        shutil.copyfile(nat["gtf"], os.path.join(native_dir, "flair.isoforms.gtf"))
    if nat["counts"]:
        shutil.copyfile(nat["counts"], os.path.join(native_dir, os.path.basename(nat["counts"])))

    seen, rows, dups = set(), [], []
    for f in native_rows:
        if f[3] in seen:
            dups.append(f[3])
            continue
        seen.add(f[3])
        rows.append(norm_bed12(f))
    if dups:
        log.warn("%d duplicate isoform names in FLAIR's BED; kept the first row of each: %s"
                 % (len(dups), ",".join(dups[:5])))
    rows.sort(key=lambda f: (f[0], int(f[1]), int(f[2]), f[3]))
    names = [f[3] for f in rows]
    with open(os.path.join(out, "isoforms.bed"), "w") as fh:
        for f in rows:
            fh.write("\t".join(f) + "\n")

    genes, sources = assign_genes(names, native_rows, nat["gtf"])
    write_gtf(rows, genes, os.path.join(out, "isoforms.gtf"))
    log.data["outputs"]["gene_id_sources"] = sources

    if not nat["fa"]:
        raise StepFailed("FLAIR wrote no isoforms.fa")
    seqs = read_fasta(nat["fa"])
    missing = [n for n in names if n not in seqs]
    extra = len(set(seqs) - set(names))
    if missing:
        log.warn("%d BED isoforms have no FASTA record: %s" % (len(missing), ",".join(missing[:5])))
    if extra:
        log.warn("%d FASTA records are not in the BED and were dropped" % extra)
    with open(os.path.join(out, "isoforms.fa"), "w") as fh:
        for n in names:
            if n in seqs:
                fh.write(">%s\n%s\n" % (n, seqs[n]))

    n_map = 0
    if nat["map"]:
        name_set, dropped = set(names), 0
        with open_text(nat["map"]) as fh, open(os.path.join(out, "read_map.tsv"), "w") as outm:
            outm.write("read_id\tisoform_id\n")
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if len(f) < 2 or (f[0] == "name" and f[1] == "reads"):
                    continue
                if f[0] not in name_set:
                    dropped += 1
                    continue
                for r in f[1].split(","):
                    if r:
                        outm.write("%s\t%s\n" % (r, f[0]))
                        n_map += 1
        if dropped:
            log.warn("%d read-map rows name isoforms that are not in the BED (dropped)" % dropped)
    else:
        log.note("this FLAIR wrote no read map; read_map.tsv omitted (optional in the contract)")
    log.data["outputs"].update({"isoforms": len(rows), "genes": len(set(genes.values())),
                                "read_map_rows": n_map, "fasta_records": len(names) - len(missing)})


def native_bed_for(a, log, *paths):
    """FLAIR's own isoform BED for quantify/productivity: --isoforms-bed if it has FLAIR's extra
    columns, else isoforms.flair.bed written next to the assemble outputs."""
    if getattr(a, "isoforms_bed", None) and os.path.exists(a.isoforms_bed):
        with open_text(a.isoforms_bed) as fh:
            for line in fh:
                if line.strip() and not line.startswith(("#", "track")):
                    if len(line.split("\t")) >= 13:
                        return os.path.abspath(a.isoforms_bed)
                    break
    return native_sibling([getattr(a, "isoforms_bed", None)] + list(paths))


def step_quantify(a, log, prof, flair, work):
    mode = prof.get("quantify_mode")
    if mode is None:
        raise Skip("this FLAIR (%s) has no `flair quantify`" % prof["flair_version"])
    fa = require_file(a.isoforms_fa, "--isoforms-fa")
    log.note("sample %r is quantified as FLAIR sample id %s (names are not carried into counts.tsv)"
             % (a.sample, QUANT_SID))
    qp = prof["quantify"]
    prefix = os.path.join(work, "quant")
    man = os.path.join(work, "manifest.tsv")
    tmp = os.path.join(work, "tmp")
    os.makedirs(tmp, exist_ok=True)
    nbed = native_bed_for(a, log, fa) if qp["needs_isoform_bed"] else None
    if mode == "bam_manifest":
        missing = [f for f, v in (("--bam", a.bam), ("--genome", a.genome)) if not v]
        if missing or not nbed:
            what = missing + ([] if nbed else ["FLAIR's isoform BED (isoforms.flair.bed)"])
            raise Skip("this FLAIR (%s) quantifies from a genome BAM, which needs %s"
                       % (prof["flair_version"], ", ".join(what)))
        bam = require_file(a.bam, "--bam")
        genome = require_file(a.genome, "--genome")
        with open(man, "w") as fh:
            fh.write("%s\t%s\t%s\t%s\n" % (QUANT_SID, QUANT_COND, QUANT_BATCH, bam))
        cmd = [flair, "quantify", "--manifest", man, "-g", genome, "--isoform_bed", nbed,
               "-t", str(a.threads), "-o", prefix]
    else:
        reads = require_file(a.reads, "--reads")
        with open(man, "w") as fh:
            fh.write("%s\t%s\t%s\t%s\n" % (QUANT_SID, QUANT_COND, QUANT_BATCH, reads))
        cmd = [flair, "quantify", "-r", man, "-i", fa, "-t", str(a.threads), "-o", prefix]
        if qp["needs_isoform_bed"]:
            if not nbed:
                raise Skip("this FLAIR (%s) quantify requires --isoform_bed with FLAIR's own BED, "
                           "and none was given or found next to --isoforms-fa" % prof["flair_version"])
            cmd += ["--isoform_bed", nbed]
    if qp["temp_dir"]:
        cmd += ["--temp_dir", tmp + os.sep]
    log.run(cmd, cwd=work, label="flair_quantify")
    # 2.x/3.0.1/TED/dev d033600: <o>.counts.tsv; dev 389748e: <o>.isoforms.counts.tsv
    counts_file = first_existing([prefix + s for s in (".counts.tsv", ".isoforms.counts.tsv",
                                                       ".counts_matrix.tsv", ".counts.txt")])
    if not counts_file:
        base = os.path.basename(prefix)
        cands = sorted(f for f in os.listdir(work) if f.startswith(base + ".") and "counts" in f
                       and "tpm" not in f and f.endswith((".tsv", ".txt")))
        if cands:
            counts_file = os.path.join(work, cands[0])
            log.warn("counts file found by name pattern: %s" % cands[0])
    if not counts_file:
        raise StepFailed("flair quantify wrote no %s.counts.tsv" % prefix)
    shutil.copyfile(counts_file, os.path.join(a.out, "flair.counts.tsv"))
    counts = parse_counts(counts_file)
    names = fasta_names(fa)
    write_counts(names, counts, os.path.join(a.out, "counts.tsv"), log)


def parse_counts(path):
    """{isoform_id: count} from FLAIR's counts matrix (first sample column)."""
    out = {}
    with open_text(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        id_col = header.index("isoform_id") if "isoform_id" in header else 0
        val_col = len(header) - 1
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) <= max(id_col, val_col):
                continue
            out[f[id_col]] = out.get(f[id_col], 0.0) + float(f[val_col])
    return out


def write_counts(names, counts, path, log):
    name_set = set(names)
    unmatched = {k: v for k, v in counts.items() if k not in name_set}
    if unmatched:
        # FLAIR 2.x/3.0.1 can name a counts row '<isoform>_<gene>' when the FASTA has '<isoform>'
        fixed = {}
        for k, v in unmatched.items():
            base = k.rsplit("_", 1)[0]
            if base in name_set and base not in counts:
                fixed[base] = v
        if fixed:
            counts = {k: v for k, v in counts.items() if k in name_set}
            counts.update(fixed)
        still = [k for k in unmatched if k.rsplit("_", 1)[0] not in fixed]
        if still:
            log.warn("%d counts rows name no isoform in the FASTA (dropped): %s"
                     % (len(still), ",".join(still[:5])))
    zero = [n for n in names if n not in counts]
    if zero:
        log.note("%d isoforms had no counts row (FLAIR omits zero-read isoforms); written as 0" % len(zero))
    total = sum(counts.get(n, 0.0) for n in names)
    with open(path, "w") as fh:
        fh.write("isoform_id\tcount\tcpm\n")
        for n in names:
            c = counts.get(n, 0.0)
            cpm = c / total * 1e6 if total > 0 else 0.0
            fh.write("%s\t%s\t%.4f\n" % (n, fmt_num(c), cpm))
    log.data["outputs"].update({"isoforms": len(names), "assigned_reads": fmt_num(total),
                                "zero_count_isoforms": sum(1 for n in names if counts.get(n, 0.0) == 0)})


def fmt_num(x):
    return str(int(x)) if float(x).is_integer() else ("%.4f" % x)


def step_productivity(a, log, prof, flair, work):
    bed = require_file(a.isoforms_bed, "--isoforms-bed")
    names = [f[3] for f in read_bed(bed)]
    nbed = native_bed_for(a, log, bed)
    if nbed:
        inline = productivity_from_native(nbed)
        if inline and len(inline) >= max(1, int(0.5 * len(names))):
            log.note("productivity read from FLAIR's own isoform BED (column 18, computed by "
                     "flair transcriptome itself): %s" % nbed)
            log.data["outputs"]["method"] = "inline"
            write_productivity(names, inline, os.path.join(a.out, "productivity.tsv"), log)
            return
    pp = prof.get("predict_productivity")
    if not pp:
        raise Skip("this FLAIR (%s) has no predictProductivity and its isoform BED carries no "
                   "productivity column" % prof["flair_version"])
    if not pp["append_column"]:
        raise Skip("this FLAIR's predictProductivity has no --append_column and no inline "
                   "productivity was found")
    genome = require_file(a.genome, "--genome")
    gtf = require_file(a.gtf, "--gtf")
    # copy the BED under a plain name: TED's predictProductivity splits paths on the first '.bed'
    pdir = os.path.join(work, "productivity")
    os.makedirs(pdir, exist_ok=True)
    shutil.copyfile(bed, os.path.join(pdir, "isoforms.bed"))
    cmd = [pp["bin"], "-i", "isoforms.bed", "-g", gtf, "-f", genome, "-o", "pred", "--append_column"]
    # FLAIR's own pipeline (transcriptome --predictCDS, 2.x collapse) calls it with --longestORF,
    # and 2.x/3.0.1 require one of --firstTIS/--longestORF.
    if "--longestORF" in pp["flags"]:
        cmd.append("--longestORF")
    log.run(cmd, cwd=pdir, label="predictProductivity")
    outbed = first_existing([os.path.join(pdir, "pred.bed"), os.path.join(pdir, "pred.productivity.bed")])
    if not outbed:
        cands = [os.path.join(pdir, x) for x in os.listdir(pdir) if x.startswith("pred") and x.endswith(".bed")]
        outbed = cands[0] if cands else None
    if not outbed:
        raise StepFailed("predictProductivity wrote no BED")
    shutil.copyfile(outbed, os.path.join(a.out, "flair.productivity.bed"))
    calls = {}
    for f in read_bed(outbed, min_cols=13):
        cls = f[12].strip()
        calls[f[3]] = (f[6], f[7], cls)
    log.data["outputs"]["method"] = "predictProductivity"
    write_productivity(names, calls, os.path.join(a.out, "productivity.tsv"), log)


def productivity_from_native(path):
    calls = {}
    for f in read_bed(path):
        if len(f) >= 18 and f[17].strip() in PRODUCTIVITY_CLASSES:
            calls.setdefault(f[3], (f[6], f[7], f[17].strip()))
    return calls


def write_productivity(names, calls, path, log):
    bad = sorted({c[2] for c in calls.values() if c[2] not in PRODUCTIVITY_CLASSES})
    if bad:
        log.warn("unexpected productivity classes kept verbatim: %s" % ",".join(bad))
    missing = [n for n in names if n not in calls]
    if missing:
        log.warn("%d isoforms have no productivity call (written as NA): %s"
                 % (len(missing), ",".join(missing[:5])))
    counts = {}
    with open(path, "w") as fh:
        fh.write("isoform_id\tcds_start\tcds_end\tclass\n")
        for n in names:
            s, e, c = calls.get(n, ("NA", "NA", "NA"))
            counts[c] = counts.get(c, 0) + 1
            fh.write("%s\t%s\t%s\t%s\n" % (n, s, e, c))
    log.data["outputs"].update({"isoforms": len(names), "classes": counts})


STEPS = {"align": step_align, "assemble": step_assemble, "quantify": step_quantify,
         "productivity": step_productivity}


def run_step(a):
    out = os.path.abspath(a.out)
    a.out = out
    os.makedirs(out, exist_ok=True)
    for stale in ("skipped.json",):
        if os.path.exists(os.path.join(out, stale)):
            os.remove(os.path.join(out, stale))
    log = RunLog(a.step, out, a)
    work = os.path.join(out, ".work")
    status, code = "failed", 1
    try:
        flair = find_flair(a.flair_bin)
        log.data["flair_bin"] = flair
        prof = detect_profile(flair)
        log.data["profile"] = prof
        log.data["flair_version"] = prof["flair_version"]
        os.makedirs(work, exist_ok=True)
        STEPS[a.step](a, log, prof, flair, work)
        status, code = "ok", 0
    except Skip as e:
        status, code = "skipped", 0
        with open(os.path.join(out, "skipped.json"), "w") as fh:
            json.dump({"schema": SKIP_SCHEMA, "step": a.step, "reason": str(e),
                       "flair_version": log.data["flair_version"], "created": now_iso()}, fh, indent=2)
            fh.write("\n")
        sys.stderr.write("flair_adapter.py: %s skipped: %s\n" % (a.step, e))
    except AdapterError as e:
        log.data["error"] = str(e)
        sys.stderr.write("flair_adapter.py: error: %s\n" % e)
        code = 2
    except StepFailed as e:
        log.data["error"] = str(e)
        sys.stderr.write("flair_adapter.py: %s failed: %s\n" % (a.step, e))
        code = e.code
    finally:
        log.finish(status)
        if not a.keep_work and status != "failed" and os.path.isdir(work):
            shutil.rmtree(work, ignore_errors=True)
    return code


def run_detect(a):
    try:
        prof = detect_profile(find_flair(a.flair_bin))
    except AdapterError as e:
        sys.stderr.write("flair_adapter.py: error: %s\n" % e)
        return 2
    if a.json:
        json.dump(prof, sys.stdout, indent=2)
        sys.stdout.write("\n")
    else:
        for k in ("flair_version", "subcommands", "assemble_cli", "align", "quantify_mode"):
            print("%s\t%s" % (k, prof.get(k)))
        print("predict_productivity\t%s" % bool(prof.get("predict_productivity")))
    return 0


def positive_int(s):
    v = int(s)
    if v < 1:
        raise argparse.ArgumentTypeError("must be >= 1")
    return v


def build_parser():
    p = argparse.ArgumentParser(prog="flair_adapter.py", description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="step")

    def common(sp, threads=True):
        sp.add_argument("--out", required=True, help="output directory")
        if threads:
            sp.add_argument("--threads", type=positive_int, default=1)
        sp.add_argument("--flair-bin", default=None, help="flair executable (default: PATH)")
        sp.add_argument("--keep-work", action="store_true", help="keep DIR/.work (FLAIR's raw run)")

    s = sub.add_parser("align", help="flair align -> aligned.bam(.bai)")
    s.add_argument("--reads", required=True)
    s.add_argument("--genome", required=True)
    s.add_argument("--gtf", default=None)
    common(s)

    s = sub.add_parser("assemble", help="transcriptome (>=3) or correct+collapse (2.x)")
    s.add_argument("--bam", required=True)
    s.add_argument("--reads", default=None, help="needed by 2.x collapse")
    s.add_argument("--genome", required=True)
    s.add_argument("--gtf", required=True)
    s.add_argument("--args", default="", help="appended verbatim to transcriptome/collapse")
    s.add_argument("--regions", default=None, help="BED: restrict the BAM to these regions")
    s.add_argument("--sample", default="flaireval", help="--sample_name for FLAIR versions that need one")
    common(s)

    s = sub.add_parser("quantify", help="flair quantify -> counts.tsv")
    s.add_argument("--isoforms-fa", required=True)
    s.add_argument("--reads", default=None)
    s.add_argument("--sample", required=True)
    s.add_argument("--isoforms-bed", default=None, help="optional; FLAIR's own BED is found next to the FASTA")
    s.add_argument("--bam", default=None, help="optional; needed by FLAIR dev quantify")
    s.add_argument("--genome", default=None, help="optional; needed by FLAIR dev quantify")
    common(s)

    s = sub.add_parser("productivity", help="predictProductivity (or FLAIR's inline CDS) -> productivity.tsv")
    s.add_argument("--isoforms-bed", required=True)
    s.add_argument("--genome", required=True)
    s.add_argument("--gtf", required=True)
    common(s, threads=False)

    s = sub.add_parser("detect", help="print the detected FLAIR CLI profile")
    s.add_argument("--flair-bin", default=None)
    s.add_argument("--json", action="store_true")
    return p


def main(argv=None):
    p = build_parser()
    a = p.parse_args(argv)
    if not a.step:
        p.print_help()
        return 2
    if a.step == "detect":
        return run_detect(a)
    if a.step == "quantify" and not a.reads and not a.bam:
        p.error("quantify needs --reads (or --bam for FLAIR dev)")
    return run_step(a)


if __name__ == "__main__":
    sys.exit(main())

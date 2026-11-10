"""Command-line entry point: `flair-eval init | check | run | report`."""

from __future__ import annotations

import argparse
import glob
import json
import subprocess
import sys
from pathlib import Path

from . import DEFAULT_REF, PIPELINE, __version__, config, doctor, env_check
from .command import Command, Options, opt, build_command, files_command, quick_command, \
    report_command, split_extra
from .tokens import TokenError, parse_tokens, runs_from_tokens
from .validate import Issue, check_files, check_paths, errors, validate_bundle, \
    validate_runs, validate_samples
from .wizard import Abort, Wizard, write_configs

EPILOG = f"""New here? Three steps:
  flair-eval doctor      check Java, Nextflow, conda and disk space
  flair-eval init        answer a few questions; writes samples.csv and runs.yaml
  flair-eval run --outdir /big/storage/flair-eval-results

Try it on the bundled test data first: flair-eval run --test
Guide: https://github.com/{PIPELINE}/blob/v2/docs/usage_cli.md"""


# ---------------------------------------------------------------- argument parsing
def _quick_sample_args(p: argparse.ArgumentParser) -> None:
    g = p.add_argument_group("one sample from flags (instead of samples.csv)")
    x = g.add_mutually_exclusive_group()
    x.add_argument("--reads", help="reads FASTQ/FASTA(.gz)")
    x.add_argument("--bam", help="sorted, indexed BAM")
    g.add_argument("--genome", help="genome FASTA")
    g.add_argument("--gtf", help="reference annotation GTF")
    g.add_argument("--sample", help="sample name (default: reads/BAM file name)")
    g.add_argument("--platform", choices=["pacbio", "ont"])
    g.add_argument("--library", choices=["cDNA", "dRNA"])
    r = p.add_argument_group("runs from flags (instead of runs.yaml)")
    r.add_argument("--flair", help="comma-separated tokens, e.g. 3.0.1,dev,env:/path/to/env")
    r.add_argument("--flair_args", "--flair-args", default="",
                   help="flair transcriptome arguments for every --flair run")


def _launch_args(p: argparse.ArgumentParser) -> None:
    g = p.add_argument_group("launch")
    g.add_argument("-r", "--ref", default=DEFAULT_REF,
                   help=f"pipeline revision: branch, tag or SHA (default {DEFAULT_REF})")
    g.add_argument("--profile", "-profile",
                   help="Nextflow profile(s), e.g. slurm,conda (default: picked for this machine, "
                        "see `flair-eval doctor`)")
    g.add_argument("--outdir", default=None,
                   help="output folder (default ./flair-eval-results); use big storage")
    g.add_argument("--work-dir", "-w", dest="work_dir", help="Nextflow work directory")
    g.add_argument("--resume", "-resume", action="store_true", help="resume a previous launch")
    g.add_argument("--dry-run", "-n", action="store_true",
                   help="print the command without running it")
    g.add_argument("--skip-checks", action="store_true",
                   help="do not check that java and nextflow are installed")


def _pipeline_params(p: argparse.ArgumentParser) -> None:
    g = p.add_argument_group("pipeline options")
    g.add_argument("--baseline", help="run id that report differences are computed against")
    g.add_argument("--report_title", "--report-title", help='report heading')
    g.add_argument("--align_flair", "--align-flair",
                   help="FLAIR release for the shared align step (pipeline default 3.0.1)")
    g.add_argument("--regions", help="BED file to restrict assembly and evaluation")
    g.add_argument("--bundles", help="earlier bundle folders/globs to add to the report")
    g.add_argument("--env_cache", "--env-cache", help="where FLAIR envs are built and cached")


def make_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(
        prog="flair-eval",
        description="Compare FLAIR runs (versions, branches, settings) with the flair-eval "
                    "Nextflow pipeline.",
        epilog=EPILOG, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--version", action="version", version=f"flair-eval {__version__}")
    sub = ap.add_subparsers(dest="cmd", metavar="{doctor,init,check,run,report}")

    p = sub.add_parser("doctor", help="check this machine (java, nextflow, conda, disk) with fixes",
                       description="Check that this machine can run flair-eval, and suggest a "
                                   "Nextflow profile.")
    p.add_argument("--outdir", help="where results will go (checked for free space)")

    p = sub.add_parser("init", help="interactive wizard that writes samples.csv and runs.yaml",
                       description="Answer a few questions; writes samples.csv and runs.yaml. "
                                   "Give --reads/--bam, --genome, --gtf and --flair to skip "
                                   "the questions.")
    p.add_argument("--dir", default=".", help="where to write the files (default: here)")
    p.add_argument("--force", action="store_true", help="overwrite existing files")
    _quick_sample_args(p)

    p = sub.add_parser("check", help="validate samples.csv, runs.yaml and bundles",
                       description="Validate config files against the pipeline's schemas. "
                                   "Missing input files are warnings.")
    p.add_argument("--samples", "--input", dest="samples", help="samples.csv")
    p.add_argument("--runs", help="runs.yaml")
    p.add_argument("--bundles", nargs="+", help="bundle folders to check bundle.json in")
    p.add_argument("--no-files", action="store_true", help="skip the input-file checks")
    p.add_argument("--strict", action="store_true", help="treat warnings as errors")

    p = sub.add_parser("run", help="check the configs and launch the pipeline",
                       description="Launch `nextflow run hdheath/flair-eval`. Use samples.csv "
                                   "+ runs.yaml, or the one-command flags. Anything after "
                                   "`--` goes to nextflow unchanged.")
    p.add_argument("--samples", "--input", dest="samples", help="samples.csv")
    p.add_argument("--runs", help="runs.yaml")
    p.add_argument("--test", action="store_true",
                   help="run the bundled test data (-profile test); no config files needed")
    _quick_sample_args(p)
    _pipeline_params(p)
    _launch_args(p)

    p = sub.add_parser("report", help="build the report from existing bundles only",
                       description="Report-only mode: no FLAIR runs, just report.html from "
                                   "the bundles given.")
    p.add_argument("bundles", nargs="+", help="bundle folders or globs (e.g. out1 out2/bundles)")
    p.add_argument("--title", dest="report_title", help="report heading")
    p.add_argument("--baseline", help="run id that differences are computed against")
    _launch_args(p)
    return ap


DEFAULT_OUTDIR = "./flair-eval-results"
TEST_OUTDIR = "./flair-eval-test"


# ---------------------------------------------------------------- helpers
def _print_issues(issues: list[Issue], out=None) -> None:
    for i in issues:
        print(f"  {i}", file=out or sys.stderr)


def _quick_sample(a) -> dict | None:
    path = a.reads or a.bam
    if not (path or a.genome or a.gtf):
        return None
    s = {"sample": a.sample or (config.default_sample_name(path) if path else "")}
    for k in ("reads", "bam", "genome", "gtf", "platform", "library"):
        if getattr(a, k, None):
            s[k] = getattr(a, k)
    return config.clean_sample(s)


def _sample_flag_args(s: dict) -> list:
    args = ["--reads", s["reads"]] if s.get("reads") else ["--bam", s["bam"]]
    args += ["--genome", s["genome"], "--gtf", s["gtf"], "--sample", s["sample"]]
    for k in ("platform", "library"):
        if s.get(k):
            args += [f"--{k}", s[k]]
    return args


def _options(a, extra: list) -> Options:
    params = {k: getattr(a, k, None) for k in
              ("bundles", "baseline", "report_title", "align_flair", "regions", "env_cache")}
    return Options(ref=a.ref, profile=a.profile, outdir=a.outdir or DEFAULT_OUTDIR,
                   params=params,
                   resume=a.resume, work_dir=a.work_dir, extra=extra,
                   samples_path=getattr(a, "samples", None) or "samples.csv",
                   runs_path=getattr(a, "runs", None) or "runs.yaml")


def _auto_profile(a, extra: str = "") -> None:
    """Fill a.profile from this machine when --profile was not given."""
    if a.profile:
        if extra and extra not in a.profile.split(","):
            a.profile = f"{extra},{a.profile}"
        return
    found = doctor.suggest_profile()
    a.profile = ",".join(p for p in (extra, found) if p) or None
    if found:
        print(f"Using -profile {a.profile} (picked for this machine; set --profile to change).")


def _launch(cmd, a) -> int:
    print("Command:\n" + cmd.pretty() + "\n")
    if a.dry_run:
        return 0
    if not a.skip_checks:
        errs, warns = env_check.check_tools()
        for w in warns:
            print("warning: " + w, file=sys.stderr)
        if errs:
            for e in errs:
                print("error: " + e, file=sys.stderr)
            print("Fix the above, or use --dry-run to only print the command.", file=sys.stderr)
            return 2
    try:
        return subprocess.call(cmd.argv)
    except KeyboardInterrupt:
        return 130


# ---------------------------------------------------------------- commands
def cmd_init(a) -> int:
    d = Path(a.dir)
    if a.flair or a.reads or a.bam:
        s = _quick_sample(a)
        if not s or not a.flair:
            print("error: without questions, init needs --reads or --bam, --genome, --gtf "
                  "and --flair", file=sys.stderr)
            return 2
        try:
            runs = runs_from_tokens(a.flair, a.flair_args)
        except TokenError as e:
            print(f"error: {e}", file=sys.stderr)
            return 2
        samples = [s]
    else:
        try:
            answers = Wizard().run()
        except (Abort, KeyboardInterrupt) as e:
            print(f"\nStopped ({e or 'interrupted'}); nothing written.", file=sys.stderr)
            return 1
        if answers is None:
            print("\nNo files needed for the test data. Check this machine, then launch:")
            print("  flair-eval doctor")
            print("  flair-eval run --test")
            print("Results land in ./flair-eval-test (report: report/report.html).")
            return 0
        samples, runs = answers
    issues = errors(validate_samples(samples) + validate_runs({"runs": runs}))
    if issues:
        print("error: these answers do not make a valid config:", file=sys.stderr)
        _print_issues(issues)
        return 1
    try:
        sp, rp = write_configs(d, samples, runs, force=a.force)
    except FileExistsError as e:
        print(f"error: {e} already exists; use --force to overwrite or --dir elsewhere",
              file=sys.stderr)
        return 1
    print(f"\nWrote {sp} ({len(samples)} sample{'s' * (len(samples) != 1)}) and "
          f"{rp} ({len(runs)} run{'s' * (len(runs) != 1)}).")
    for w in check_paths(samples):
        print(f"  {w}")
    flags = [] if d == Path(".") else ["--samples", str(sp), "--runs", str(rp)]
    print("\nNext:")
    print("  flair-eval doctor")
    print("  " + " ".join(["flair-eval check"] + flags))
    print("  " + " ".join(["flair-eval run"] + flags + ["--outdir", "<folder on big storage>"]))
    cmd = build_command(samples, runs, Options(samples_path=str(sp), runs_path=str(rp)))
    print("\nOr run Nextflow directly:\n" + cmd.pretty())
    return 0


def cmd_check(a) -> int:
    samples, runs = a.samples, a.runs
    if samples is None and runs is None and not a.bundles:
        samples = "samples.csv" if Path("samples.csv").exists() else None
        runs = "runs.yaml" if Path("runs.yaml").exists() else None
        if samples is None and runs is None:
            print("error: nothing to check. Give --samples and/or --runs (or --bundles), "
                  "or run this where samples.csv/runs.yaml are.", file=sys.stderr)
            return 2
    issues = check_files(samples, runs, files=not a.no_files)
    for b in a.bundles or []:
        found = sorted(set(glob.glob(str(Path(b) / "**" / "bundle.json"), recursive=True)))
        if not found:
            issues.append(Issue(("bundles",), "missing_file", f"no bundle.json under {b}"))
        for f in found:
            try:
                doc = json.loads(Path(f).read_text())
            except ValueError as e:
                issues.append(Issue((f,), "json", f"not valid JSON: {e}"))
                continue
            issues += [Issue((f,) + i.path[1:], i.code, i.message, i.level)
                       for i in validate_bundle(doc)]
            print(f"checked {f}")
    for name, path in (("samples", samples), ("runs", runs)):
        if path is not None:
            bad = [i for i in issues if i.path[:1] == (name,) and i.level == "error"]
            print(f"{path}: {'OK' if not bad else f'{len(bad)} problem(s)'}")
    _print_issues(issues, sys.stdout)
    n_err = len(errors(issues))
    n_warn = len(issues) - n_err
    if n_err or (a.strict and n_warn):
        print(f"\n{n_err} error(s), {n_warn} warning(s).")
        return 1
    print("\nAll good." if not n_warn else f"\nNo errors, {n_warn} warning(s).")
    return 0


def files_command_test(opts: Options) -> Command:
    """`nextflow run ... -profile test,...`: the test profile sets --input and --runs."""
    argv = files_command(opts).argv
    i = argv.index("--input")
    del argv[i:i + 4]
    return Command(argv, "test", [])


def cmd_doctor(a) -> int:
    return doctor.report(doctor.run_checks(a.outdir), doctor.suggest_profile())


def cmd_run(a, extra: list) -> int:
    if a.test:
        if a.samples or a.runs or a.flair or _quick_sample(a) is not None:
            print("error: --test brings its own samples and runs; drop the other inputs",
                  file=sys.stderr)
            return 2
        _auto_profile(a, "test")
        a.outdir = a.outdir or TEST_OUTDIR
        return _launch(files_command_test(_options(a, extra)), a)
    quick = _quick_sample(a)
    if quick is not None and a.samples:
        print("error: use either --samples or --reads/--bam/--genome/--gtf, not both",
              file=sys.stderr)
        return 2
    if a.flair and a.runs:
        print("error: use either --runs or --flair, not both", file=sys.stderr)
        return 2
    if quick is None and not a.samples and Path("samples.csv").exists():
        a.samples = "samples.csv"
    if not a.flair and not a.runs and Path("runs.yaml").exists():
        a.runs = "runs.yaml"
    if quick is None and not a.samples:
        print("error: no samples. Give --samples samples.csv, or --reads/--bam with --genome "
              "and --gtf (or run `flair-eval init`).", file=sys.stderr)
        return 2
    if not a.flair and not a.runs:
        print("error: no runs. Give --runs runs.yaml or --flair 3.0.1,dev "
              "(or run `flair-eval init`).", file=sys.stderr)
        return 2

    issues: list[Issue] = []
    if quick is not None:
        issues += validate_samples([quick]) + check_paths([quick])
    else:
        issues += check_files(a.samples, None)
    run_ids: list[str] = []
    if a.flair:
        try:
            run_ids = [r["id"] for r in runs_from_tokens(a.flair, a.flair_args)]
        except TokenError as e:
            issues.append(Issue(("flair",), "token", str(e)))
    else:
        issues += check_files(None, a.runs)
        try:
            doc = config.read_runs(a.runs) if Path(a.runs).exists() else {}
            run_ids = [r.get("id") for r in (doc.get("runs") or []) if isinstance(r, dict)]
        except (config.YamlError, AttributeError):
            pass
    if a.baseline and run_ids and a.baseline not in run_ids:
        issues.append(Issue(("baseline",), "baseline",
                            f"--baseline {a.baseline!r} is not a run id ({', '.join(run_ids)})"))
    if issues:
        _print_issues(issues)
    if errors(issues):
        print("Fix the errors above first (`flair-eval check` explains them).", file=sys.stderr)
        return 1

    _auto_profile(a)
    opts = _options(a, extra)
    if quick is not None and a.flair:
        cmd = quick_command(_sample_flag_args(quick), a.flair, a.flair_args, opts)
    elif quick is None and a.runs:
        cmd = files_command(opts)
    else:
        # mixed: one source from flags, the other from a file -> write nothing, pass both
        argv = files_command(opts).argv
        i = argv.index("--input")
        if quick is not None:
            argv[i:i + 2] = _sample_flag_args(quick)
        else:
            j = argv.index("--runs")
            argv[j:j + 2] = ["--flair", a.flair] + (opt("--flair_args", a.flair_args)
                                                     if a.flair_args else [])
        cmd = Command(argv, "mixed", [])
    return _launch(cmd, a)


def cmd_report(a, extra: list) -> int:
    missing = [b for b in a.bundles if not any(ch in b for ch in "*?[") and not Path(b).exists()]
    for b in missing:
        print(f"warning: {b} does not exist here (fine if it is a path on the cluster)",
              file=sys.stderr)
    _auto_profile(a)
    opts = Options(ref=a.ref, profile=a.profile, outdir=a.outdir or DEFAULT_OUTDIR,
                   resume=a.resume,
                   work_dir=a.work_dir, extra=extra,
                   params={"baseline": a.baseline, "report_title": a.report_title})
    return _launch(report_command(a.bundles, opts), a)


#: Options whose value is itself flags (e.g. --flair_args "--stringent"). argparse would read
#: such a value as an option, so main() glues it on as --flair_args=--stringent.
DASH_VALUE_OPTS = ("--flair_args", "--flair-args")


def _glue_dash_values(argv: list) -> list:
    out, i = [], 0
    while i < len(argv):
        if argv[i] in DASH_VALUE_OPTS and i + 1 < len(argv):
            out.append(f"{argv[i]}={argv[i + 1]}")
            i += 2
        else:
            out.append(argv[i])
            i += 1
    return out


def main(argv: list | None = None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    argv, extra = split_extra(_glue_dash_values(argv))
    ap = make_parser()
    a = ap.parse_args(argv)
    if a.cmd is None:
        ap.print_help()
        return 0
    if extra and a.cmd not in ("run", "report"):
        ap.error("arguments after `--` are only for run and report")
    if a.cmd == "doctor":
        return cmd_doctor(a)
    if a.cmd == "init":
        return cmd_init(a)
    if a.cmd == "check":
        return cmd_check(a)
    if a.cmd == "run":
        return cmd_run(a, extra)
    return cmd_report(a, extra)


if __name__ == "__main__":  # pragma: no cover
    sys.exit(main())

"""`flair-eval init`: an interactive wizard that writes samples.csv and runs.yaml.

`ask` and `say` are injectable so tests can drive the wizard with scripted answers.
"""

from __future__ import annotations

import datetime
import getpass
import re
from pathlib import Path
from typing import Callable

from . import config, phoenix, presets
from .tokens import ID_RE, TokenError, parse_tokens, parse_token
from .validate import NAME_HINT

PLATFORMS = ("pacbio", "ont")
LIBRARIES = ("cDNA", "dRNA")


def guess_platform(path: str) -> str:
    """Default platform from a file name: 'pacbio', 'ont' or ''."""
    n = path.lower().rsplit("/", 1)[-1]
    if re.search(r"(^|[^a-z])(ont|nanopore|promethion|minion|gridion|r10|dorado|guppy)([^a-z]|$)", n):
        return "ont"
    if re.search(r"(pacbio|hifi|ccs|isoseq|iso-seq|kinnex|masseq|mas-seq|sequel|revio|flnc)", n):
        return "pacbio"
    return ""


def guess_library(path: str) -> str:
    n = path.lower().rsplit("/", 1)[-1]
    if re.search(r"(drna|direct.?rna|rna00[24]|sqk-rna)", n):
        return "dRNA"
    if re.search(r"(cdna|isoseq|iso-seq|kinnex|masseq|flnc|pcs1|sqk-pcs|sqk-lsk)", n):
        return "cDNA"
    return ""


class Abort(Exception):
    pass


def suggested_outdir(user: str | None = None, today: datetime.date | None = None) -> str:
    """Default --outdir offered on phoenix: groups storage, named dd_mm_yy_flair-eval."""
    try:
        user = user or getpass.getuser()
    except (KeyError, OSError):  # pragma: no cover - no passwd entry
        user = "$USER"
    stamp = (today or datetime.date.today()).strftime("%d_%m_%y")
    return f"{phoenix.LAB_ROOT}/{user}/flair-eval/{stamp}_flair-eval"


class Wizard:
    def __init__(self, ask: Callable[[str], str] | None = None,
                 say: Callable[[str], None] | None = None, on_phoenix: bool = False,
                 default_outdir: str | None = None):
        self._ask = ask or (lambda p: input(p))
        self.say = say or (lambda m: print(m))
        self.on_phoenix = on_phoenix
        self.default_outdir = default_outdir

    # ------------------------------------------------------------ prompting helpers
    def ask(self, prompt: str, default: str | None = None, required: bool = False,
            check: Callable[[str], str | None] | None = None) -> str:
        suffix = f" [{default}]" if default not in (None, "") else ""
        for _ in range(20):
            try:
                ans = self._ask(f"{prompt}{suffix}: ").strip()
            except EOFError as e:
                raise Abort("input ended") from e
            if not ans and default is not None:
                ans = default
            if not ans and required:
                self.say("  This one is required.")
                continue
            err = check(ans) if (check and ans) else None
            if err:
                self.say(f"  {err}")
                continue
            return ans
        raise Abort(f"no usable answer for: {prompt}")

    def yes(self, prompt: str, default: bool = False) -> bool:
        ans = self.ask(prompt + (" [Y/n]" if default else " [y/N]"), default="")
        return default if not ans else ans.lower().startswith("y")

    def choose(self, prompt: str, options: tuple, default: str = "") -> str:
        opts = "/".join(options)

        def check(a):
            return None if a in options or a == "-" else f"Use one of {opts}."
        hint = ", - to skip" if default else ", Enter to skip"
        ans = self.ask(f"{prompt} ({opts}{hint})", default=default, check=check)
        return "" if ans == "-" else ans

    # ------------------------------------------------------------ steps
    def sample(self, n: int, taken: set, prev: dict | None) -> dict:
        self.say(f"\nSample {n}")
        path = self.ask("  Reads (FASTQ/FASTA, .gz is fine) or a sorted, indexed BAM",
                        required=True)
        s: dict = {}
        s["bam" if path.lower().endswith(".bam") else "reads"] = path

        def name_ok(a):
            if not ID_RE.fullmatch(a):
                return f"Not a valid name: {NAME_HINT}."
            if a in taken:
                return f"Sample name {a!r} is already used."
            return None
        default = config.default_sample_name(path)
        s["sample"] = self.ask("  Sample name", default=default if ID_RE.fullmatch(default) else None,
                               required=True, check=name_ok)
        genome_default = (prev or {}).get("genome") or ("hg38" if self.on_phoenix else None)
        hint = f" (or a lab preset: {phoenix.preset_names()})" if self.on_phoenix else ""
        s["genome"] = self.ask(f"  Genome FASTA{hint}", default=genome_default, required=True)
        if phoenix.is_genome_preset(s["genome"]):
            note = phoenix.GENOME_PRESETS[s["genome"]]["note"]
            gtf = self.ask(f"  Reference annotation GTF (Enter to use the preset's: {note})",
                           default=(prev or {}).get("gtf") or "")
        else:
            gtf = self.ask("  Reference annotation GTF", default=(prev or {}).get("gtf"),
                           required=True)
        if gtf:
            s["gtf"] = gtf
        plat = self.choose("  Platform", PLATFORMS, guess_platform(path))
        if plat:
            s["platform"] = plat
        lib = self.choose("  Library", LIBRARIES, guess_library(path))
        if lib:
            s["library"] = lib
        return {k: s[k] for k in config.ALL_SAMPLE_COLUMNS if k in s}

    def settings(self, run_id: str) -> str:
        self.say(f"  Settings for {run_id}:")
        for i, p in enumerate(presets.PRESETS, 1):
            args = f"  ({p['args']})" if p["args"] else ""
            self.say(f"    {i}) {p['label']}{args}")
        ans = self.ask("  Pick a number, or type your own flair transcriptome arguments",
                       default="1")
        if ans.isdigit() and 1 <= int(ans) <= len(presets.PRESETS):
            return presets.PRESETS[int(ans) - 1]["args"]
        return ans

    def runs(self) -> list[dict]:
        self.say("\nRuns: which FLAIR builds do you want to compare?")
        self.say("  Tokens: a release (3.0.1), dev, main, git:<branch>, git:<url>@<ref>,")
        self.say("  env:/path/to/conda/env or container:<image>. Separate them with commas.")

        def tokens_ok(a):
            try:
                parse_tokens(a)
            except TokenError as e:
                return str(e)
            return None
        value = self.ask("  FLAIR", default="3.0.1", check=tokens_ok)
        runs = []
        for t in parse_tokens(value):
            runs.append({"id": t.id, "flair": dict(t.flair),
                         "transcriptome_args": self.settings(t.id)})

        while self.yes("\nAdd another run (for example the same FLAIR with other settings)?"):
            taken = {r["id"] for r in runs}

            def token_ok(a):
                try:
                    parse_token(a)
                except TokenError as e:
                    return str(e)
                return None
            tok = parse_token(self.ask("  FLAIR", default=_token_of(runs[-1]),
                                       check=token_ok))
            args = self.settings(tok.id)
            key = next((p["key"] for p in presets.PRESETS if p["args"] == args), "custom")
            base = tok.id if tok.id not in taken else f"{tok.id}-{key}"
            default_id, n = base, 2
            while default_id in taken:
                default_id, n = f"{base}-{n}", n + 1

            def id_ok(a):
                if not ID_RE.fullmatch(a):
                    return f"Not a valid run id: {NAME_HINT}."
                return f"Run id {a!r} is already used." if a in taken else None
            rid = self.ask("  Run id", default=default_id, check=id_ok)
            run = {"id": rid, "flair": dict(tok.flair), "transcriptome_args": args}
            label = self.ask("  Label shown in the report (Enter to use the id)", default="")
            if label:
                run["label"] = label
            runs.append(run)
        return runs

    def outdir(self) -> str:
        self.say("\nOutput folder for this comparison: bundles, the report and (on phoenix) "
                 "Nextflow's work files.")
        self.say("  Use groups or warm-archive storage, not your home folder (small quota).")
        return self.ask("  Output folder (--outdir)", default=self.default_outdir, required=True)

    def run(self) -> tuple[list[dict], list[dict], str] | None:
        """Return (samples, runs, outdir), or None when the user picks the test data instead."""
        self.say("flair-eval init: set up a FLAIR comparison. Press Enter to accept [defaults].")
        if self.yes("First time? Try flair-eval on its bundled test data instead (no files needed)?"):
            return None
        samples: list[dict] = []
        while True:
            samples.append(self.sample(len(samples) + 1, {s["sample"] for s in samples},
                                       samples[-1] if samples else None))
            if not self.yes("Add another sample?"):
                break
        runs = self.runs()
        return samples, runs, self.outdir()


def _token_of(run: dict) -> str:
    from .tokens import token_for_flair
    return token_for_flair(run["flair"]) or "3.0.1"


def write_configs(directory: Path, samples: list[dict], runs: list[dict],
                  force: bool = False) -> tuple[Path, Path]:
    directory.mkdir(parents=True, exist_ok=True)
    sp, rp = directory / "samples.csv", directory / "runs.yaml"
    existing = [p for p in (sp, rp) if p.exists()]
    if existing and not force:
        raise FileExistsError(", ".join(str(p) for p in existing))
    sp.write_text(config.samples_csv(samples), encoding="utf-8")
    rp.write_text(config.runs_yaml(runs), encoding="utf-8")
    return sp, rp

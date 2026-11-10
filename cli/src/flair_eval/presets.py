"""Common `transcriptome_args` presets, shown by `flair-eval init` and the builder page.

Flags are from `flair transcriptome --help` (FLAIR 3.0.1). Keep this list
in step with PRESETS in frontend/builder/logic.js (a test checks they match).
"""

PRESETS = [
    {"key": "default", "label": "FLAIR defaults", "args": "",
     "note": "No extra arguments."},
    {"key": "stringent", "label": "Stringent", "args": "--stringent",
     "note": "Supporting reads must be full-length (span 25 bp of the first and last exons)."},
    {"key": "check_splice", "label": "Check splice sites", "args": "--check_splice",
     "note": "Reads must cover 4 of 6 bp around each splice site, with no insertion over 3 bp."},
    {"key": "stringent_splice", "label": "Stringent + check splice sites",
     "args": "--stringent --check_splice", "note": "Both of the above."},
    {"key": "support5", "label": "More read support (5)", "args": "--support 5",
     "note": "At least 5 supporting reads per isoform (default 3)."},
    {"key": "noaligntoannot", "label": "Skip annotation alignment", "args": "--noaligntoannot",
     "note": "Detect transcripts from the genome alignment only. Faster, less accurate."},
]


def by_key(key: str) -> dict | None:
    return next((p for p in PRESETS if p["key"] == key), None)

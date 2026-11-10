"""flair-eval: a friendly front end for the flair-eval Nextflow pipeline.

The pipeline compares FLAIR runs (versions, branches, settings) on the same data and
writes one self-contained HTML report. This package writes and checks its config files
(samples.csv, runs.yaml) and builds the `nextflow run` command line.
Interfaces follow docs/CONTRACTS.md in the pipeline repository.
"""

__version__ = "0.1.0"

#: The pipeline as Nextflow pulls it.
PIPELINE = "hdheath/flair-eval"
#: Pipeline revision (`-r`) used when none is given. v2 is the integration branch until
#: the first v2 release tag exists.
DEFAULT_REF = "v2"

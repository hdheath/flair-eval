# FLAIR Evaluation Pipeline

A Nextflow pipeline for evaluating FLAIR isoform predictions using comprehensive metrics including splice junction analysis, transcript classification, and transcript end distance (TED) evaluation.

## Quick Start

```bash
nextflow run flair.test.suite.nf --input samples.csv --test_name "SMARCA4_locus_comprehensive"
```

## Usage

### Required Parameters

- `--input`: Path to CSV samplesheet containing sample information (required)
- `--test_name`: Name for the test run (default: `flair_test_suite`)

### Samplesheet Format

The input samplesheet must be a CSV file with the following columns:

| Column      | Required | Description                                    |
|-------------|----------|------------------------------------------------|
| sample_id   | Yes      | Unique sample identifier                       |
| genome      | Yes      | Path to genome FASTA file                      |
| gtf         | Yes      | Path to GTF annotation file                    |
| bam         | Yes      | Path to aligned BAM file                       |
| reads       | No       | Path to raw reads file (for FlairAlign)        |
| cage        | No       | Path to CAGE peaks file                        |
| drna    | No       | Path to dRNA peaks file                    |
| junction_tab| No       | STAR SJ.out.tab (or use legacy `junctions`)   |
| cage_signal_plus   | No | Optional CAGE plus-strand bedGraph signal      |
| cage_signal_minus  | No | Optional CAGE minus-strand bedGraph signal     |
| drna_signal_plus  | No | Optional dRNA plus-strand bedGraph signal |
| drna_signal_minus | No | Optional dRNA minus-strand bedGraph signal |

**Example samplesheet (`samples.csv`):**

```csv
sample_id,genome,gtf,bam,reads,cage,drna,junction_tab,cage_signal_plus,cage_signal_minus,drna_signal_plus,drna_signal_minus
B1A_kd_induced_rep1,/path/to/genome.fa,/path/to/annotation.gtf,/path/to/sample1.bam,,/path/to/cage.bed,/path/to/drna.bed,/path/to/SJ.out.tab,/path/to/cage_plus.bg,/path/to/cage_minus.bg,/path/to/quant_plus.bg,/path/to/quant_minus.bg
B1A_kd_induced_rep2,/path/to/genome.fa,/path/to/annotation.gtf,/path/to/sample2.bam,,/path/to/cage.bed,/path/to/drna.bed,,,,,
```

Note: `junctions` is still accepted as a legacy alias of `junction_tab`.

### Running the Pipeline

**Basic usage with custom test name:**

```bash
nextflow run flair.test.suite.nf \
    --input samples.csv \
    --test_name "my_experiment"
```

**Using default test name (will show warning):**

```bash
nextflow run flair.test.suite.nf --input samples.csv
```

This will use the default test name `flair_test_suite` and display a warning.

## Output

Results are organized in the `results/` directory:

```
results/
├── align/                     # Alignment outputs (if reads provided)
├── partition/                 # Partitioned regions
├── transcriptome/             # FLAIR transcriptome outputs
├── evaluations/
│   ├── individual/           # Per-sample evaluation metrics
│   └── summary/              # Combined evaluation summary
└── isoform_plots/            # Visualization outputs
```

### Key Output Files

- `results/evaluations/summary/<test_name>_evaluation_summary.tsv` - Comprehensive evaluation metrics for all samples
- Individual TED and FLAIR evaluation files in `results/evaluations/individual/`

## Configuration

Edit `nextflow.config` to modify:

- Executor settings (SLURM configuration)
- Resource allocations (CPU, memory, time)
- Conda environment paths

## Pipeline Modes

The pipeline supports different operational modes configured in the workflow:

- **Alignment modes**: Default alignment settings
- **Partition modes**: Region-specific partitioning (e.g., SMARCA4 locus: chr19:10900001-11100000)
- **Transcriptome modes**: With or without GTF guidance

## Requirements

- Nextflow >= 20.0
- Conda/Mamba
- FLAIR
- bedtools
- samtools
- Python 3.7+

## Example Run

See `samples.csv` for a complete example with 8 SMARCA4 samples.

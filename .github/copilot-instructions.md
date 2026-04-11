# FLAIR-Eval Pipeline

Nextflow pipeline for benchmarking long-read isoform assemblers (FLAIR, Bambu, IsoQuant, etc.) using orthogonal signal data (CAGE, QuantSeq).

## Architecture

- `main.nf` — Pipeline entry point (DSL2)
- `modules/` — Nextflow processes grouped by stage (align, partition, evaluation, etc.)
- `subworkflows/` — Composed workflows combining processes
- `bin/evaluation/` — Python package (~58 modules) for metrics, plotting, analysis
- `bin/evaluation/signal_utils.py` — **Canonical shared parser** for BED12/GTF/signal data
- `bin/evaluation/pub_style.py` — Publication-quality plot styling (Nature Portfolio specs)
- `params/` — Per-dataset YAML config files
- `tests/` — pytest test suite

## Code Style

### Python (`bin/evaluation/`)

- **Shared parsing**: Always import from `signal_utils.py` — never write inline BED12/GTF parsers
  - `parse_bed12(path)` → `List[dict]` with keys: chrom, start, end, name, score, strand, tss, tts, junctions, n_exons, spliced_len
  - `parse_bed12_by_name(path)` → `Dict[str, dict]`
  - `group_by_junction_chain(isoforms)` → `Dict[tuple, List[dict]]`
  - `tss_tts(start, end, strand)` → `(tss, tts)`
  - Junctions are 2-tuples: `(donor, acceptor)` — NOT 4-tuples
- **Type annotations**: Required on all function signatures. Use `str | Path` (3.10+ union syntax), not `Union`
- **Naming**: `snake_case` for functions/variables, prefix private functions with `_`
- **Imports**: Relative imports within package (`from .utils import get_logger`). Standalone scripts use try/except pattern for flat vs package imports
- **Optional deps**: Wrap matplotlib/numpy in try/except with `HAS_MATPLOTLIB` flag
- **Logging**: `logger = get_logger()` per module
- **Timing**: Use `with timed_section("label"):` for performance-critical sections

### Plotting

- Import from `pub_style.py`: `style_ax`, `savefig`, `W1` (89mm), `W2` (183mm), `GOLDEN_RATIO` (1.618), `PALETTE`, `MODE_COLORS`
- Okabe-Ito colorblind-safe palette; 8pt labels, 7pt ticks, 300 DPI PNG + SVG
- Standalone plot scripts use `matplotlib.use('Agg')` before importing pyplot
- CLI pattern: `--bed label:path [label:path ...]`, `--output dir/`

### Nextflow

- Process names: CamelCase (`FlairAlign`, `Evaluation`, `GeneVariationPlot`)
- Always include `errorStrategy 'ignore'` for visualization processes
- Preserve metadata tuples through channels: `(test_name, dataset_name, ...)`
- Placeholder files: `NO_ISOFORMS_BED`, `NO_CAGE` — check with `.name != 'NO_FILE'`
- Scripts called as: `python ${projectDir}/bin/evaluation/script.py`

## Build & Run

```bash
# Activate Nextflow environment
conda activate nextflow_env

# Run pipeline
nextflow run main.nf -params-file params/config.yaml -resume

# Run Python analysis (standalone)
conda activate nextflow_env
python bin/evaluation/script.py --bed label:file.bed --output outdir/

# Run tests
conda activate nextflow_env
pytest tests/
```

## Conventions

- Results go to `results/<test_name>/` with subfolders: transcriptome/, evaluations/, summary/
- Evaluation TSVs are the primary metrics output; plots are secondary
- BED12 isoform files carry gene info in the name field: `ENST..._ENSG...`
- Gene extraction: `gene_from_name(iso_name)` parses ENSG IDs

## Key References

- See `docs/` for pipeline documentation
- See `bin/evaluation/signal_utils.py` docstring for full shared API
- See `bin/evaluation/pub_style.py` for plot styling constants

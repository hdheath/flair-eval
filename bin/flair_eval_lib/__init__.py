"""Shared Python for flair-eval v2 (structure metrics + report).

Stdlib only, so it runs in any Python 3.10+ (the report stream may use pandas/numpy,
but nothing here needs them).
"""

__version__ = "2.0.0.dev0"

BUNDLE_SCHEMA_ID = "flair-eval/bundle/1"

CATEGORIES = ("FSM", "ISM", "NIC", "NNC", "genic", "genic_intron",
              "antisense", "intergenic", "fusion")
FSM_SUBCATEGORIES = ("CIM", "FSM5", "FSM3", "FSMTE")

# Fixed histogram edges (CONTRACTS.md section 3). The last bin is open-ended.
HISTOGRAM_EDGES = {
    "isoforms_per_gene": [1, 2, 3, 4, 5, 6, 8, 11, 16, 21],
    "transcript_length": [0, 250, 500, 750, 1000, 1500, 2000, 3000, 4000, 6000, 10000],
    "exon_length": [0, 50, 100, 150, 200, 300, 500, 1000, 2000, 5000],
    "exons_per_isoform": [1, 2, 3, 4, 5, 6, 8, 11, 16, 21, 31],
}

import importlib.util
from pathlib import Path


_SCRIPT = Path(__file__).parents[1] / "bin" / "evaluation" / "peak_read_coverage.py"
_SPEC = importlib.util.spec_from_file_location("peak_read_coverage", _SCRIPT)
_MODULE = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(_MODULE)


def test_count_near_uses_half_open_right_boundary():
    # BED peak [1000, 1001) with W=50 supports [950, 1051): both edge
    # coordinates 950 and 1050 count, while 1051 is one base too far right.
    ends = [949, 950, 1000, 1050, 1051]
    assert _MODULE.count_near(ends, 950, 1051) == 3

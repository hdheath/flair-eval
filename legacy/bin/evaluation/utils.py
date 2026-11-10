"""
Shared utility functions for the evaluation package.

Provides timing, subprocess execution, logging, and common helper functions.
"""

import logging
import subprocess
import time
from collections import defaultdict
from contextlib import contextmanager
from pathlib import Path
from shutil import which as shutil_which
from typing import Optional

# Performance timing storage (module-level)
_TIMING_DATA = defaultdict(list)

# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format='[%(levelname)s] %(message)s'
)
logger = logging.getLogger(__name__)


def get_logger():
    """Get the module logger."""
    return logger


def get_timing_data():
    """Get the timing data dictionary for report generation."""
    return _TIMING_DATA


@contextmanager
def timed_section(section_name: str):
    """Context manager to time code sections and store results."""
    start = time.time()
    try:
        yield
    finally:
        elapsed = time.time() - start
        _TIMING_DATA[section_name].append(elapsed)
        logger.debug(f"[TIMING] {section_name}: {elapsed:.2f}s")


def which(exe: str) -> bool:
    """Check if executable exists on PATH."""
    return shutil_which(exe) is not None


def run(cmd: list) -> subprocess.CompletedProcess:
    """Run command and return result."""
    logger.debug(f"RUN: {' '.join(cmd)}")
    try:
        return subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, check=True)
    except subprocess.CalledProcessError as e:
        logger.error(f"Command failed ({' '.join(cmd)}): {e.stderr.strip()}")
        raise


def count_lines(path: Path) -> int:
    """Count lines in a file."""
    if not path.exists():
        return 0
    try:
        with open(path) as f:
            return sum(1 for _ in f)
    except Exception:
        return 0


def safe_f1(p: Optional[float], r: Optional[float]) -> Optional[float]:
    """Calculate F1 score from precision and recall."""
    if p is None or r is None:
        return None
    s = p + r
    return (2 * p * r / s) if s > 0 else 0.0


def write_timing_report(output_path: Path) -> None:
    """Write a timing report to the specified file."""
    with open(output_path, 'w') as tf:
        tf.write("TED Performance Timing Report\n")
        tf.write("=" * 80 + "\n\n")

        # Sort by total time (sum of all calls)
        sorted_sections = sorted(
            _TIMING_DATA.items(),
            key=lambda x: sum(x[1]),
            reverse=True
        )

        total_time = sum(sum(times) for times in _TIMING_DATA.values())
        tf.write(f"Total instrumented time: {total_time:.2f}s\n\n")

        tf.write(f"{'Section':<50} {'Calls':>8} {'Total (s)':>12} {'Avg (s)':>12} {'% of Total':>12}\n")
        tf.write("-" * 100 + "\n")

        for section, times in sorted_sections:
            n_calls = len(times)
            total = sum(times)
            avg = total / n_calls
            pct = (total / total_time * 100) if total_time > 0 else 0
            tf.write(f"{section:<50} {n_calls:>8} {total:>12.2f} {avg:>12.2f} {pct:>11.1f}%\n")

    logger.info(f"Wrote timing report to {output_path}")

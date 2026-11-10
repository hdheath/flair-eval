from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[2]


@pytest.fixture(scope="session")
def repo() -> Path:
    return REPO


@pytest.fixture
def sample() -> dict:
    return {"sample": "HG002", "reads": "/data/HG002.fq.gz", "genome": "/ref/hg38.fa",
            "gtf": "/ref/gencode.gtf"}


@pytest.fixture
def two_runs() -> list:
    return [
        {"id": "flair-3.0.1", "flair": {"release": "3.0.1"}, "transcriptome_args": ""},
        {"id": "flair-dev", "flair": {"git": "https://github.com/BrooksLabUCSC/flair",
                                       "ref": "dev"}, "transcriptome_args": ""},
    ]

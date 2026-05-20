"""
Shared pytest fixtures and sys.path wiring for the strain-deformation test suite.
"""
import sys
from pathlib import Path

import numpy as np
import pytest

# Make the project root importable without installing a package
ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(ROOT))


# ---------------------------------------------------------------------------
# Reusable numpy helpers
# ---------------------------------------------------------------------------

def assert_matrix_close(a: np.ndarray, b: np.ndarray, rtol=1e-8, atol=1e-10) -> None:
    """Thin wrapper around np.testing.assert_allclose with descriptive labels."""
    np.testing.assert_allclose(a, b, rtol=rtol, atol=atol,
                               err_msg=f"\nGot:\n{a}\nExpected:\n{b}")


# ---------------------------------------------------------------------------
# POSCAR / CONTCAR fixture
# ---------------------------------------------------------------------------

BCC_FE_CONTCAR = """\
Fe BCC
1.0
  2.8700000000000000   0.0000000000000000   0.0000000000000000
  0.0000000000000000   2.8700000000000000   0.0000000000000000
  0.0000000000000000   0.0000000000000000   2.8700000000000000
Fe
1
Direct
  0.0000000000000000   0.0000000000000000   0.0000000000000000
"""


@pytest.fixture
def bcc_fe_contcar(tmp_path: Path) -> Path:
    """Write a minimal BCC Fe CONTCAR to a temp directory and return its path."""
    p = tmp_path / "CONTCAR"
    p.write_text(BCC_FE_CONTCAR)
    return p


@pytest.fixture
def bcc_fe_axis() -> np.ndarray:
    """3×3 lattice matrix for BCC Fe (a = 2.87 Å)."""
    return np.eye(3) * 2.87

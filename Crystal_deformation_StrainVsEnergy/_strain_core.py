"""
Shared numerical kernels for Lagrangian-strain crystal deformation.

All functions operate on plain NumPy arrays and have no I/O side-effects,
making them straightforward to unit-test.
"""
from __future__ import annotations

import numpy as np


def build_voigt_strains(dc: str, eta: float) -> np.ndarray:
    """
    Map a 6-character deformation code to Voigt strain components.

    Characters
    ----------
    'E'  → +eta   (positive strain)
    'e'  → -eta   (negative strain; used for volume-conserving shear partners)
    '0'  → 0
    '2'  → 2*eta  (doubled off-diagonal, so e[i]/2 = eta in the tensor)

    Parameters
    ----------
    dc : str
        6-character code, e.g. 'EEE000', 'Ee0000', '00E002'.
    eta : float
        Strain amplitude.

    Returns
    -------
    np.ndarray, shape (6,)
        Voigt components [e1, e2, e3, e4, e5, e6].
    """
    if len(dc) != 6:
        raise ValueError(f"Deformation code must be 6 characters, got '{dc}'")
    mapping: dict[str, float] = {'E': eta, 'e': -eta, '0': 0.0, '2': 2.0 * eta}
    try:
        return np.array([mapping[ch] for ch in dc], dtype=float)
    except KeyError as exc:
        raise ValueError(f"Unknown character {exc} in deformation code '{dc}'") from exc


def voigt_to_eta_matrix(e: np.ndarray) -> np.ndarray:
    """
    Assemble the symmetric 3×3 Lagrangian strain tensor η from Voigt components.

        η = | e[0]    e[5]/2  e[4]/2 |
            | e[5]/2  e[1]    e[3]/2 |
            | e[4]/2  e[3]/2  e[2]   |

    Parameters
    ----------
    e : np.ndarray, shape (6,)

    Returns
    -------
    np.ndarray, shape (3, 3), symmetric.
    """
    return np.array([
        [e[0],    e[5]/2, e[4]/2],
        [e[5]/2,  e[1],   e[3]/2],
        [e[4]/2,  e[3]/2, e[2]  ],
    ], dtype=float)


def lagrangian_to_linear(eta: np.ndarray, tol: float = 1e-10) -> np.ndarray:
    """
    Recover the linear (infinitesimal) strain ε from the Lagrangian strain η.

    Solves the implicit relation  η = ε + ½ ε²  iteratively:

        x_{n+1} = η − ½ εₙ²,   until  ‖x − ε‖ < tol.

    Parameters
    ----------
    eta : np.ndarray, shape (3, 3)
        Lagrangian strain tensor.
    tol : float
        Convergence criterion on the Frobenius-norm increment.

    Returns
    -------
    np.ndarray, shape (3, 3)
        Linear strain tensor ε such that η = ε + ½ ε·ε  (to tolerance).

    Raises
    ------
    ValueError  If ‖η‖ > 0.7 (deformation too large for the scheme).
    RuntimeError  If the iteration does not converge within 1000 steps.
    """
    if np.linalg.norm(eta) > 0.7:
        raise ValueError(
            f"Strain norm {np.linalg.norm(eta):.4f} > 0.7; reduce deformation amplitude."
        )
    eps = eta.copy()
    for _ in range(1000):
        x = eta - 0.5 * eps @ eps
        if np.linalg.norm(x - eps) < tol:
            return x
        eps = x
    raise RuntimeError("Lagrangian→linear strain iteration did not converge in 1000 steps.")


def deform_cell(axis_matrix: np.ndarray, def_matrix: np.ndarray) -> np.ndarray:
    """
    Apply a deformation matrix D to unit-cell basis vectors.

    Each row of *axis_matrix* is a lattice vector.  The transformation is:

        R'ᵀ = D · Rᵀ   →   R' = (D · Rᵀ)ᵀ

    Parameters
    ----------
    axis_matrix : np.ndarray, shape (3, 3)
        Rows are lattice vectors a₁, a₂, a₃.
    def_matrix : np.ndarray, shape (3, 3)
        Deformation matrix D = I + ε.

    Returns
    -------
    np.ndarray, shape (3, 3)
        Deformed lattice vectors (rows).
    """
    return (def_matrix @ axis_matrix.T).T

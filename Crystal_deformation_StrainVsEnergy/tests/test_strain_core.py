"""
Unit tests for _strain_core.py — the shared numerical kernels.

Every test is deterministic and requires only NumPy; no VASP files needed.
"""
import numpy as np
import pytest
from conftest import assert_matrix_close

from _strain_core import (
    build_voigt_strains,
    deform_cell,
    lagrangian_to_linear,
    voigt_to_eta_matrix,
)


# ===========================================================================
# build_voigt_strains
# ===========================================================================

class TestBuildVoigtStrains:

    def test_volumetric_EEE000(self):
        e = build_voigt_strains("EEE000", 0.05)
        np.testing.assert_allclose(e, [0.05, 0.05, 0.05, 0.0, 0.0, 0.0])

    def test_shear_partner_Ee0000(self):
        # 'e' gives -eta, confirming sign flip
        e = build_voigt_strains("Ee0000", 0.05)
        np.testing.assert_allclose(e, [0.05, -0.05, 0.0, 0.0, 0.0, 0.0])

    def test_doubled_offdiagonal_00E002(self):
        # '2' → 2*eta; used so that e[5]/2 = eta in the tensor
        e = build_voigt_strains("00E002", 0.05)
        np.testing.assert_allclose(e, [0.0, 0.0, 0.05, 0.0, 0.0, 0.10])

    def test_zero_strain(self):
        e = build_voigt_strains("EEEEEE", 0.0)
        np.testing.assert_allclose(e, np.zeros(6))

    def test_negative_eta(self):
        # 'E' → eta (negative), 'e' → -eta (positive)
        e = build_voigt_strains("Ee0000", -0.03)
        np.testing.assert_allclose(e, [-0.03, 0.03, 0.0, 0.0, 0.0, 0.0])

    def test_invalid_character_raises(self):
        with pytest.raises(ValueError, match="Unknown character"):
            build_voigt_strains("EEE00X", 0.05)

    def test_wrong_length_raises(self):
        with pytest.raises(ValueError, match="6 characters"):
            build_voigt_strains("EEE", 0.05)


# ===========================================================================
# voigt_to_eta_matrix
# ===========================================================================

class TestVoigtToEtaMatrix:

    def test_isotropic_gives_diagonal(self):
        e = np.array([0.05, 0.05, 0.05, 0.0, 0.0, 0.0])
        eta = voigt_to_eta_matrix(e)
        assert_matrix_close(eta, np.eye(3) * 0.05)

    def test_xy_shear_only(self):
        # e6 = 2*η → off-diag e[5]/2 = η
        e = np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.10])
        eta = voigt_to_eta_matrix(e)
        expected = np.array([
            [0.0, 0.05, 0.0],
            [0.05, 0.0, 0.0],
            [0.0,  0.0, 0.0],
        ])
        assert_matrix_close(eta, expected)

    def test_symmetry(self):
        e = np.array([0.02, 0.03, 0.04, 0.05, 0.06, 0.07])
        eta = voigt_to_eta_matrix(e)
        assert_matrix_close(eta, eta.T)

    def test_diagonal_mapping(self):
        e = np.array([0.01, 0.02, 0.03, 0.0, 0.0, 0.0])
        eta = voigt_to_eta_matrix(e)
        assert eta[0, 0] == pytest.approx(0.01)
        assert eta[1, 1] == pytest.approx(0.02)
        assert eta[2, 2] == pytest.approx(0.03)

    def test_voigt_index_e4_yz(self):
        e = np.array([0.0, 0.0, 0.0, 0.08, 0.0, 0.0])
        eta = voigt_to_eta_matrix(e)
        # e4 = e[3] → η[1,2] = η[2,1] = e[3]/2
        assert eta[1, 2] == pytest.approx(0.04)
        assert eta[2, 1] == pytest.approx(0.04)


# ===========================================================================
# lagrangian_to_linear
# ===========================================================================

class TestLagrangianToLinear:

    def test_zero_strain(self):
        eta = np.zeros((3, 3))
        eps = lagrangian_to_linear(eta)
        assert_matrix_close(eps, np.zeros((3, 3)))

    def test_small_strain_identity_limit(self):
        # For tiny η, ε ≈ η (higher-order corrections negligible)
        eta = 1e-4 * np.eye(3)
        eps = lagrangian_to_linear(eta)
        assert_matrix_close(eps, eta, rtol=1e-4)

    def test_reconstruction_roundtrip(self):
        # Must satisfy η = ε + ½ ε²  to machine precision
        eta_in = np.diag([0.05, -0.03, 0.02])
        eps    = lagrangian_to_linear(eta_in)
        eta_re = eps + 0.5 * eps @ eps
        assert_matrix_close(eta_re, eta_in, atol=1e-9)

    def test_shear_reconstruction(self):
        e   = build_voigt_strains("000EEE", 0.04)
        eta = voigt_to_eta_matrix(e)
        eps = lagrangian_to_linear(eta)
        assert_matrix_close(eps + 0.5 * eps @ eps, eta, atol=1e-9)

    def test_too_large_strain_raises(self):
        with pytest.raises(ValueError, match="0.7"):
            lagrangian_to_linear(0.8 * np.eye(3))

    def test_result_is_symmetric_for_symmetric_input(self):
        eta = np.diag([0.04, -0.02, 0.01])
        eta += eta.T
        eta /= 2  # ensure symmetric
        eps = lagrangian_to_linear(eta)
        assert_matrix_close(eps, eps.T, atol=1e-12)


# ===========================================================================
# deform_cell
# ===========================================================================

class TestDeformCell:

    def test_identity_deformation_unchanged(self):
        axis = np.array([[3.0, 0.0, 0.0],
                         [0.0, 4.0, 0.0],
                         [0.0, 0.0, 5.0]])
        new  = deform_cell(axis, np.eye(3))
        assert_matrix_close(new, axis)

    def test_uniform_dilation(self):
        axis = np.eye(3) * 3.0
        D    = 1.1 * np.eye(3)
        new  = deform_cell(axis, D)
        assert_matrix_close(new, np.eye(3) * 3.3)

    def test_volume_scales_with_det_D(self):
        axis = np.eye(3) * 2.87   # BCC Fe
        e    = build_voigt_strains("EEE000", 0.05)
        eta  = voigt_to_eta_matrix(e)
        eps  = lagrangian_to_linear(eta)
        D    = np.eye(3) + eps
        new  = deform_cell(axis, D)
        ratio = abs(np.linalg.det(new)) / abs(np.linalg.det(axis))
        assert ratio == pytest.approx(abs(np.linalg.det(D)), rel=1e-8)

    def test_shear_preserves_first_vector(self):
        # Pure xy shear: only row-1 and row-0 interact
        axis = np.eye(3)
        D    = np.eye(3)
        D[0, 1] = 0.1   # shear in x due to y
        new  = deform_cell(axis, D)
        # Row 0 of axis is [1,0,0]; D·[1,0,0]ᵀ = [1, 0, 0]ᵀ
        assert_matrix_close(new[0], [1.0, 0.0, 0.0])

    def test_volumetric_strain_increases_volume(self):
        axis = np.eye(3) * 2.87
        eta  = 0.05
        e    = build_voigt_strains("EEE000", eta)
        eps  = lagrangian_to_linear(voigt_to_eta_matrix(e))
        new  = deform_cell(axis, np.eye(3) + eps)
        V0   = np.linalg.det(axis)
        V1   = np.linalg.det(new)
        assert V1 > V0

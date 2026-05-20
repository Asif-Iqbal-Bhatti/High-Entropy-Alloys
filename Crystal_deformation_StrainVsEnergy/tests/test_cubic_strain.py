"""
Tests for cubic_strain.py — cubic-specific deformation matrices.
"""
import numpy as np
import pytest
from conftest import assert_matrix_close

from cubic_strain import DEFORMATION_CODES, build_eta_matrix_cubic
from _strain_core import build_voigt_strains, lagrangian_to_linear


# ===========================================================================
# DEFORMATION_CODES (cubic)
# ===========================================================================

class TestCubicDeformationCodes:

    def test_only_three_codes(self):
        assert set(DEFORMATION_CODES.keys()) == {0, 1, 2}

    def test_code_0_is_EEE000(self):
        dc, _ = DEFORMATION_CODES[0]
        assert dc == "EEE000"

    def test_code_1_is_EeE000(self):
        dc, _ = DEFORMATION_CODES[1]
        assert dc == "EeE000"

    def test_code_2_is_00E002(self):
        dc, _ = DEFORMATION_CODES[2]
        assert dc == "00E002"


# ===========================================================================
# build_eta_matrix_cubic
# ===========================================================================

class TestBuildEtaMatrixCubic:

    # ---- Code 0: volumetric ------------------------------------------------

    def test_code0_diagonal(self):
        e   = build_voigt_strains("EEE000", 0.05)
        eta = build_eta_matrix_cubic(e, 0)
        assert_matrix_close(eta, np.eye(3) * 0.05)

    def test_code0_no_volume_conservation_constraint(self):
        # Code 0 does NOT enforce det(I+η)=1; volume changes by (1+η)^3
        e   = build_voigt_strains("EEE000", 0.05)
        eta = build_eta_matrix_cubic(e, 0)
        det = np.linalg.det(np.eye(3) + eta)
        assert det == pytest.approx((1.05)**3, rel=1e-8)

    # ---- Code 1: orthorhombic volume-conserving ----------------------------

    def test_code1_volume_conserving_to_second_order(self):
        # det(I + η) should be 1 to O(η²) for any small η
        for eta_val in [0.01, 0.03, 0.05, 0.08]:
            e   = build_voigt_strains("EeE000", eta_val)
            eta = build_eta_matrix_cubic(e, 1)
            det = np.linalg.det(np.eye(3) + eta)
            assert det == pytest.approx(1.0, abs=2e-4), \
                f"Volume not conserved at η={eta_val}: det={det:.8f}"

    def test_code1_zz_element_formula(self):
        eta_val = 0.05
        e       = build_voigt_strains("EeE000", eta_val)
        eta_mat = build_eta_matrix_cubic(e, 1)
        expected_zz = 1.0 / (1.0 - eta_val**2) - 1.0
        assert eta_mat[2, 2] == pytest.approx(expected_zz, rel=1e-10)

    def test_code1_xy_off_diagonals_zero(self):
        e   = build_voigt_strains("EeE000", 0.05)
        eta = build_eta_matrix_cubic(e, 1)
        assert eta[0, 1] == pytest.approx(0.0)
        assert eta[1, 0] == pytest.approx(0.0)

    # ---- Code 2: monoclinic shear volume-conserving ------------------------

    def test_code2_volume_conserving_to_second_order(self):
        for eta_val in [0.01, 0.03, 0.05, 0.08]:
            e   = build_voigt_strains("00E002", eta_val)
            eta = build_eta_matrix_cubic(e, 2)
            det = np.linalg.det(np.eye(3) + eta)
            assert det == pytest.approx(1.0, abs=2e-4), \
                f"Volume not conserved at η={eta_val}: det={det:.8f}"

    def test_code2_has_xy_shear(self):
        # dc='00E002': e[5] = 2*eta → η[0,1] = η[1,0] = eta
        eta_val = 0.05
        e       = build_voigt_strains("00E002", eta_val)
        eta     = build_eta_matrix_cubic(e, 2)
        assert eta[0, 1] == pytest.approx(eta_val, rel=1e-10)
        assert eta[1, 0] == pytest.approx(eta_val, rel=1e-10)

    def test_code2_diagonal_xx_yy_zero(self):
        e   = build_voigt_strains("00E002", 0.05)
        eta = build_eta_matrix_cubic(e, 2)
        assert eta[0, 0] == pytest.approx(0.0)
        assert eta[1, 1] == pytest.approx(0.0)

    def test_code2_zz_element_formula(self):
        eta_val = 0.05
        e       = build_voigt_strains("00E002", eta_val)
        eta_mat = build_eta_matrix_cubic(e, 2)
        expected_zz = 1.0 / (1.0 - eta_val**2) - 1.0
        assert eta_mat[2, 2] == pytest.approx(expected_zz, rel=1e-10)

    # ---- Codes 1 and 2 differ (different dc strings) -----------------------

    def test_code1_and_code2_matrices_differ(self):
        e1  = build_voigt_strains("EeE000", 0.05)
        e2  = build_voigt_strains("00E002", 0.05)
        m1  = build_eta_matrix_cubic(e1, 1)
        m2  = build_eta_matrix_cubic(e2, 2)
        assert not np.allclose(m1, m2), "Codes 1 and 2 should produce different matrices"

    # ---- Roundtrip: η = ε + ½ε² -----------------------------------------

    @pytest.mark.parametrize("code,dc,eta_val", [
        (0, "EEE000", 0.05),
        (1, "EeE000", 0.05),
        (2, "00E002", 0.05),
        (0, "EEE000", 0.01),
        (1, "EeE000", 0.08),
    ])
    def test_lagrangian_roundtrip(self, code, dc, eta_val):
        e   = build_voigt_strains(dc, eta_val)
        eta = build_eta_matrix_cubic(e, code)
        eps = lagrangian_to_linear(eta)
        eta_re = eps + 0.5 * eps @ eps
        assert_matrix_close(eta_re, eta, atol=1e-9)

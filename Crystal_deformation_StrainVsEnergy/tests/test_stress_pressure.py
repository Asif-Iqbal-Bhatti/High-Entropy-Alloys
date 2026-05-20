"""
Tests for stress_pressure.py — hydrostatic pressure and deviatoric stress.

All tests use analytic inputs; no vasprun.xml file required.
"""
import numpy as np
import pytest
from conftest import assert_matrix_close

from stress_pressure import (
    deviatoric_stress_matrix,
    deviatoric_stress_voigt,
    hydrostatic_pressure,
    manual_kbar_example,
)


# ===========================================================================
# hydrostatic_pressure
# ===========================================================================

class TestHydrostaticPressure:

    def test_isotropic_compression(self):
        # σ = σ₀ I with σ₀ < 0 → P = -Tr(σ)/3 = -σ₀ > 0
        sigma0 = -2.0  # eV/Å³, compressive
        stress = sigma0 * np.eye(3)
        p      = hydrostatic_pressure(stress)
        assert p == pytest.approx(-sigma0, rel=1e-12)

    def test_isotropic_tension(self):
        stress = 3.0 * np.eye(3)
        p      = hydrostatic_pressure(stress)
        assert p == pytest.approx(-3.0, rel=1e-12)

    def test_zero_stress(self):
        assert hydrostatic_pressure(np.zeros((3, 3))) == pytest.approx(0.0)

    def test_general_asymmetric_stress(self):
        # Hand-computed: Tr = 1 + 5 + 9 = 15 → P = -5
        stress = np.array([
            [1.0, 2.0, 3.0],
            [4.0, 5.0, 6.0],
            [7.0, 8.0, 9.0],
        ], dtype=float)
        assert hydrostatic_pressure(stress) == pytest.approx(-5.0, rel=1e-12)

    def test_units_scaling(self):
        # Doubling the stress should double the pressure
        stress = np.diag([1.0, 2.0, 3.0])
        p1 = hydrostatic_pressure(stress)
        p2 = hydrostatic_pressure(2.0 * stress)
        assert p2 == pytest.approx(2.0 * p1, rel=1e-12)


# ===========================================================================
# deviatoric_stress_matrix
# ===========================================================================

class TestDeviatoricStressMatrix:

    def test_isotropic_gives_zero_deviatoric(self):
        stress = -2.0 * np.eye(3)
        dev    = deviatoric_stress_matrix(stress)
        assert_matrix_close(dev, np.zeros((3, 3)), atol=1e-14)

    def test_traceless(self):
        stress = np.array([
            [1.0, 0.5, 0.0],
            [0.5, 2.0, 0.3],
            [0.0, 0.3, 3.0],
        ])
        dev = deviatoric_stress_matrix(stress)
        assert np.trace(dev) == pytest.approx(0.0, abs=1e-12)

    def test_decomposition_recovers_original(self):
        # σ = s + P·I, so s + (P·I) == σ
        stress = np.array([
            [-0.1,  0.02, 0.01],
            [ 0.02, -0.3, 0.005],
            [ 0.01, 0.005, -0.2],
        ])
        p   = hydrostatic_pressure(stress)
        dev = deviatoric_stress_matrix(stress)
        reconstructed = dev - p * np.eye(3)
        assert_matrix_close(reconstructed, stress, atol=1e-14)

    def test_pure_shear_unchanged(self):
        # Off-diagonal shear has zero trace → deviatoric = original
        stress = np.array([
            [0.0, 0.05, 0.0],
            [0.05, 0.0,  0.0],
            [0.0,  0.0,  0.0],
        ])
        dev = deviatoric_stress_matrix(stress)
        assert_matrix_close(dev, stress, atol=1e-14)

    def test_symmetry_preserved(self):
        stress = np.array([
            [1.0, 0.3, 0.1],
            [0.3, 2.0, 0.2],
            [0.1, 0.2, 3.0],
        ])
        dev = deviatoric_stress_matrix(stress)
        assert_matrix_close(dev, dev.T, atol=1e-14)


# ===========================================================================
# deviatoric_stress_voigt
# ===========================================================================

class TestDeviatoricStressVoigt:

    def test_isotropic_diagonal_zero(self):
        # [σ, σ, σ, 0, 0, 0] → dev[:3] = 0
        sv  = np.array([-2.0, -2.0, -2.0, 0.0, 0.0, 0.0])
        dev, p = deviatoric_stress_voigt(sv)
        np.testing.assert_allclose(dev[:3], 0.0, atol=1e-14)

    def test_pressure_sign_convention(self):
        # Compressive diagonal (tension-positive convention): P > 0
        sv  = np.array([-1.0, -2.0, -3.0, 0.0, 0.0, 0.0])
        _, p = deviatoric_stress_voigt(sv)
        assert p == pytest.approx(2.0, rel=1e-12)  # p = -mean([-1,-2,-3]) = 2

    def test_off_diagonal_unchanged(self):
        sv      = np.array([1.0, 1.0, 1.0, 0.05, 0.03, 0.02])
        dev, _  = deviatoric_stress_voigt(sv)
        np.testing.assert_allclose(dev[3:], [0.05, 0.03, 0.02])

    def test_diagonal_only_modified(self):
        # Only dev[:3] should differ from input[:3]
        sv     = np.array([-0.1, -0.2, -0.3, 0.01, 0.02, 0.03])
        dev, p = deviatoric_stress_voigt(sv)
        # dev[i] = sv[i] + p for i in 0,1,2
        expected_diag = sv[:3] + p
        np.testing.assert_allclose(dev[:3], expected_diag, atol=1e-14)

    def test_deviatoric_diagonal_sum_zero(self):
        sv     = np.array([-0.1, -0.2, -0.3, 0.05, 0.04, 0.03])
        dev, _ = deviatoric_stress_voigt(sv)
        assert np.sum(dev[:3]) == pytest.approx(0.0, abs=1e-12)

    def test_consistency_with_matrix_version(self):
        # The two implementations must agree on pressure
        sv = np.array([-0.05, -0.08, -0.07, 0.01, 0.02, 0.005])
        # Build 3×3 from Voigt (symmetric)
        stress_3x3 = np.array([
            [sv[0], sv[5], sv[4]],
            [sv[5], sv[1], sv[3]],
            [sv[4], sv[3], sv[2]],
        ])
        p_matrix = hydrostatic_pressure(stress_3x3)
        _, p_voigt = deviatoric_stress_voigt(sv)
        assert p_voigt == pytest.approx(p_matrix, rel=1e-10)


# ===========================================================================
# manual_kbar_example (smoke test — just checks it runs without exceptions)
# ===========================================================================

def test_manual_kbar_example_runs(capsys):
    manual_kbar_example()
    out = capsys.readouterr().out
    assert "pressure" in out.lower()
    assert "GPa" in out

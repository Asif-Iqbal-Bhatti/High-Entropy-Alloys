"""
Tests for lagrangian_strain.py — I/O and deformation-code mapping.
"""
import textwrap
from pathlib import Path

import numpy as np
import pytest
from conftest import BCC_FE_CONTCAR, assert_matrix_close

from lagrangian_strain import (
    DEFORMATION_CODES,
    read_contcar,
    write_poscar,
)
from _strain_core import build_voigt_strains, voigt_to_eta_matrix, lagrangian_to_linear, deform_cell


# ===========================================================================
# DEFORMATION_CODES catalogue
# ===========================================================================

class TestDeformationCodes:

    def test_all_codes_have_6char_dc(self):
        for code, (dc, _) in DEFORMATION_CODES.items():
            assert len(dc) == 6, f"Code {code}: dc '{dc}' is not 6 chars"

    def test_all_dc_characters_valid(self):
        valid = {'E', 'e', '0', '2'}
        for code, (dc, _) in DEFORMATION_CODES.items():
            for ch in dc:
                assert ch in valid, f"Code {code}: invalid char '{ch}' in '{dc}'"

    def test_code_0_is_volumetric(self):
        dc, _ = DEFORMATION_CODES[0]
        e = build_voigt_strains(dc, 0.05)
        assert np.all(e[:3] == pytest.approx(0.05))
        assert np.all(e[3:] == pytest.approx(0.0))

    def test_code_9_xy_shear_antisymmetric(self):
        dc, _ = DEFORMATION_CODES[9]   # 'Ee0000'
        e = build_voigt_strains(dc, 0.03)
        assert e[0] == pytest.approx(0.03)
        assert e[1] == pytest.approx(-0.03)
        assert np.all(e[2:] == pytest.approx(0.0))

    def test_all_codes_produce_valid_strains(self):
        for code, (dc, _) in DEFORMATION_CODES.items():
            e   = build_voigt_strains(dc, 0.05)
            eta = voigt_to_eta_matrix(e)
            eps = lagrangian_to_linear(eta)
            # Roundtrip: η = ε + ½ε²
            eta_re = eps + 0.5 * eps @ eps
            assert_matrix_close(eta_re, eta, atol=1e-9,
                                ), f"Roundtrip failed for code {code}"


# ===========================================================================
# read_contcar
# ===========================================================================

class TestReadContcar:

    def test_reads_title(self, bcc_fe_contcar):
        d = read_contcar(bcc_fe_contcar)
        assert "Fe" in d["title"]

    def test_reads_scale(self, bcc_fe_contcar):
        d = read_contcar(bcc_fe_contcar)
        assert d["scale"] == pytest.approx(1.0)

    def test_reads_cubic_axis_matrix(self, bcc_fe_contcar):
        d = read_contcar(bcc_fe_contcar)
        expected = np.eye(3) * 2.87
        assert_matrix_close(d["axis_matrix"], expected, atol=1e-10)

    def test_reads_nat(self, bcc_fe_contcar):
        d = read_contcar(bcc_fe_contcar)
        assert d["nat"] == [1]

    def test_tail_lines_start_with_direct(self, bcc_fe_contcar):
        d = read_contcar(bcc_fe_contcar)
        assert "direct" in d["tail_lines"][0].lower()

    def test_rejects_short_file(self, tmp_path):
        short = tmp_path / "CONTCAR"
        short.write_text("only\ntwo\nlines\n")
        with pytest.raises(ValueError, match="fewer than 8"):
            read_contcar(short)

    def test_volume_from_axis(self, bcc_fe_contcar):
        d    = read_contcar(bcc_fe_contcar)
        vol  = abs(np.linalg.det(d["axis_matrix"]) * d["scale"]**3)
        # BCC Fe: a = 2.87 Å → V = 2.87³ ≈ 23.64 Å³
        assert vol == pytest.approx(2.87**3, rel=1e-6)


# ===========================================================================
# write_poscar / round-trip
# ===========================================================================

class TestWritePoscar:

    def test_roundtrip_scale_preserved(self, bcc_fe_contcar, tmp_path):
        d   = read_contcar(bcc_fe_contcar)
        out = tmp_path / "POSCAR_out"
        write_poscar(out, d["title"], d["scale"], d["axis_matrix"],
                     d["element_line"], d["atom_number_line"], d["tail_lines"])
        d2 = read_contcar(out)
        assert d2["scale"] == pytest.approx(d["scale"])

    def test_roundtrip_axis_preserved(self, bcc_fe_contcar, tmp_path):
        d   = read_contcar(bcc_fe_contcar)
        out = tmp_path / "POSCAR_out"
        write_poscar(out, d["title"], d["scale"], d["axis_matrix"],
                     d["element_line"], d["atom_number_line"], d["tail_lines"])
        d2 = read_contcar(out)
        assert_matrix_close(d2["axis_matrix"], d["axis_matrix"])

    def test_deformed_axes_written_correctly(self, bcc_fe_contcar, tmp_path):
        d    = read_contcar(bcc_fe_contcar)
        axis = d["axis_matrix"]
        # Apply 5 % volumetric strain
        e    = build_voigt_strains("EEE000", 0.05)
        eps  = lagrangian_to_linear(voigt_to_eta_matrix(e))
        D    = np.eye(3) + eps
        new  = deform_cell(axis, D)

        out  = tmp_path / "POSCAR_def"
        write_poscar(out, d["title"], d["scale"], new,
                     d["element_line"], d["atom_number_line"], d["tail_lines"])
        d2   = read_contcar(out)
        assert_matrix_close(d2["axis_matrix"], new, atol=1e-12)

    def test_output_has_correct_number_of_lines(self, bcc_fe_contcar, tmp_path):
        d   = read_contcar(bcc_fe_contcar)
        out = tmp_path / "POSCAR_out"
        write_poscar(out, d["title"], d["scale"], d["axis_matrix"],
                     d["element_line"], d["atom_number_line"], d["tail_lines"])
        lines = out.read_text().splitlines()
        # title + scale + 3 vectors + element + natom + coord_type + 1 atom = 9
        assert len(lines) >= 9

    def test_tail_lines_unchanged(self, bcc_fe_contcar, tmp_path):
        d     = read_contcar(bcc_fe_contcar)
        out   = tmp_path / "POSCAR_out"
        write_poscar(out, d["title"], d["scale"], d["axis_matrix"],
                     d["element_line"], d["atom_number_line"], d["tail_lines"])
        d2    = read_contcar(out)
        # Coordinate-type line should still be preserved
        assert "direct" in d2["tail_lines"][0].lower()

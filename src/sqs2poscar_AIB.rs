//! sqs2poscar — Rust port of sqs2poscar_modified_AIB.cpp
//!
//! Usage:  sqs2poscar <bestsqs.out>
//!
//! Reads a standard ATAT 5.0 bestsqs.out file and writes a VASP POSCAR
//! named after the reduced composition (e.g. Li3PS4Cl1.vasp).
//!
//! Original C++ by Changning Niu (2013); modified by Asif Iqbal Bhatti.
//! Rust translation preserves all algorithmic behaviour.

use std::env;
use std::fmt::Write as FmtWrite;
use std::fs;
use std::io::{self, BufRead, Write as IoWrite};
use std::process;

// ─── helpers ────────────────────────────────────────────────────────────────

fn gcd(mut a: usize, mut b: usize) -> usize {
    while b != 0 {
        let t = a % b;
        a = b;
        b = t;
    }
    a
}

fn gcd_all(v: &[usize]) -> usize {
    v.iter().copied().fold(0, gcd).max(1)
}

/// Build a composition string like "Li3_P1_S4_Cl1" (underscore-separated).
/// If `reduced` is true, divide all counts by their GCD first.
fn build_composition_string(elems: &[String], counts: &[usize], reduced: bool) -> String {
    let g = if reduced { gcd_all(counts) } else { 1 };
    let mut out = String::new();
    for (i, (e, c)) in elems.iter().zip(counts.iter()).enumerate() {
        if i > 0 {
            out.push('_');
        }
        write!(out, "{}{}", e, c / g).unwrap();
    }
    out
}

// ─── matrix helpers ─────────────────────────────────────────────────────────

type Mat3 = [[f64; 3]; 3];

#[allow(dead_code)]
fn mat_mul_vec(m: &Mat3, v: [f64; 3]) -> [f64; 3] {
    [
        m[0][0] * v[0] + m[0][1] * v[1] + m[0][2] * v[2],
        m[1][0] * v[0] + m[1][1] * v[1] + m[1][2] * v[2],
        m[2][0] * v[0] + m[2][1] * v[1] + m[2][2] * v[2],
    ]
}

/// Compute the inverse of a 3×3 matrix.  Returns None if singular.
fn inverse_matrix(m: &Mat3) -> Option<Mat3> {
    let det = m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
        - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
        + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);
    if det == 0.0 {
        return None;
    }
    let inv_det = 1.0 / det;
    let mut inv = [[0f64; 3]; 3];
    inv[0][0] = (m[1][1] * m[2][2] - m[1][2] * m[2][1]) * inv_det;
    inv[0][1] = (m[0][2] * m[2][1] - m[0][1] * m[2][2]) * inv_det;
    inv[0][2] = (m[0][1] * m[1][2] - m[0][2] * m[1][1]) * inv_det;
    inv[1][0] = (m[1][2] * m[2][0] - m[1][0] * m[2][2]) * inv_det;
    inv[1][1] = (m[0][0] * m[2][2] - m[0][2] * m[2][0]) * inv_det;
    inv[1][2] = (m[0][2] * m[1][0] - m[0][0] * m[1][2]) * inv_det;
    inv[2][0] = (m[1][0] * m[2][1] - m[1][1] * m[2][0]) * inv_det;
    inv[2][1] = (m[0][1] * m[2][0] - m[0][0] * m[2][1]) * inv_det;
    inv[2][2] = (m[0][0] * m[1][1] - m[0][1] * m[1][0]) * inv_det;
    Some(inv)
}

// ─── bestsqs parser ─────────────────────────────────────────────────────────

struct BestSqs {
    vec1: Mat3,                   // basis coordinate system (lines 1-3)
    vec2: Mat3,                   // lattice vectors scaled by vec1 (lines 4-6)
    atom_coords: Vec<[f64; 3]>,   // atomic positions (lines 7+)
    atom_species: Vec<String>,    // element labels matching atom_coords
    elem_names: Vec<String>,      // unique element names (insertion order)
    elem_counts: Vec<usize>,      // counts matching elem_names
}

fn read_bestsqs(path: &str) -> io::Result<BestSqs> {
    let file = fs::File::open(path)?;
    let mut lines = io::BufReader::new(file).lines();

    // Helper: parse one "x y z" line
    let read_vec3 = |lines: &mut dyn Iterator<Item = io::Result<String>>| -> io::Result<[f64; 3]> {
        let line = lines
            .next()
            .ok_or_else(|| io::Error::new(io::ErrorKind::UnexpectedEof, "unexpected EOF"))??;
        let mut it = line.split_whitespace();
        let parse = |s: Option<&str>| -> io::Result<f64> {
            s.ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "missing field"))?
                .parse::<f64>()
                .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
        };
        Ok([parse(it.next())?, parse(it.next())?, parse(it.next())?])
    };

    // Lines 1-3: basis
    let mut vec1 = [[0f64; 3]; 3];
    for row in vec1.iter_mut() {
        *row = read_vec3(&mut lines)?;
    }

    // Lines 4-6: supercell lattice vectors (in basis coordinates)
    let mut vec2 = [[0f64; 3]; 3];
    for row in vec2.iter_mut() {
        *row = read_vec3(&mut lines)?;
    }

    // Lines 7+: atom positions + species
    let mut atom_coords: Vec<[f64; 3]> = Vec::new();
    let mut atom_species: Vec<String> = Vec::new();

    for line in lines {
        let line = line?;
        let parts: Vec<&str> = line.split_whitespace().collect();
        if parts.len() < 4 {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!("malformed atom line: '{}'", line),
            ));
        }
        let x: f64 = parts[0].parse().map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
        let y: f64 = parts[1].parse().map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
        let z: f64 = parts[2].parse().map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
        atom_coords.push([x, y, z]);
        atom_species.push(parts[3].to_string());
    }

    // Build unique element list (preserving first-seen order, matching C++ behaviour)
    let mut elem_names: Vec<String> = Vec::new();
    let mut elem_counts: Vec<usize> = Vec::new();

    for sp in &atom_species {
        if let Some(idx) = elem_names.iter().position(|e| e == sp) {
            elem_counts[idx] += 1;
        } else {
            elem_names.push(sp.clone());
            elem_counts.push(1);
        }
    }

    Ok(BestSqs {
        vec1,
        vec2,
        atom_coords,
        atom_species,
        elem_names,
        elem_counts,
    })
}

// ─── main ────────────────────────────────────────────────────────────────────

fn main() {
    let args: Vec<String> = env::args().collect();
    if args.len() != 2 {
        eprintln!("Usage: sqs2poscar <FileName>");
        process::exit(1);
    }
    let file = &args[1];

    let sqs = match read_bestsqs(file) {
        Ok(s) => s,
        Err(e) => {
            eprintln!("Error reading {}: {}", file, e);
            process::exit(1);
        }
    };

    // ── Cartesian lattice vectors: lat[i] = vec2[i] · vec1  ─────────────────
    // (vec2 rows are in basis-coordinate space; vec1 rows are basis vectors)
    // lat[i][j] = sum_k  vec2[i][k] * vec1[k][j]
    let mut lat_vec: Mat3 = [[0f64; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            lat_vec[i][j] = sqs.vec2[i][0] * sqs.vec1[0][j]
                + sqs.vec2[i][1] * sqs.vec1[1][j]
                + sqs.vec2[i][2] * sqs.vec1[2][j];
        }
    }

    // ── Cartesian atomic positions: pos_cart[i] = vec1^T · coord[i]  ────────
    // Original C++ multiplies [x,y,z] · vec1 where vec1 rows are basis vecs,
    // which gives  pos[j] = sum_k coord[k] * vec1[k][j].
    let atom_cart: Vec<[f64; 3]> = sqs
        .atom_coords
        .iter()
        .map(|c| {
            [
                c[0] * sqs.vec1[0][0] + c[1] * sqs.vec1[1][0] + c[2] * sqs.vec1[2][0],
                c[0] * sqs.vec1[0][1] + c[1] * sqs.vec1[1][1] + c[2] * sqs.vec1[2][1],
                c[0] * sqs.vec1[0][2] + c[1] * sqs.vec1[1][2] + c[2] * sqs.vec1[2][2],
            ]
        })
        .collect();

    // ── Sort atoms by element (same order as elem_names) ────────────────────
    let mut atom_sorted: Vec<[f64; 3]> = Vec::with_capacity(atom_cart.len());
    for elem in &sqs.elem_names {
        for (j, sp) in sqs.atom_species.iter().enumerate() {
            if sp == elem {
                atom_sorted.push(atom_cart[j]);
            }
        }
    }
    assert_eq!(atom_sorted.len(), atom_cart.len(), "sort count mismatch");

    // ── Convert to fractional (Direct) coordinates ───────────────────────────
    // frac[i] = lat_vec_inv^T · cart[i]  (same transposed-inverse convention as C++)
    let lat_inv = match inverse_matrix(&lat_vec) {
        Some(m) => m,
        None => {
            eprintln!("Singular lattice matrix — cannot convert to Direct coordinates.");
            process::exit(1);
        }
    };

    // The C++ uses the same transposed-multiply convention for fractional coords:
    // frac[j] = sum_k cart[k] * lat_inv[k][j]
    let atom_frac: Vec<[f64; 3]> = atom_sorted
        .iter()
        .map(|c| {
            [
                c[0] * lat_inv[0][0] + c[1] * lat_inv[1][0] + c[2] * lat_inv[2][0],
                c[0] * lat_inv[0][1] + c[1] * lat_inv[1][1] + c[2] * lat_inv[2][1],
                c[0] * lat_inv[0][2] + c[1] * lat_inv[1][2] + c[2] * lat_inv[2][2],
            ]
        })
        .collect();

    // ── Write POSCAR ─────────────────────────────────────────────────────────
    let outname =
        build_composition_string(&sqs.elem_names, &sqs.elem_counts, true) + ".vasp";

    let mut out = match fs::File::create(&outname) {
        Ok(f) => io::BufWriter::new(f),
        Err(e) => {
            eprintln!("Cannot create {}: {}", outname, e);
            process::exit(1);
        }
    };

    writeln!(out, "POSCAR").unwrap();
    writeln!(out, "3.0").unwrap();

    for row in &lat_vec {
        writeln!(
            out,
            "{:12.8}{:12.8}{:12.8}",
            row[0], row[1], row[2]
        )
        .unwrap();
    }

    // Element names line
    for e in &sqs.elem_names {
        write!(out, "{:>4}", e).unwrap();
    }
    writeln!(out).unwrap();

    // Element counts line
    for c in &sqs.elem_counts {
        write!(out, "{:>4}", c).unwrap();
    }
    writeln!(out).unwrap();

    writeln!(out, "Direct").unwrap();

    for pos in &atom_frac {
        writeln!(
            out,
            "{:12.8}{:12.8}{:12.8}",
            pos[0], pos[1], pos[2]
        )
        .unwrap();
    }

    println!("Written: {}", outname);
}

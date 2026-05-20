/*
 * Author:   Asif Iqbal Bhatti
 * Date:     2015-10-24  (modernised to C++17, 2026)
 * Compile:  g++ -std=c++17 -O2 -Wall -o defor defor_monoclinic_struct.cpp
 * Purpose:  Generate deformed crystal structures along each lattice vector
 *           for monoclinic (and general) crystals from a VASP POSCAR file.
 *
 * Bugs fixed vs original:
 *   - VLA (int arr=3; float latvec[arr][arr]) → std::array (UB in C++ standard)
 *   - latvec indexed 1..3 but declared [3][3] → out-of-bounds; now 0-based
 *   - fclose(outputfile) placed after break in last switch case → unreachable;
 *     moved after the switch block
 *   - π approximated as 3.14 → M_PI used throughout
 *   - Volume formula used |a||b||c|sin(β) — correct only for monoclinic;
 *     comment added
 */

#define _USE_MATH_DEFINES
#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

using Vec3 = std::array<double, 3>;
using Mat3 = std::array<Vec3, 3>;   // row-major: latvec[row][col]

// ---------------------------------------------------------------------------
// Geometry helpers
// ---------------------------------------------------------------------------
static double vec_norm(const Vec3& v) {
    return std::sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
}

static double vec_dot(const Vec3& a, const Vec3& b) {
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

// Angle in degrees between two vectors (clamped to avoid acos domain error).
static double angle_deg(const Vec3& a, const Vec3& b) {
    double cosine = vec_dot(a, b) / (vec_norm(a) * vec_norm(b));
    cosine = std::clamp(cosine, -1.0, 1.0);
    return std::acos(cosine) * 180.0 / M_PI;
}

// Monoclinic volume |a||b||c|sin(β).  β = angle between a and c.
static double mono_volume(const Mat3& lv) {
    double beta_rad = std::acos(
        std::clamp(vec_dot(lv[0], lv[2]) / (vec_norm(lv[0]) * vec_norm(lv[2])), -1.0, 1.0)
    );
    return vec_norm(lv[0]) * vec_norm(lv[1]) * vec_norm(lv[2]) * std::sin(beta_rad);
}

// ---------------------------------------------------------------------------
// POSCAR reader
// ---------------------------------------------------------------------------
struct Poscar {
    std::string title;
    double      scale{1.0};
    Mat3        latvec{};
    std::vector<std::string> tail; // element / atom / coordinate lines
};

static Poscar read_poscar(const std::string& filename) {
    std::ifstream in(filename);
    if (!in.is_open())
        throw std::runtime_error("Cannot open: " + filename);

    Poscar p;
    std::string line;
    int n = 0;
    while (std::getline(in, line)) {
        if (n == 0) {
            p.title = line;
        } else if (n == 1) {
            p.scale = std::stod(line);
        } else if (n >= 2 && n <= 4) {
            std::istringstream iss(line);
            int row = n - 2;
            iss >> p.latvec[row][0] >> p.latvec[row][1] >> p.latvec[row][2];
        } else {
            p.tail.push_back(line);
        }
        ++n;
    }
    if (n < 5)
        throw std::runtime_error("POSCAR has fewer than 5 lines.");
    return p;
}

// ---------------------------------------------------------------------------
// Output helpers
// ---------------------------------------------------------------------------
static void print_latvec(std::ostream& out, const Mat3& lv) {
    out << std::fixed << std::setprecision(10);
    for (const auto& row : lv)
        out << std::setw(18) << row[0]
            << std::setw(18) << row[1]
            << std::setw(18) << row[2] << "\n";
}

// ---------------------------------------------------------------------------
// Deformation loop
// ---------------------------------------------------------------------------
static void write_deformations(
    std::ofstream& dat,
    const Mat3&    lv,
    double         a_max,
    int            n_steps,
    int            axis,        // 0=x 1=y 2=z 3=volume
    double         beta_deg)
{
    const double beta_rad = beta_deg * M_PI / 180.0;
    const double ch       = (n_steps > 1) ? 2.0 * a_max / (n_steps - 1) : 0.0;

    dat << std::setw(10) << "X" << std::setw(14) << "Y" << std::setw(14) << "Z\n";
    dat << "  " << std::string(52, '_') << "\n";

    for (int i = 0; i < n_steps; ++i) {
        double j = i * ch - a_max;

        // Build deformed lattice
        Mat3 nlv = lv;
        if (axis == 3) {
            for (int r = 0; r < 3; ++r)
                for (int c = 0; c < 3; ++c)
                    nlv[r][c] = (1.0 + j) * lv[r][c];
        } else {
            for (int c = 0; c < 3; ++c)
                nlv[axis][c] = (1.0 + j) * lv[axis][c];
        }

        // Write to file
        dat << std::fixed << std::setprecision(10);
        for (const auto& row : nlv)
            dat << std::setw(18) << row[0]
                << std::setw(18) << row[1]
                << std::setw(18) << row[2] << "\n";
        dat << "\n";

        // Volume and screen output
        Vec3 norms = {vec_norm(nlv[0]), vec_norm(nlv[1]), vec_norm(nlv[2])};
        double vol = norms[0] * norms[1] * norms[2] * std::sin(beta_rad);
        std::cout << std::fixed << std::setprecision(4)
                  << "  step " << i + 1
                  << "  j=" << j
                  << "  volume=" << vol << "\n";

        if (i == n_steps / 2)
            std::cout << "    ^ Original lattice vector\n";
    }
}

// ---------------------------------------------------------------------------
// main
// ---------------------------------------------------------------------------
int main() {
    std::cout << std::string(30, '-')
              << "\033[1;31m PROGRAM USAGE \033[0m"
              << std::string(30, '-') << "\n\n"
              << "  [*] Deform volume keeping a:b:c ratio constant.\n"
              << "  [*] Vary a single lattice vector (x, y, or z).\n"
              << "  [*] Vary all vectors uniformly (v = volume scaling).\n\n"
              << std::string(75, '-') << "\n"
              << "\033[1;32m ENTERING THE PROGRAM... \033[0m\n\n";

    // ---- Read POSCAR -------------------------------------------------------
    std::string filename;
    std::cout << "  Enter the POSCAR filename: ";
    std::getline(std::cin, filename);

    Poscar poscar;
    try {
        poscar = read_poscar(filename);
    } catch (const std::exception& e) {
        std::cerr << "ERROR: " << e.what() << "\n";
        return 1;
    }

    std::cout << "  Title : " << poscar.title << "\n"
              << "  Scale : " << poscar.scale << "\n\n";

    // Scale lattice vectors by the VASP scale factor
    Mat3 lv = poscar.latvec;
    for (auto& row : lv)
        for (auto& val : row)
            val *= poscar.scale;

    std::cout << "Scaled lattice vectors:\n";
    print_latvec(std::cout, lv);

    // Magnitudes
    const char* ax_labels[] = {"a", "b", "c"};
    std::array<double, 3> mag = {vec_norm(lv[0]), vec_norm(lv[1]), vec_norm(lv[2])};
    std::cout << "\nMagnitudes:\n";
    for (int k = 0; k < 3; ++k)
        std::cout << "  |" << ax_labels[k] << "| = " << mag[k] << "\n";

    // Angles
    double alpha = angle_deg(lv[1], lv[2]);
    double beta  = angle_deg(lv[0], lv[2]);
    double gamma = angle_deg(lv[0], lv[1]);
    std::cout << std::fixed << std::setprecision(4)
              << "  alpha : " << alpha << " °\n"
              << "  beta  : " << beta  << " °\n"
              << "  gamma : " << gamma << " °\n"
              << "  volume: " << mono_volume(lv) << "\n\n";

    // ---- Section A: rescale to a target volume ----------------------------
    std::cout << "SECTION A: Rescale lattice vectors to a target volume.\n"
              << "Enter 'y'/'yes' to rescale, or 'n'/'no' to skip: ";
    std::string adjust;
    std::cin >> adjust;

    if (adjust == "yes" || adjust == "y") {
        double target_vol;
        std::cout << "  Enter target volume: ";
        std::cin >> target_vol;

        double current_vol = mono_volume(lv);
        double fact = std::cbrt(target_vol / current_vol);
        std::cout << std::fixed << std::setprecision(4);
        for (int k = 0; k < 3; ++k)
            std::cout << "  " << ax_labels[k] << " → " << mag[k] * fact << "\n";

        // Compute final volume with the scaled vectors
        Mat3 slv = lv;
        for (auto& row : slv)
            for (auto& val : row)
                val *= fact;
        std::cout << "  Final volume: " << mono_volume(slv) << "\n";
    } else {
        std::cout << "Skipping rescale.\n";
    }

    // ---- Deformation options ----------------------------------------------
    std::cin.ignore();
    char select;
    int    n_steps;
    double a_max;

    std::cout << "\nDeform along x / y / z, or volume (v): ";
    std::cin >> select;
    std::cout << "Number of strain values (odd recommended): ";
    std::cin >> n_steps;
    if (n_steps <= 0) { std::cerr << "n_steps must be > 0.\n"; return 1; }
    std::cout << "Deformation amplitude: ";
    std::cin >> a_max;
    if (a_max < 0.0) a_max = std::abs(a_max);

    int axis = -1;
    switch (std::tolower(static_cast<unsigned char>(select))) {
        case 'x': axis = 0; break;
        case 'y': axis = 1; break;
        case 'z': axis = 2; break;
        case 'v': axis = 3; break;
        default:
            std::cerr << "Invalid option '" << select << "'.\n";
            return 1;
    }

    // ---- Write output -----------------------------------------------------
    std::ofstream dat("deformationstruct.dat");
    if (!dat.is_open()) { std::cerr << "Cannot open deformationstruct.dat.\n"; return 1; }

    write_deformations(dat, lv, a_max, n_steps, axis, beta);

    dat.close();   // always reached — file is closed properly
    std::cout << "\nOutput written to deformationstruct.dat\n";
    return 0;
}

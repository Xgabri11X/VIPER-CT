/*
 * VIPER-CT
 * Voronoi Integration and Partitioning of Electron Redistribution for Charge Transfer
 *
 * Grid-based integration of Gaussian CUBE scalar fields over atom-centered,
 * optionally vdW-weighted Voronoi regions.
 *
 * Supported workflows
 * -------------------
 * 1) Single-CUBE mode
 *      viper_ct [options] field.cube
 *
 *    The scalar field stored in field.cube is integrated directly.
 *
 * 2) Three-CUBE difference mode
 *      viper_ct [options] combined.cube fragment1.cube fragment2.cube
 *
 *    VIPER-CT constructs, voxel by voxel,
 *
 *      delta(r) = combined(r) - fragment1(r) - fragment2(r)
 *
 *    and integrates delta(r) over atom-centered regions defined from the
 *    geometry stored in combined.cube.
 *
 * Partitioning
 * ------------
 * Unweighted mode (default):
 *      owner(P) = argmin_a |P - R_a|
 *
 * vdW-weighted mode (--weighted vdw):
 *      owner(P) = argmin_a |P - R_a| / r_vdw(a)
 *
 * After the nearest/weighted-nearest atom is selected, a voxel is retained only if
 *      |P - R_owner| <= cutoff_scale * r_vdw(owner)
 *
 * Optional periodic boundary conditions use the minimum-image distance in the
 * triclinic cell reconstructed from the CUBE grid vectors.
 *
 * Output
 * ------
 * CSV written to stdout:
 *      idx,Z,x,y,z,Integral,Positive,Negative,AbsIntegral
 *
 * where AbsIntegral is the integral of the absolute value of the voxel
 * contributions and is NOT simply abs(Integral).
 *
 * Build example
 * -------------
 *      g++ -O3 -march=native -std=c++17 -fopenmp src/viper_ct.cpp -o viper_ct
 */

#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <locale>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace std;

// ============================================================================
// Basic data structures
// ============================================================================

struct Atom {
    int Z;
    double x, y, z;
};

struct Cube {
    int natoms = 0;
    array<double, 3> origin{};
    int Nx = 0, Ny = 0, Nz = 0;
    array<double, 3> dX{}, dY{}, dZ{};
    vector<Atom> atoms;
    vector<double> grid;
};

// 1 Angstrom in Bohr.
static constexpr double ANG2BOHR = 1.8897259886;

// ============================================================================
// van der Waals radii
// ============================================================================

/*
 * van der Waals radii used by the weighted partition and by the cutoff.
 * Values are stored in Angstrom and converted internally to Bohr.
 *
 * To add a new element, add one entry to this table using
 *      {atomic_number, radius_in_angstrom}
 *
 * Example syntax only:
 *      {Z, radius_A},
 *
 * Unknown elements fall back to 1.50 Angstrom. VIPER-CT prints a warning to
 * stderr when this fallback is used so that published calculations do not
 * silently depend on an unintended radius.
 */
static const unordered_map<int, double>& vdw_radii_angstrom() {
    static const unordered_map<int, double> radii = {
        {1, 1.20},   // H
        {6, 1.70},   // C
        {7, 1.55},   // N
        {8, 1.52},   // O
        {9, 1.47},   // F
        {11, 2.27},  // Na
        {12, 1.73},  // Mg
        {14, 2.10},  // Si
        {15, 1.80},  // P
        {16, 1.80},  // S
        {17, 1.75},  // Cl
        {19, 2.75},  // K
        {20, 2.31},  // Ca
        {35, 1.85},  // Br
        {53, 1.98}   // I
    };
    return radii;
}

static bool has_vdw_radius(int Z) {
    return vdw_radii_angstrom().find(Z) != vdw_radii_angstrom().end();
}

static double vdw_radius_bohr(int Z) {
    const auto& radii = vdw_radii_angstrom();
    const auto it = radii.find(Z);
    const double radius_angstrom = (it == radii.end()) ? 1.50 : it->second;
    return radius_angstrom * ANG2BOHR;
}

// ============================================================================
// Small vector helpers
// ============================================================================

static inline array<double, 3> add(const array<double, 3>& a,
                                   const array<double, 3>& b) {
    return {a[0] + b[0], a[1] + b[1], a[2] + b[2]};
}

static inline array<double, 3> smul(const array<double, 3>& a, double s) {
    return {a[0] * s, a[1] * s, a[2] * s};
}

static inline double dot3(const array<double, 3>& a,
                          const array<double, 3>& b) {
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

static inline array<double, 3> cross3(const array<double, 3>& a,
                                      const array<double, 3>& b) {
    return {
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0]
    };
}

static inline double norm(const array<double, 3>& a) {
    return sqrt(dot3(a, a));
}

// Cartesian position of grid node (i,j,k).
static inline array<double, 3> node_pos(const Cube& cube, int i, int j, int k) {
    return add(cube.origin,
               add(smul(cube.dX, i),
                   add(smul(cube.dY, j), smul(cube.dZ, k))));
}

// Volume associated with one grid voxel for a general triclinic grid.
static inline double voxel_volume(const Cube& cube) {
    const auto dX_cross_dY = cross3(cube.dX, cube.dY);
    return fabs(dot3(dX_cross_dY, cube.dZ));
}

// ============================================================================
// Gaussian CUBE I/O
// ============================================================================

static Cube read_cube(const string& path) {
    ifstream fin(path);
    if (!fin) {
        throw runtime_error("Cannot open CUBE file: " + path);
    }

    // Standard CUBE files begin with two comment lines.
    string line;
    getline(fin, line);
    getline(fin, line);

    Cube cube;
    fin >> cube.natoms >> cube.origin[0] >> cube.origin[1] >> cube.origin[2];
    fin >> cube.Nx >> cube.dX[0] >> cube.dX[1] >> cube.dX[2];
    fin >> cube.Ny >> cube.dY[0] >> cube.dY[1] >> cube.dY[2];
    fin >> cube.Nz >> cube.dZ[0] >> cube.dZ[1] >> cube.dZ[2];

    if (!fin) {
        throw runtime_error("Malformed CUBE header in: " + path);
    }
    if (cube.natoms < 0) {
        throw runtime_error(
            "Negative atom counts / orbital-style CUBE files are not supported: " + path
        );
    }
    if (cube.Nx <= 0 || cube.Ny <= 0 || cube.Nz <= 0) {
        throw runtime_error(
            "VIPER-CT expects positive CUBE grid counts and coordinates in Bohr: " + path
        );
    }

    cube.atoms.resize(cube.natoms);
    for (int i = 0; i < cube.natoms; ++i) {
        int Z;
        double nuclear_charge, x, y, z;
        fin >> Z >> nuclear_charge >> x >> y >> z;
        if (!fin) {
            throw runtime_error("Malformed atom list in: " + path);
        }
        cube.atoms[i] = {Z, x, y, z};
    }

    const size_t expected = static_cast<size_t>(cube.Nx)
                          * static_cast<size_t>(cube.Ny)
                          * static_cast<size_t>(cube.Nz);

    cube.grid.reserve(expected);
    double value;
    while (fin >> value) {
        cube.grid.push_back(value);
    }

    if (cube.grid.size() != expected) {
        throw runtime_error(
            "Unexpected number of grid values in " + path +
            ": expected " + to_string(expected) +
            ", read " + to_string(cube.grid.size())
        );
    }

    return cube;
}

static inline bool approx_equal(double a, double b, double eps = 1e-8) {
    return fabs(a - b) < eps;
}

static inline bool approx_equal_vec(const array<double, 3>& a,
                                    const array<double, 3>& b,
                                    double eps = 1e-8) {
    return approx_equal(a[0], b[0], eps)
        && approx_equal(a[1], b[1], eps)
        && approx_equal(a[2], b[2], eps);
}

/*
 * The three fields must live on exactly the same real-space grid.
 * Fragment CUBE files do not need to contain the same atom list because the
 * partitioning geometry is taken exclusively from the first (combined) CUBE.
 */
static void ensure_compatible(const Cube& A, const Cube& B, const Cube& C) {
    auto incompatible = [](const string& message) {
        throw runtime_error("The three CUBE grids are incompatible: " + message);
    };

    if (A.Nx != B.Nx || A.Ny != B.Ny || A.Nz != B.Nz
        || A.Nx != C.Nx || A.Ny != C.Ny || A.Nz != C.Nz) {
        incompatible("different Nx/Ny/Nz dimensions");
    }

    if (!approx_equal_vec(A.origin, B.origin)
        || !approx_equal_vec(A.origin, C.origin)) {
        incompatible("different grid origins");
    }

    if (!approx_equal_vec(A.dX, B.dX)
        || !approx_equal_vec(A.dY, B.dY)
        || !approx_equal_vec(A.dZ, B.dZ)
        || !approx_equal_vec(A.dX, C.dX)
        || !approx_equal_vec(A.dY, C.dY)
        || !approx_equal_vec(A.dZ, C.dZ)) {
        incompatible("different grid vectors");
    }
}

/*
 * Write a CUBE file containing arbitrary double values on the reference grid.
 * For owner maps, voxel values are atom_index+1, while 0 means unassigned.
 */
static void write_cube(const string& path,
                       const Cube& ref,
                       const vector<double>& data,
                       const string& title1 = "Owners (atom_index+1; 0=unassigned)",
                       const string& title2 = "Voronoi partition on CUBE grid") {
    ofstream fout(path, ios::binary);
    if (!fout) {
        throw runtime_error("Cannot write CUBE file: " + path);
    }
    if (data.size() != ref.grid.size()) {
        throw runtime_error("Internal error: owner-map size does not match CUBE grid.");
    }

    fout.imbue(locale::classic());
    fout << title1 << '\n' << title2 << '\n';

    fout.setf(ios::fixed);
    fout << setw(5) << ref.natoms
         << setw(13) << setprecision(6) << ref.origin[0]
         << setw(13) << setprecision(6) << ref.origin[1]
         << setw(13) << setprecision(6) << ref.origin[2] << '\n';

    fout << setw(5) << ref.Nx
         << setw(13) << setprecision(6) << ref.dX[0]
         << setw(13) << setprecision(6) << ref.dX[1]
         << setw(13) << setprecision(6) << ref.dX[2] << '\n';

    fout << setw(5) << ref.Ny
         << setw(13) << setprecision(6) << ref.dY[0]
         << setw(13) << setprecision(6) << ref.dY[1]
         << setw(13) << setprecision(6) << ref.dY[2] << '\n';

    fout << setw(5) << ref.Nz
         << setw(13) << setprecision(6) << ref.dZ[0]
         << setw(13) << setprecision(6) << ref.dZ[1]
         << setw(13) << setprecision(6) << ref.dZ[2] << '\n';

    for (const auto& atom : ref.atoms) {
        fout << setw(5) << atom.Z
             << setw(13) << setprecision(6) << 0.0
             << setw(13) << setprecision(6) << atom.x
             << setw(13) << setprecision(6) << atom.y
             << setw(13) << setprecision(6) << atom.z << '\n';
    }

    fout.setf(ios::scientific);
    fout << setprecision(5);

    constexpr int values_per_line = 6;
    int column = 0;
    for (double value : data) {
        fout << setw(13) << value;
        if (++column == values_per_line) {
            fout << '\n';
            column = 0;
        }
    }
    if (column != 0) {
        fout << '\n';
    }
}

// ============================================================================
// Periodic boundary conditions
// ============================================================================

struct PBC {
    // Direct lattice vectors.
    array<double, 3> A, B, C;

    // Reciprocal rows used to convert Cartesian displacement -> fractional.
    array<double, 3> Astar, Bstar, Cstar;

    double det = 0.0;
    bool enabled = false;
};

static inline array<double, 3> frac_wrap(const array<double, 3>& s) {
    // Wrap each fractional coordinate to the minimum-image interval.
    return {
        s[0] - nearbyint(s[0]),
        s[1] - nearbyint(s[1]),
        s[2] - nearbyint(s[2])
    };
}

static inline array<double, 3> frac2cart(const PBC& p,
                                         const array<double, 3>& s) {
    return add(smul(p.A, s[0]), add(smul(p.B, s[1]), smul(p.C, s[2])));
}

static PBC make_pbc_from_cube(const Cube& cube, bool enable) {
    PBC p;
    p.enabled = enable;
    if (!enable) {
        return p;
    }

    // Reconstruct the full cell vectors from the CUBE voxel vectors.
    p.A = smul(cube.dX, cube.Nx);
    p.B = smul(cube.dY, cube.Ny);
    p.C = smul(cube.dZ, cube.Nz);

    const auto BxC = cross3(p.B, p.C);
    p.det = dot3(p.A, BxC);
    if (fabs(p.det) < 1e-12) {
        throw runtime_error("PBC error: nearly singular cell matrix (det ~ 0).");
    }

    // Reciprocal basis rows: s = (A*·r, B*·r, C*·r).
    p.Astar = {BxC[0] / p.det, BxC[1] / p.det, BxC[2] / p.det};

    const auto CxA = cross3(p.C, p.A);
    p.Bstar = {CxA[0] / p.det, CxA[1] / p.det, CxA[2] / p.det};

    const auto AxB = cross3(p.A, p.B);
    p.Cstar = {AxB[0] / p.det, AxB[1] / p.det, AxB[2] / p.det};

    return p;
}

static double pbc_distance(const PBC& p, const array<double, 3>& displacement) {
    if (!p.enabled) {
        return norm(displacement);
    }

    const array<double, 3> fractional = {
        dot3(p.Astar, displacement),
        dot3(p.Bstar, displacement),
        dot3(p.Cstar, displacement)
    };

    return norm(frac2cart(p, frac_wrap(fractional)));
}

// ============================================================================
// Command-line interface
// ============================================================================

struct Options {
    string weighted = "none";   // none | vdw
    double cutoff_scale = 2.0;
    string owners_path;
    vector<string> cubes;        // exactly 1 or 3 CUBE files
    bool pbc = false;
    int nthreads = 0;            // 0 -> OpenMP/system default
    bool help = false;
};

static void print_help(const char* program) {
    cout
        << "VIPER-CT - Voronoi Integration and Partitioning of Electron Redistribution for Charge Transfer\n\n"
        << "Usage:\n"
        << "  " << program << " [options] field.cube\n"
        << "  " << program << " [options] combined.cube fragment1.cube fragment2.cube\n\n"
        << "Modes:\n"
        << "  1 CUBE   Integrate the field stored in field.cube directly.\n"
        << "  3 CUBEs  Build combined - fragment1 - fragment2 voxel by voxel, then integrate it.\n\n"
        << "Options:\n"
        << "  --weighted MODE      Partition metric: none (default) or vdw.\n"
        << "  --cutoff-scale S     Retain voxels within S * r_vdw(owner). Default: 2.0.\n"
        << "  --owners PATH        Write an owner-map CUBE (atom_index+1; 0 = unassigned).\n"
        << "  --pbc                Use triclinic minimum-image periodic distances.\n"
        << "  --threads N          Number of OpenMP threads; 0 uses the runtime default.\n"
        << "  -h, --help           Show this help message.\n\n"
        << "CSV output columns:\n"
        << "  idx,Z,x,y,z,Integral,Positive,Negative,AbsIntegral\n";
}

static Options parse_cli(int argc, char** argv) {
    Options opt;
    vector<string> positional;

    for (int i = 1; i < argc; ++i) {
        const string arg = argv[i];

        if (arg == "--weighted") {
            if (i + 1 >= argc) {
                throw runtime_error("Missing value after --weighted.");
            }
            opt.weighted = argv[++i];
        } else if (arg == "--cutoff-scale") {
            if (i + 1 >= argc) {
                throw runtime_error("Missing value after --cutoff-scale.");
            }
            opt.cutoff_scale = stod(argv[++i]);
        } else if (arg == "--owners") {
            if (i + 1 >= argc) {
                throw runtime_error("Missing path after --owners.");
            }
            opt.owners_path = argv[++i];
        } else if (arg == "--pbc") {
            opt.pbc = true;
        } else if (arg == "--threads") {
            if (i + 1 >= argc) {
                throw runtime_error("Missing value after --threads.");
            }
            opt.nthreads = stoi(argv[++i]);
        } else if (arg == "-h" || arg == "--help") {
            opt.help = true;
        } else if (!arg.empty() && arg[0] == '-') {
            throw runtime_error("Unknown option: " + arg);
        } else {
            positional.push_back(arg);
        }
    }

    if (opt.help) {
        return opt;
    }

    if (!(positional.size() == 1 || positional.size() == 3)) {
        throw runtime_error(
            "Expected either one CUBE file or three CUBE files. Use --help for usage."
        );
    }

    opt.cubes = positional;

    transform(opt.weighted.begin(), opt.weighted.end(), opt.weighted.begin(),
              [](unsigned char c) { return static_cast<char>(tolower(c)); });

    if (opt.weighted != "none" && opt.weighted != "vdw") {
        throw runtime_error("Invalid --weighted value. Use: none | vdw.");
    }

    if (opt.cutoff_scale <= 0.0) {
        throw runtime_error("--cutoff-scale must be > 0.");
    }

    if (opt.nthreads < 0) {
        throw runtime_error("--threads must be >= 0.");
    }

    return opt;
}

// ============================================================================
// Main integration workflow
// ============================================================================

int main(int argc, char** argv) {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    try {
        const Options opt = parse_cli(argc, argv);

        if (opt.help) {
            print_help(argv[0]);
            return 0;
        }

        // The first CUBE is always the reference geometry and reference grid.
        const Cube reference = read_cube(opt.cubes[0]);
        const size_t N = static_cast<size_t>(reference.Nx)
                       * static_cast<size_t>(reference.Ny)
                       * static_cast<size_t>(reference.Nz);

        vector<double> field(N);
        const bool single_cube_mode = (opt.cubes.size() == 1);

        if (single_cube_mode) {
            // Directly integrate the scalar field already stored in the CUBE.
            field = reference.grid;
        } else {
            // Construct combined - fragment1 - fragment2 on the common grid.
            const Cube fragment1 = read_cube(opt.cubes[1]);
            const Cube fragment2 = read_cube(opt.cubes[2]);
            ensure_compatible(reference, fragment1, fragment2);

            for (size_t i = 0; i < N; ++i) {
                field[i] = reference.grid[i]
                         - fragment1.grid[i]
                         - fragment2.grid[i];
            }
        }

        // Build periodic-cell information only when requested.
        const PBC pbc = make_pbc_from_cube(reference, opt.pbc);

        // Precompute atom-specific vdW radii and absolute cutoff distances.
        const int Nat = reference.natoms;
        vector<double> r_vdw(Nat), r_cut(Nat);

        for (int a = 0; a < Nat; ++a) {
            const int Z = reference.atoms[a].Z;
            if (!has_vdw_radius(Z)) {
                cerr << "[VIPER-CT warning] No vdW radius is defined for atomic number "
                     << Z << ". Using the 1.50 A fallback.\n";
            }
            r_vdw[a] = vdw_radius_bohr(Z);
            r_cut[a] = opt.cutoff_scale * r_vdw[a];
        }

        const double dV = voxel_volume(reference);

        // Optional owner map: 0 = unassigned, atom index + 1 = assigned atom.
        vector<double> owners;
        const bool write_owner_map = !opt.owners_path.empty();
        if (write_owner_map) {
            owners.assign(N, 0.0);
        }

        // Cache atomic positions for the tight voxel loop.
        vector<array<double, 3>> atom_positions(Nat);
        for (int a = 0; a < Nat; ++a) {
            atom_positions[a] = {
                reference.atoms[a].x,
                reference.atoms[a].y,
                reference.atoms[a].z
            };
        }

        // Global integration buffers. Long double reduces summation noise.
        vector<long double> integral_total(Nat, 0.0L);
        vector<long double> integral_positive(Nat, 0.0L);
        vector<long double> integral_negative(Nat, 0.0L);
        vector<long double> integral_absolute(Nat, 0.0L);

        // Configure OpenMP and allocate one accumulator per thread and atom.
        int thread_count = 1;
#ifdef _OPENMP
        if (opt.nthreads > 0) {
            omp_set_num_threads(opt.nthreads);
        }
        thread_count = omp_get_max_threads();
#endif

        vector<vector<long double>> total_by_thread(
            thread_count, vector<long double>(Nat, 0.0L));
        vector<vector<long double>> positive_by_thread(
            thread_count, vector<long double>(Nat, 0.0L));
        vector<vector<long double>> negative_by_thread(
            thread_count, vector<long double>(Nat, 0.0L));
        vector<vector<long double>> absolute_by_thread(
            thread_count, vector<long double>(Nat, 0.0L));

        // --------------------------------------------------------------------
        // Main voxel loop
        // --------------------------------------------------------------------
#ifdef _OPENMP
#pragma omp parallel for collapse(2) schedule(static)
#endif
        for (int ix = 0; ix < reference.Nx; ++ix) {
            for (int iy = 0; iy < reference.Ny; ++iy) {
                for (int iz = 0; iz < reference.Nz; ++iz) {
                    const size_t idx = static_cast<size_t>(iz)
                                     + static_cast<size_t>(reference.Nz)
                                     * (static_cast<size_t>(iy)
                                     + static_cast<size_t>(reference.Ny)
                                     * static_cast<size_t>(ix));

                    const array<double, 3> P = node_pos(reference, ix, iy, iz);

                    int best_atom = -1;
                    double best_metric = numeric_limits<double>::infinity();
                    double best_distance = numeric_limits<double>::infinity();

                    // Find the Voronoi owner of this voxel.
                    for (int a = 0; a < Nat; ++a) {
                        const array<double, 3> displacement = {
                            P[0] - atom_positions[a][0],
                            P[1] - atom_positions[a][1],
                            P[2] - atom_positions[a][2]
                        };

                        const double distance = pbc_distance(pbc, displacement);
                        const double metric = (opt.weighted == "vdw")
                            ? distance / r_vdw[a]
                            : distance;

                        if (metric < best_metric) {
                            best_metric = metric;
                            best_atom = a;
                            best_distance = distance;
                        }
                    }

                    // The cutoff is applied after selecting the Voronoi owner.
                    const bool assigned = (
                        best_atom >= 0 && best_distance <= r_cut[best_atom]
                    );

                    if (assigned) {
                        const long double contribution =
                            static_cast<long double>(field[idx])
                            * static_cast<long double>(dV);

                        int tid = 0;
#ifdef _OPENMP
                        tid = omp_get_thread_num();
#endif

                        total_by_thread[tid][best_atom] += contribution;
                        absolute_by_thread[tid][best_atom] += fabsl(contribution);

                        if (contribution >= 0.0L) {
                            positive_by_thread[tid][best_atom] += contribution;
                        } else {
                            negative_by_thread[tid][best_atom] += contribution;
                        }
                    }

                    if (write_owner_map) {
                        // Thread-safe: each loop iteration writes a unique voxel.
                        owners[idx] = assigned ? static_cast<double>(best_atom + 1) : 0.0;
                    }
                }
            }
        }

        // Reduce thread-local integration buffers into the final per-atom values.
        for (int t = 0; t < thread_count; ++t) {
            for (int a = 0; a < Nat; ++a) {
                integral_total[a] += total_by_thread[t][a];
                integral_positive[a] += positive_by_thread[t][a];
                integral_negative[a] += negative_by_thread[t][a];
                integral_absolute[a] += absolute_by_thread[t][a];
            }
        }

        // CSV output to stdout.
        cout.setf(ios::fixed);
        cout << setprecision(10);
        cout << "idx,Z,x,y,z,Integral,Positive,Negative,AbsIntegral\n";

        for (int a = 0; a < Nat; ++a) {
            cout << a << ','
                 << reference.atoms[a].Z << ','
                 << reference.atoms[a].x << ','
                 << reference.atoms[a].y << ','
                 << reference.atoms[a].z << ','
                 << integral_total[a] << ','
                 << integral_positive[a] << ','
                 << integral_negative[a] << ','
                 << integral_absolute[a] << '\n';
        }

        // Optional CUBE map of the atom owning each assigned voxel.
        if (write_owner_map) {
            string description = string("Voronoi ")
                + (opt.weighted == "vdw" ? "[vdW-weighted]" : "[unweighted]")
                + ", cutoff_scale=" + to_string(opt.cutoff_scale) + " * r_vdw"
                + (opt.pbc ? " (PBC)" : "");

            write_cube(
                opt.owners_path,
                reference,
                owners,
                "VIPER-CT owner map (atom_index+1; 0=unassigned)",
                description
            );
        }

        return 0;

    } catch (const exception& e) {
        cerr << "VIPER-CT error: " << e.what() << '\n';
        return 2;
    }
}

# VIPER-CT

**VIPER-CT** — **V**oronoi **I**ntegration and **P**artitioning of **E**lectron **R**edistribution for **C**harge **T**ransfer — is a lightweight C++17 command-line program for integrating scalar fields stored in Gaussian CUBE files over atom-centered Voronoi regions.

The program was designed for real-space charge-transfer / electron-density-difference analysis. It supports both a precomputed single CUBE field and an on-the-fly three-CUBE density difference, optional van der Waals (vdW) weighting, a vdW-based spatial cutoff, triclinic periodic boundary conditions, OpenMP parallelization, and optional export of the resulting voxel-to-atom partition as a CUBE file.

## What VIPER-CT computes

For every grid voxel at position **P**, VIPER-CT assigns an atomic owner using either an ordinary Voronoi metric

```text
owner(P) = argmin_a |P - R_a|
```

or a multiplicatively vdW-weighted metric

```text
owner(P) = argmin_a |P - R_a| / r_vdW(a)
```

when `--weighted vdw` is enabled.

After the owner is identified, the voxel is included only when

```text
|P - R_owner| <= cutoff_scale * r_vdW(owner)
```

where the default `cutoff_scale` is `2.0`.

For every atom `a`, the program then evaluates the discrete real-space integral

```text
I_a = sum_{voxels owned by a} field(P) * dV
```

and separately accumulates its positive, negative, and absolute contributions.

> **Important:** `AbsIntegral` is `sum |field(P) * dV|` over the atom-owned voxels. It is therefore not the same quantity as `abs(Integral)`.

---

## Repository layout

```text
viper-ct/
├── README.md
├── Makefile
└── src/
    └── viper_ct.cpp
```

---

## Requirements

- A C++17 compiler (`g++` is recommended)
- OpenMP support for multithreaded execution
- Gaussian CUBE files with:
  - positive grid counts,
  - atomic coordinates and grid vectors expressed in Bohr,
  - one scalar value per grid point.

The code has no external library dependencies.

---

## Compilation

### Using the included Makefile

```bash
make
```

This creates

```text
./viper_ct
```

To remove the executable:

```bash
make clean
```

### Manual compilation

```bash
g++ -O3 -march=native -std=c++17 -fopenmp src/viper_ct.cpp -o viper_ct
```

If OpenMP is not available, the source can also be compiled without `-fopenmp`; the calculation will then run serially.

---

# Usage

```bash
./viper_ct [options] field.cube
```

or

```bash
./viper_ct [options] combined.cube fragment1.cube fragment2.cube
```

Display the built-in help with

```bash
./viper_ct --help
```

## Command-line options

| Option | Meaning |
|---|---|
| `--weighted none` | Ordinary Voronoi partition. This is the default. |
| `--weighted vdw` | Multiplicatively weighted Voronoi partition using `distance / r_vdW`. |
| `--cutoff-scale S` | Keep a voxel only if its distance from its owner is `<= S * r_vdW`. Default: `2.0`. |
| `--owners PATH` | Write a CUBE map containing the owner of each grid point. |
| `--pbc` | Use minimum-image periodic distances in the triclinic cell reconstructed from the CUBE grid. |
| `--threads N` | Set the number of OpenMP threads. `0` uses the OpenMP/runtime default. |
| `-h`, `--help` | Print the help message. |

---

# Workflow 1 — starting from a single CUBE

Use this mode when the scalar field you want to partition and integrate is already available as one CUBE file.

Typical examples include:

- a precomputed electron-density difference `delta_rho.cube`;
- a charge-density difference;
- another scalar field for which atom-resolved real-space integrals are meaningful.

The program does **not** modify the field in single-CUBE mode. It simply reads the grid values and integrates them over the atomic partitions.

### Minimal example

```bash
./viper_ct delta_rho.cube > per_atom.csv
```

### vdW-weighted partition

```bash
./viper_ct \
    --weighted vdw \
    --cutoff-scale 2.0 \
    delta_rho.cube \
    > per_atom_weighted.csv
```

### With PBC and 8 OpenMP threads

```bash
./viper_ct \
    --weighted vdw \
    --cutoff-scale 2.0 \
    --pbc \
    --threads 8 \
    delta_rho.cube \
    > per_atom_weighted_pbc.csv
```

### Also save the atom-owner map

```bash
./viper_ct \
    --weighted vdw \
    --cutoff-scale 2.0 \
    --owners owners.cube \
    delta_rho.cube \
    > per_atom.csv
```

In `owners.cube`:

```text
0        = voxel not assigned because it lies outside the cutoff
1        = atom with zero-based index 0
2        = atom with zero-based index 1
3        = atom with zero-based index 2
...
```

The owner map can be inspected in CUBE-compatible visualization software to verify the spatial partition.

---

# Workflow 2 — starting from three CUBE files

Use this mode when the density difference has not yet been generated explicitly.

VIPER-CT reads three scalar fields on the same real-space grid and constructs

```text
delta(r) = combined(r) - fragment1(r) - fragment2(r)
```

voxel by voxel before carrying out the Voronoi integration.

For an electron-density calculation this would usually correspond to

```text
Delta rho(r) = rho_AB(r) - rho_A(r) - rho_B(r)
```

provided that all three densities were generated consistently.

### Minimal example

```bash
./viper_ct \
    complex.cube \
    fragment_A.cube \
    fragment_B.cube \
    > per_atom_delta_rho.csv
```

### Recommended weighted example

```bash
./viper_ct \
    --weighted vdw \
    --cutoff-scale 2.0 \
    --threads 8 \
    complex.cube \
    fragment_A.cube \
    fragment_B.cube \
    > per_atom_delta_rho_weighted.csv
```

### Periodic system with owner-map output

```bash
./viper_ct \
    --weighted vdw \
    --cutoff-scale 2.0 \
    --pbc \
    --threads 8 \
    --owners owners.cube \
    complex.cube \
    fragment_A.cube \
    fragment_B.cube \
    > per_atom_delta_rho_weighted.csv
```

## Three-CUBE compatibility requirements

The three CUBE files must have the same:

- `Nx`, `Ny`, and `Nz` grid dimensions;
- grid origin;
- `dX`, `dY`, and `dZ` voxel vectors.

VIPER-CT checks these quantities before calculating the difference.

The atomic geometry used to define all Voronoi regions is taken **only from the first CUBE file** (`combined.cube`). The two fragment CUBE files may therefore contain different atom lists as long as their scalar fields are defined on the same grid.

For physically meaningful density-difference calculations, the three fields should also have been generated with mutually consistent electronic-structure settings, geometry, simulation cell, and real-space grid.

---

# Output

VIPER-CT writes a CSV table to standard output:

```text
idx,Z,x,y,z,Integral,Positive,Negative,AbsIntegral
```

The columns are:

| Column | Meaning |
|---|---|
| `idx` | Zero-based atom index in the first/reference CUBE. |
| `Z` | Atomic number. |
| `x,y,z` | Atomic coordinates as stored in the CUBE file. |
| `Integral` | Signed integral of the field over all assigned voxels belonging to that atom. |
| `Positive` | Sum of positive voxel contributions. |
| `Negative` | Sum of negative voxel contributions. |
| `AbsIntegral` | Sum of the absolute value of every voxel contribution. |

Because the CSV is printed to standard output, redirect it directly to a file:

```bash
./viper_ct ... > results.csv
```

Warnings and errors are printed to standard error, so they do not contaminate the CSV output.

## Units

VIPER-CT does not impose a physical unit on the scalar field. The integrated unit is simply

```text
[field unit] x [CUBE coordinate unit]^3
```

For example, if the CUBE contains an electron number density in `e / bohr^3`, the integrated quantities are in electrons (`e`).

The current implementation assumes CUBE coordinates/grid vectors are expressed in **Bohr**, because the internal vdW radii are converted from Angstrom to Bohr before distances are compared.

---

# Sign convention for charge-transfer analysis

In three-CUBE mode VIPER-CT always constructs

```text
combined - fragment1 - fragment2
```

If the CUBE values are **electron number densities**, then:

- positive values of `Delta rho` indicate electron accumulation;
- negative values indicate electron depletion.

The atomic `Integral` therefore measures the net integrated redistribution inside that atom's assigned spatial region.

If your CUBE contains a quantity with a different sign convention, interpret the resulting integrals accordingly. VIPER-CT performs the numerical operation above without changing signs.

---

# Available atomic vdW radii

The following vdW radii are currently hard-coded in `src/viper_ct.cpp`.

| Element | Symbol | Atomic number `Z` | vdW radius (A) |
|---|---:|---:|---:|
| Hydrogen | H | 1 | 1.20 |
| Carbon | C | 6 | 1.70 |
| Nitrogen | N | 7 | 1.55 |
| Oxygen | O | 8 | 1.52 |
| Fluorine | F | 9 | 1.47 |
| Sodium | Na | 11 | 2.27 |
| Magnesium | Mg | 12 | 1.73 |
| Silicon | Si | 14 | 2.10 |
| Phosphorus | P | 15 | 1.80 |
| Sulfur | S | 16 | 1.80 |
| Chlorine | Cl | 17 | 1.75 |
| Potassium | K | 19 | 2.75 |
| Calcium | Ca | 20 | 2.31 |
| Bromine | Br | 35 | 1.85 |
| Iodine | I | 53 | 1.98 |

For any atomic number that is not listed, VIPER-CT currently uses a fallback radius of **1.50 A** and prints a warning to standard error.

For reproducible published calculations, it is strongly recommended to add the intended radius explicitly instead of relying on the fallback.

---

# Adding a new atom / element

Open

```text
src/viper_ct.cpp
```

and locate the function

```cpp
static const unordered_map<int, double>& vdw_radii_angstrom()
```

The table has the form

```cpp
static const unordered_map<int, double> radii = {
    {1, 1.20},   // H
    {6, 1.70},   // C
    {7, 1.55},   // N
    ...
};
```

Add the new element using

```cpp
{atomic_number, radius_in_angstrom},  // Symbol
```

For example, schematically:

```cpp
{Z, radius_A},  // ElementSymbol
```

where:

- `Z` is the element's atomic number;
- `radius_A` is the vdW radius in Angstrom.

Do not convert the value to Bohr manually. VIPER-CT performs the Angstrom-to-Bohr conversion internally.

After editing the table, recompile:

```bash
make clean
make
```

The new radius is then used both for:

1. the normalized distance `distance / r_vdW` in `--weighted vdw` mode;
2. the cutoff radius `cutoff_scale * r_vdW` in both weighted and unweighted modes.

---

# Periodic boundary conditions

With

```bash
--pbc
```

VIPER-CT reconstructs the direct lattice vectors from the CUBE grid:

```text
A = Nx * dX
B = Ny * dY
C = Nz * dZ
```

and evaluates atom-voxel distances using a triclinic minimum-image convention.

Use this option when the scalar field represents a periodic simulation cell and atoms near opposite cell boundaries should compete for the same voxels through periodic images.

Do not enable `--pbc` for isolated/nonperiodic CUBE calculations unless that is explicitly the intended partition.

---

# Practical examples

## A precomputed density-difference CUBE for an isolated dimer

```bash
./viper_ct \
    --weighted vdw \
    --cutoff-scale 2.0 \
    --threads 8 \
    delta_rho.cube \
    > dimer_per_atom.csv
```

## Build the density difference from complex + two isolated fragments

```bash
./viper_ct \
    --weighted vdw \
    --cutoff-scale 2.0 \
    --threads 8 \
    dimer.cube monomer_A.cube monomer_B.cube \
    > dimer_per_atom.csv
```

## Periodic interface calculation

```bash
./viper_ct \
    --weighted vdw \
    --cutoff-scale 2.0 \
    --pbc \
    --threads 16 \
    --owners interface_owners.cube \
    interface.cube fragment_A.cube fragment_B.cube \
    > interface_per_atom.csv
```

---

# Numerical notes

- Grid integration uses the triclinic voxel volume

  ```text
  dV = |(dX x dY) . dZ|
  ```

- The main voxel loop is parallelized with OpenMP.
- Integration accumulators use `long double` to reduce summation noise.
- The owner map is written only when requested, avoiding the additional grid-sized memory allocation otherwise.
- In three-CUBE mode, the difference field is generated in memory; no intermediate `Delta rho` CUBE is required.
- Changing `--cutoff-scale`, vdW radii, grid spacing, or the weighting scheme changes the real-space partition and can therefore change the integrated atomic values. These settings should be reported alongside published results.

---

# Reproducibility checklist

For a calculation intended to accompany a publication, record at least:

```text
VIPER-CT source/version or commit
single-CUBE or three-CUBE workflow
weighted mode: none / vdw
cutoff scale
enabled/disabled PBC
vdW radii used for all elements
grid dimensions and spacing
electronic-structure method used to generate the CUBE files
sign convention of the scalar field
```

For three-CUBE calculations, also verify that all three fields were evaluated on exactly the same grid.

---

# License and citation

Before public release, add the license you want to distribute VIPER-CT under and the final citation for the accompanying paper.

A typical repository can then include a `LICENSE` file and a `CITATION.cff` file alongside this README.

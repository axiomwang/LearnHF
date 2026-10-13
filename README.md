# LearnHF

LearnHF is a compact C++23 Hartree–Fock program for molecular electronic-structure calculations. It is written as a readable research and learning code: the implementation exposes the main numerical steps of a self-consistent-field (SCF) calculation instead of hiding them behind a large chemistry framework.

The executable is named `hftest`. It reads a TOML input file, loads a Gaussian basis set in Basis Set Exchange (BSE) Gaussian format, and runs restricted or unrestricted Hartree–Fock (RHF/UHF).

## Capabilities

- **RHF and UHF SCF.** The driver selects RHF for a singlet with an even electron count; other valid spin/electron combinations use UHF. UHF keeps separate alpha and beta densities and reports the exact and computed $\langle \hat S^2\rangle$ values.
- **Gaussian-basis integrals.** The McMurchie–Davidson engine evaluates overlap, kinetic, nuclear-attraction, and electron-repulsion integrals over Cartesian Gaussian functions. Boys functions use small-$x$ series and large-$x$ recurrences with an asymptotic fallback.
- **Direct SCF.** Electron-repulsion integrals are recomputed during each SCF iteration rather than stored on disk. The J/K build uses ERI permutational symmetry plus shell- and function-level Schwarz screening.
- **SCF numerics.** The solver uses a core-Hamiltonian initial guess, canonical orthogonalization of the overlap matrix (discarding eigenvalues below $10^{-8}$), and Pulay DIIS. Convergence requires both the energy change and RMS density change to be below the configured threshold.
- **Basis-set input.** The parser accepts BSE Gaussian-format files, splits `SP` shells into separate S and P shells, handles Fortran `D` exponents, normalizes primitive contributions on load, and applies shell normalization during integral construction. The repository includes Pople, Dunning correlation-consistent, def2, and LANL-family files.
- **LibTorch linear algebra.** Dense matrices, eigendecompositions, and the DIIS linear solve use `torch::Tensor`/LibTorch.

The implementation is intentionally small: `src/main.cpp` handles input and dispatch, `src/basis_parser.cpp` loads basis data, `src/boys.cpp` implements Boys functions, `src/integrals.cpp` builds AO integrals and J/K matrices, and `src/scf.cpp` contains the RHF/UHF solvers. Public data structures and solver APIs are in `include/`.

## Build

### Requirements

- CMake **3.20** or newer (the minimum in `CMakeLists.txt`)
- A compiler with C++23 support (the existing usage docs list AppleClang 17+, GCC 13+, and Clang 17+ as supported targets)
- LibTorch 2.0 or newer; standalone LibTorch is preferred. CMake can also use a compatible Python `torch` installation when its CMake package is discoverable. The local check below used standalone LibTorch 2.11.0.
- A C++23-capable build tool such as Make or Ninja

The TOML parser is vendored as `external/toml.hpp`; no separate TOML installation is required.

```bash
git clone https://github.com/axiomwang/LearnHF.git
cd LearnHF

cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
# If LibTorch is not found automatically, configure with:
# cmake -S . -B build -DCMAKE_BUILD_TYPE=Release \
#   -DTORCH_HINT=/path/to/libtorch/share/cmake/Torch
cmake --build build --parallel
```

The executable is written to `build/hftest`.

## Running a calculation

Run the bundled water example from the repository root:

```bash
./build/hftest input.toml
```

Input coordinates are Cartesian **angstroms** and are converted to bohr internally. The minimal input sections are:

```toml
[system]
charge       = 0
multiplicity = 1       # 1 with an even electron count -> RHF; otherwise UHF
basis        = "sto-3g"

[[atoms]]
symbol = "H"
x = 0.0
y = 0.0
z = 0.0

[scf]
max_iter    = 100
convergence = 1e-8
diis_size   = 8
screen_tol  = 1e-10
```

`basis` is matched case-insensitively to a filename in `basis/`. The program searches for that file relative to the input file, beside the executable, and in `./basis`. See [`docs/usage.rst`](docs/usage.rst) for the complete key reference and open-shell example.

## Reproducibility checks

The following smoke tests were run locally on macOS with AppleClang 21, CMake 4.3.3, and the standalone LibTorch installation found at `/usr/local/share/cmake/Torch`:

| Input | Observed result |
| --- | --- |
| Bundled `input.toml` (H₂O/STO-3G, RHF) | Converged in 12 iterations; total energy `-75.0171735279 Ha` |
| H atom, `multiplicity = 2`, STO-3G (UHF) | Converged in 2 iterations; total energy `-0.4665818504 Ha`; $\langle \hat S^2\rangle = 0.7500$ |

These are local smoke-test results, not a cross-platform benchmark or a validation against a production quantum-chemistry package. The bundled input reproduces the total energy printed in the existing Sphinx documentation, but the geometry in `docs/usage.rst` differs from the committed input and `docs/index.rst` reports different HOMO/LUMO values; those older values are not treated as a reference here.

## Scope and limitations

- ECP blocks in basis files are skipped; effective-core-potential integrals and the associated core potentials are **not implemented**. LANL files therefore provide only their parsed valence Gaussian data.
- This is a molecular, single-determinant Hartree–Fock implementation under the Born–Oppenheimer approximation. It does not implement periodic boundary conditions or solid-state calculations.
- The code performs direct SCF and has no disk-backed ERI storage or density-fitting backend, geometry optimization, excited-state method, correlation method, or analytic gradients.
- The executable uses the default LibTorch tensor device and has no explicit GPU/device-selection path.
- There is no automated unit-test suite or CI configuration in the repository. Build output currently includes a few compiler warnings from `src/integrals.cpp`; the smoke tests above completed successfully.
- The included documentation is the detailed technical reference. Start with [`docs/index.rst`](docs/index.rst), then see [`docs/principles.rst`](docs/principles.rst), [`docs/algorithms.rst`](docs/algorithms.rst), [`docs/usage.rst`](docs/usage.rst), and the API pages under [`docs/api/`](docs/api/).

## License

No `LICENSE` file is currently present in the repository. Check with the project owner before redistributing the code.

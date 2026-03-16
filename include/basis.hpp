#pragma once
#include "molecule.hpp"
#include <vector>
#include <array>
#include <string>
#include <filesystem>

// One Gaussian primitive: exponent α, contraction coefficient (for S),
// and optional p-coefficient (used only when this primitive came from an SP shell).
struct Primitive {
    double exp;       // exponent α
    double coeff;     // contraction coefficient (already normalised × prim_norm)
};

// A contracted Gaussian shell:
//   χ_μ(r) = N_μ * Σ_k c_k * g(α_k, l, r-A)
// All primitives share the same angular momentum and center.
struct Shell {
    int         atom_idx;               // index into Molecule::atoms
    std::array<double,3> center;        // Bohr
    int         l;                      // angular momentum (0=S,1=P,2=D,3=F,4=G)
    std::vector<Primitive> prims;       // exponents + coeffs (normalized)

    // Number of Cartesian functions in this shell
    int n_cart() const { return (l+1)*(l+2)/2; }
};

// A complete basis set for a molecule
struct BasisSet {
    std::vector<Shell> shells;

    int n_basis() const {
        int n = 0;
        for (auto& sh : shells) n += sh.n_cart();
        return n;
    }

    // Index of the first basis function belonging to shell i
    int shell_start(int i) const {
        int s = 0;
        for (int j = 0; j < i; ++j) s += shells[j].n_cart();
        return s;
    }
};

// Cartesian function ordering within a shell of angular momentum l.
// Returns list of (nx, ny, nz) tuples with nx+ny+nz == l.
// Ordering: iterate over nx descending, for each nx iterate ny descending,
// nz = l - nx - ny.  This matches the convention used by Gaussian09 / Psi4.
inline std::vector<std::array<int,3>> cart_components(int l) {
    std::vector<std::array<int,3>> comps;
    for (int nx = l; nx >= 0; --nx)
        for (int ny = l-nx; ny >= 0; --ny)
            comps.push_back({nx, ny, l-nx-ny});
    return comps;
}

// Load a basis set from the given BSE-format file for all atoms in the molecule.
// SP shells are split into separate S and P shells.
// ECP blocks are skipped with a warning printed to stderr.
BasisSet load_basis(const Molecule& mol, const std::filesystem::path& basis_file);

#pragma once
#include "molecule.hpp"
#include "basis.hpp"
#include <torch/torch.h>
#include <vector>
#include <string>

// SCF settings parsed from TOML [scf] section
struct SCFSettings {
    int    max_iter    = 100;
    double convergence = 1e-8;
    int    diis_size   = 8;
    double screen_tol  = 1e-10;
};

// ── Results ──────────────────────────────────────────────────────────────────

struct RHFResult {
    bool   converged    = false;
    int    n_iter       = 0;
    double E_elec       = 0.0;
    double E_nuc        = 0.0;
    double E_total      = 0.0;
    torch::Tensor C;            // MO coefficient matrix [nbas × nbas]
    torch::Tensor eps;          // orbital energies [nbas]
    torch::Tensor P;            // density matrix [nbas × nbas]
    torch::Tensor F;            // final Fock matrix [nbas × nbas]
};

struct UHFResult {
    bool   converged    = false;
    int    n_iter       = 0;
    double E_elec       = 0.0;
    double E_nuc        = 0.0;
    double E_total      = 0.0;
    torch::Tensor Ca;           // alpha MO coefficients
    torch::Tensor Cb;           // beta  MO coefficients
    torch::Tensor eps_a;        // alpha orbital energies
    torch::Tensor eps_b;        // beta  orbital energies
    torch::Tensor Pa;           // alpha density matrix
    torch::Tensor Pb;           // beta  density matrix
    double S2_exact     = 0.0;  // <S²> exact = S(S+1) where S=(na-nb)/2
    double S2_uhf       = 0.0;  // <S²> UHF (with spin contamination)
};

// ── Solvers ─────────────────────────────────────────────────────────────────

RHFResult run_rhf(const Molecule& mol, const BasisSet& bs,
                  const SCFSettings& settings);

UHFResult run_uhf(const Molecule& mol, const BasisSet& bs,
                  const SCFSettings& settings);

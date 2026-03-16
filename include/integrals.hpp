#pragma once
#include "basis.hpp"
#include "molecule.hpp"
#include <torch/torch.h>
#include <ATen/ops/linalg_eigh.h>
#include <ATen/ops/linalg_solve.h>

// ── 1-electron integral matrices ────────────────────────────────────────────

// Overlap matrix S_μν = <μ|ν>
torch::Tensor build_overlap(const BasisSet& bs);

// Kinetic energy matrix T_μν = <μ|-½∇²|ν>
torch::Tensor build_kinetic(const BasisSet& bs);

// Nuclear attraction matrix V_μν = Σ_A <μ| -Z_A / |r-R_A| |ν>
torch::Tensor build_nuclear(const BasisSet& bs, const Molecule& mol);

// Core Hamiltonian
inline torch::Tensor build_hcore(const BasisSet& bs, const Molecule& mol) {
    return build_kinetic(bs) + build_nuclear(bs, mol);
}

// ── 2-electron integrals (Direct SCF) ───────────────────────────────────────

// Compute the Schwarz screening matrix Q_μν = sqrt((μν|μν))
// Used to prescreen shell quartets: skip if Q_μν * Q_λσ < screen_tol.
torch::Tensor build_schwarz(const BasisSet& bs);

// Build the Coulomb matrix J_μν = Σ_{λσ} P_λσ (μν|λσ)
// and exchange matrix K_μν = Σ_{λσ} P_λσ (μσ|λν) [for ONE density matrix P]
// Returns {J, K}.
// screen is the Schwarz matrix; quartets with Q_μν*Q_λσ < tol are skipped.
std::pair<torch::Tensor, torch::Tensor>
build_JK(const BasisSet& bs, const torch::Tensor& P,
         const torch::Tensor& screen, double tol = 1e-10);

// RHF and UHF SCF solver with DIIS convergence acceleration
// Reference: Szabo & Ostlund, "Modern Quantum Chemistry" (1989)
//            Pulay, Chem. Phys. Lett. 73 (1980) 393.

#include "scf.hpp"
#include "integrals.hpp"
#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <stdexcept>
#include <ATen/ops/linalg_eigh.h>
#include <ATen/ops/linalg_solve.h>

// ═══════════════════════════════════════════════════════════════════════════
//  Utilities
// ═══════════════════════════════════════════════════════════════════════════

// Canonical orthogonalization: S = U Λ U^T → X = U Λ^{-1/2}
// Discards linear-dependent vectors (eigenvalue < thresh).
// Returns X [n × n_ind] where n_ind ≤ n.
static std::pair<torch::Tensor, int>
canonical_orthogonalization(const torch::Tensor& S, double thresh = 1e-8) {
    // eigh returns (eigenvalues, eigenvectors) for real symmetric matrix
        auto [vals, vecs] = at::linalg_eigh(S, "L");  // vals ascending
    auto va = vals.accessor<double,1>();
    int n = (int)vals.size(0);
    int n_ind = 0;
    for (int i = 0; i < n; ++i)
        if (va[i] >= thresh) ++n_ind;

    if (n_ind == 0)
        throw std::runtime_error("All overlap eigenvalues below threshold — linear dependency?");

    // Keep only the n_ind largest eigenvalues (they come last in ascending order)
    int skip = n - n_ind;
    auto vals_keep = vals.slice(0, skip, n);    // [n_ind]
    auto vecs_keep = vecs.slice(1, skip, n);    // [n × n_ind]

    auto X = vecs_keep * vals_keep.pow(-0.5).unsqueeze(0);  // [n × n_ind]
    return {X, n_ind};
}

// RMS of a matrix
static double rms(const torch::Tensor& M) {
    return M.pow(2).mean().sqrt().item<double>();
}

// ═══════════════════════════════════════════════════════════════════════════
//  DIIS extrapolation (Pulay)
//  Error vector: e = F·P·S - S·P·F  (in AO basis)
// ═══════════════════════════════════════════════════════════════════════════

struct DIIS {
    int max_size;
    std::vector<torch::Tensor> F_list;
    std::vector<torch::Tensor> e_list;

    explicit DIIS(int sz) : max_size(sz) {}

    // Push new Fock and error vector
    void push(const torch::Tensor& F, const torch::Tensor& e) {
        F_list.push_back(F.clone());
        e_list.push_back(e.clone());
        if ((int)F_list.size() > max_size) {
            F_list.erase(F_list.begin());
            e_list.erase(e_list.begin());
        }
    }

    // Return DIIS-extrapolated Fock matrix (or F itself if not enough vectors)
    torch::Tensor extrapolate() {
        int m = (int)F_list.size();
        if (m < 2) return F_list.back().clone();

        // Build B matrix: B_ij = <e_i | e_j>
        // B is (m+1)×(m+1) with last row/col being Lagrange constraint
        auto B = torch::zeros({m+1, m+1}, torch::kDouble);
        auto Ba = B.accessor<double,2>();
        for (int i = 0; i < m; ++i)
            for (int j = 0; j < m; ++j)
                Ba[i][j] = (e_list[i] * e_list[j]).sum().item<double>();
        for (int i = 0; i < m; ++i) { Ba[m][i] = -1.0; Ba[i][m] = -1.0; }

        auto rhs = torch::zeros({m+1}, torch::kDouble);
        rhs[m] = -1.0;

        // Solve B·c = rhs
        torch::Tensor c;
        try {
              c = at::linalg_solve(B, rhs.unsqueeze(1)).squeeze(1);
        } catch (...) {
            // If singular (e.g., duplicate error vectors), fall back to last F
            return F_list.back().clone();
        }

        auto ca = c.accessor<double,1>();
        auto F_diis = torch::zeros_like(F_list[0]);
        for (int i = 0; i < m; ++i)
            F_diis += ca[i] * F_list[i];
        return F_diis;
    }
};

// ═══════════════════════════════════════════════════════════════════════════
//  RHF
// ═══════════════════════════════════════════════════════════════════════════

RHFResult run_rhf(const Molecule& mol, const BasisSet& bs,
                  const SCFSettings& cfg) {
    int n     = bs.n_basis();
    int nocc  = mol.n_alpha();   // = n_beta for RHF (closed-shell check below)

    if (mol.n_alpha() != mol.n_beta())
        throw std::runtime_error("RHF requires a closed-shell molecule (n_alpha == n_beta). "
                                 "Use UHF for open-shell systems.");
    if (nocc > n)
        throw std::runtime_error("More occupied orbitals than basis functions.");

    std::cout << "\n=== RHF / " << bs.n_basis() << " basis functions, "
              << nocc << " occ. orbitals ===\n";

    // Build 1e- integrals
    auto S     = build_overlap(bs);
    auto Hcore = build_hcore(bs, mol);
    auto Sch   = build_schwarz(bs);

    // Canonical orthogonalization
    auto [X, n_ind] = canonical_orthogonalization(S);
    std::cout << "  Orthogonal basis dimension: " << n_ind << "/" << n << "\n";

    // Initial guess: F0 = Hcore
    torch::Tensor F = Hcore.clone();
    torch::Tensor P = torch::zeros({n, n}, torch::kDouble);
    torch::Tensor C = torch::zeros({n, n_ind}, torch::kDouble);

    auto diis = DIIS(cfg.diis_size);
    double E_prev = 0.0;
    RHFResult result;
    result.E_nuc = mol.nuclear_repulsion();

    std::cout << "  Iter     E_elec          ΔE            RMS(ΔP)\n";
    std::cout << "  ─────────────────────────────────────────────────\n";

    for (int iter = 1; iter <= cfg.max_iter; ++iter) {
        // Transform Fock to orthogonal basis
        auto Fp = torch::mm(X.transpose(0,1), torch::mm(F, X));

        // Diagonalize
            auto [eps_p, Cp] = at::linalg_eigh(Fp, "L");

        // Back-transform
        C = torch::mm(X, Cp);                          // [n × n_ind]

        // Form density matrix  P = 2 C_occ C_occ^T
        auto C_occ = C.slice(1, 0, nocc);              // [n × nocc]
        auto P_new = 2.0 * torch::mm(C_occ, C_occ.transpose(0,1));

        // Build 2e- part
        auto [J, K] = build_JK(bs, P_new, Sch, cfg.screen_tol);
        auto G  = J - 0.5 * K;
        F = Hcore + G;

        // Electronic energy
        double E_elec = 0.5 * torch::mm(P_new, Hcore + F)
                                  .diagonal().sum().item<double>();

        // DIIS error vector: e = F P S - S P F
        auto FPS  = torch::mm(F, torch::mm(P_new, S));
        auto SPF  = torch::mm(S, torch::mm(P_new, F));
        auto err  = FPS - SPF;
        diis.push(F, err);
        F = diis.extrapolate();

        double dE     = E_elec - E_prev;
        double dP_rms = rms(P_new - P);

        std::cout << "  " << std::setw(4) << iter
                  << "  " << std::fixed << std::setprecision(10) << std::setw(16) << E_elec
                  << "  " << std::setw(12) << std::scientific << std::setprecision(4) << dE
                  << "  " << std::setw(12) << dP_rms << "\n";

        P      = P_new;
        E_prev = E_elec;

        if (iter > 1 && std::abs(dE) < cfg.convergence && dP_rms < cfg.convergence) {
            result.converged = true;
            result.n_iter    = iter;
            result.E_elec    = E_elec;
            result.E_total   = E_elec + result.E_nuc;
            result.C         = C;
            result.eps       = eps_p.slice(0, 0, n_ind);  // orbital energies
            result.P         = P;
            result.F         = F;
            break;
        }
    }

    if (!result.converged) {
        std::cerr << "  WARNING: RHF did not converge after "
                  << cfg.max_iter << " iterations.\n";
        result.n_iter  = cfg.max_iter;
        result.E_elec  = E_prev;
        result.E_total = E_prev + result.E_nuc;
        result.C       = C;
        result.P       = P;
        result.F       = F;
    }

    return result;
}

// ═══════════════════════════════════════════════════════════════════════════
//  UHF
// ═══════════════════════════════════════════════════════════════════════════

UHFResult run_uhf(const Molecule& mol, const BasisSet& bs,
                  const SCFSettings& cfg) {
    int n    = bs.n_basis();
    int na   = mol.n_alpha();
    int nb   = mol.n_beta();

    if (na < nb) throw std::runtime_error("n_alpha < n_beta — check multiplicity definition.");
    if (na + nb != mol.n_electrons())
        throw std::runtime_error("n_alpha + n_beta != n_electrons.");

    std::cout << "\n=== UHF / " << n << " basis functions, α=" << na << " β=" << nb << " ===\n";

    auto S     = build_overlap(bs);
    auto Hcore = build_hcore(bs, mol);
    auto Sch   = build_schwarz(bs);

    auto [X, n_ind] = canonical_orthogonalization(S);
    std::cout << "  Orthogonal basis dimension: " << n_ind << "/" << n << "\n";

    // Initial density matrices: guess from H_core
    auto F0     = Hcore.clone();
    auto Fpr    = torch::mm(X.transpose(0,1), torch::mm(F0, X));
        auto [e0, C0p] = at::linalg_eigh(Fpr, "L");
    auto C0     = torch::mm(X, C0p);

    torch::Tensor Ca = C0.clone();
    torch::Tensor Cb = C0.clone();
    torch::Tensor Pa = torch::zeros({n, n}, torch::kDouble);
    torch::Tensor Pb = torch::zeros({n, n}, torch::kDouble);
    torch::Tensor eps_a, eps_b;

    if (na > 0) {
        auto Ca_occ = Ca.slice(1, 0, na);
        Pa = torch::mm(Ca_occ, Ca_occ.transpose(0,1));
    }
    if (nb > 0) {
        auto Cb_occ = Cb.slice(1, 0, nb);
        Pb = torch::mm(Cb_occ, Cb_occ.transpose(0,1));
    }

    auto diis_a = DIIS(cfg.diis_size);
    auto diis_b = DIIS(cfg.diis_size);

    double E_prev = 0.0;
    UHFResult result;
    result.E_nuc    = mol.nuclear_repulsion();
    double S_val    = 0.5 * (mol.n_alpha() - mol.n_beta());
    result.S2_exact = S_val * (S_val + 1.0);

    std::cout << "  Iter     E_elec          ΔE            RMS(ΔPα+Pβ)\n";
    std::cout << "  ──────────────────────────────────────────────────────\n";

    for (int iter = 1; iter <= cfg.max_iter; ++iter) {
        auto Pt = Pa + Pb;  // total density

        // J is the same for both spins (uses total density)
        auto [J, Ka] = build_JK(bs, Pa, Sch, cfg.screen_tol);
        auto [Jb, Kb] = build_JK(bs, Pb, Sch, cfg.screen_tol);
        auto Jtot = J + Jb;  // Coulomb from total density

        auto Fa_new = Hcore + Jtot - Ka;
        auto Fb_new = Hcore + Jtot - Kb;

        // DIIS error vectors
        auto err_a = torch::mm(Fa_new, torch::mm(Pa, S)) - torch::mm(S, torch::mm(Pa, Fa_new));
        auto err_b = torch::mm(Fb_new, torch::mm(Pb, S)) - torch::mm(S, torch::mm(Pb, Fb_new));
        diis_a.push(Fa_new, err_a);
        diis_b.push(Fb_new, err_b);
        auto Fa = diis_a.extrapolate();
        auto Fb = diis_b.extrapolate();

        // Diagonalize α
        auto Fap = torch::mm(X.transpose(0,1), torch::mm(Fa, X));
            auto [eap, Cap] = at::linalg_eigh(Fap, "L");
        Ca     = torch::mm(X, Cap);
        eps_a  = eap;
        auto Ca_occ = Ca.slice(1, 0, na);
        auto Pa_new = (na > 0) ? torch::mm(Ca_occ, Ca_occ.transpose(0,1))
                               : torch::zeros({n,n}, torch::kDouble);

        // Diagonalize β
            auto Fbp = torch::mm(X.transpose(0,1), torch::mm(Fb, X));
            auto [ebp, Cbp] = at::linalg_eigh(Fbp, "L");
            Cb     = torch::mm(X, Cbp);
        eps_b  = ebp;
        auto Cb_occ = Cb.slice(1, 0, nb);
        auto Pb_new = (nb > 0) ? torch::mm(Cb_occ, Cb_occ.transpose(0,1))
                               : torch::zeros({n,n}, torch::kDouble);

        // UHF electronic energy: E = Tr[(Pα+Pβ)·H + Pα·Fα + Pβ·Fβ] / 2
        // using the Fock matrices *before* DIIS (Fa_new, Fb_new)
        double E_elec = 0.5 * (
              torch::mm(Pa_new + Pb_new, Hcore).diagonal().sum().item<double>()
            + torch::mm(Pa_new, Fa_new).diagonal().sum().item<double>()
            + torch::mm(Pb_new, Fb_new).diagonal().sum().item<double>()
        );

        double dE     = E_elec - E_prev;
        double dP_rms = rms(Pa_new - Pa) + rms(Pb_new - Pb);

        std::cout << "  " << std::setw(4) << iter
                  << "  " << std::fixed << std::setprecision(10) << std::setw(16) << E_elec
                  << "  " << std::setw(12) << std::scientific << std::setprecision(4) << dE
                  << "  " << std::setw(12) << dP_rms << "\n";

        Pa     = Pa_new;
        Pb     = Pb_new;
        E_prev = E_elec;

        if (iter > 1 && std::abs(dE) < cfg.convergence && dP_rms < cfg.convergence) {
            result.converged = true;
            result.n_iter    = iter;
            result.E_elec    = E_elec;
            result.E_total   = E_elec + result.E_nuc;
            result.Ca        = Ca;
            result.Cb        = Cb;
            result.eps_a     = eps_a;
            result.eps_b     = eps_b;
            result.Pa        = Pa;
            result.Pb        = Pb;
            break;
        }
    }

    if (!result.converged) {
        std::cerr << "  WARNING: UHF did not converge after "
                  << cfg.max_iter << " iterations.\n";
        result.n_iter  = cfg.max_iter;
        result.E_elec  = E_prev;
        result.E_total = E_prev + result.E_nuc;
        result.Ca      = Ca;
        result.Cb      = Cb;
        result.eps_a   = eps_a;
        result.eps_b   = eps_b;
        result.Pa      = Pa;
        result.Pb      = Pb;
    }

    // Compute <S²> = S(S+1) + n_beta - Σ_{ij} |<φ_i^α|φ_j^β>|²
    // where the sum is over all occupied α and β orbitals.
    // <φ_i^α|φ_j^β> = (C_i^α)^T · S · C_j^β
    if (result.converged || true) {
        double S_exact = 0.5 * (na - nb);
        double S2_exact = S_exact * (S_exact + 1.0);
        // Overlap between alpha and beta occupied MOs: M = Ca_occ^T S Cb_occ
        double S2_cont = 0.0;
        if (na > 0 && nb > 0) {
            auto Ca_occ = result.Ca.slice(1, 0, na);
            auto Cb_occ = result.Cb.slice(1, 0, nb);
            auto M = torch::mm(Ca_occ.transpose(0,1), torch::mm(S, Cb_occ)); // [na × nb]
            S2_cont = torch::mm(M, M.transpose(0,1)).diagonal().sum().item<double>();
        }
        result.S2_exact = S2_exact;
        result.S2_uhf   = S2_exact + nb - S2_cont;
    }

    return result;
}

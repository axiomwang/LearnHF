#include "output.hpp"
#include <iostream>
#include <iomanip>
#include <cmath>

static constexpr double HARTREE_TO_EV = 27.211386245988;

void print_molecule(const Molecule& mol) {
    std::cout << "\n╔══════════════════════════════════════════════════════════╗\n";
    std::cout << "║                   Molecule Information                  ║\n";
    std::cout << "╚══════════════════════════════════════════════════════════╝\n";
    std::cout << "  Charge: " << mol.charge
              << "  Multiplicity: " << mol.multiplicity
              << "  Electrons: " << mol.n_electrons() << "\n";
    std::cout << "  α electrons: " << mol.n_alpha()
              << "  β electrons: " << mol.n_beta() << "\n\n";
    std::cout << "  Atom    X (Bohr)       Y (Bohr)       Z (Bohr)\n";
    std::cout << "  ────────────────────────────────────────────────\n";
    for (auto& a : mol.atoms) {
        std::cout << "  " << std::left << std::setw(4) << a.symbol
                  << "  " << std::right << std::fixed << std::setprecision(6)
                  << std::setw(12) << a.coords[0]
                  << "  " << std::setw(12) << a.coords[1]
                  << "  " << std::setw(12) << a.coords[2] << "\n";
    }
}

void print_basis_info(const BasisSet& bs, const Molecule& mol) {
    std::cout << "\n  Basis functions: " << bs.n_basis()
              << "  Shells: " << bs.shells.size() << "\n";
}

static void print_mo_energies(const torch::Tensor& eps, const std::string& label,
                               int nocc, int n_print = 20) {
    int n = (int)eps.size(0);
    int shown = std::min(n, n_print);
    std::cout << "\n  " << label << " orbital energies (Hartree | eV):\n";
    std::cout << "  ──────────────────────────────────────────────────────────\n";
    for (int i = 0; i < shown; ++i) {
        double e = eps[i].item<double>();
        std::cout << "  " << (i < nocc ? "occ" : "vir")
                  << " " << std::setw(3) << i+1
                  << "  " << std::fixed << std::setprecision(6) << std::setw(12) << e
                  << " Ha  " << std::setw(10) << e * HARTREE_TO_EV << " eV\n";
    }
    if (shown < n)
        std::cout << "  ... (" << n - shown << " more orbitals not shown)\n";
}

void print_rhf_result(const RHFResult& res, const Molecule& mol, const BasisSet& bs) {
    std::cout << "\n╔══════════════════════════════════════════════════════════╗\n";
    std::cout << "║                    RHF Results                          ║\n";
    std::cout << "╚══════════════════════════════════════════════════════════╝\n";
    std::cout << "  Converged: " << (res.converged ? "YES" : "NO")
              << "  Iterations: " << res.n_iter << "\n\n";
    std::cout << std::fixed << std::setprecision(10);
    std::cout << "  Nuclear repulsion energy:  " << std::setw(18) << res.E_nuc  << " Ha\n";
    std::cout << "  Electronic energy:         " << std::setw(18) << res.E_elec << " Ha\n";
    std::cout << "  ─────────────────────────────────────────────────────────\n";
    std::cout << "  Total energy:              " << std::setw(18) << res.E_total << " Ha\n";
    std::cout << "                             " << std::setw(18)
              << res.E_total * HARTREE_TO_EV << " eV\n";

    int nocc = mol.n_alpha();
    int n    = (int)res.eps.size(0);
    if (nocc > 0 && nocc <= n) {
        double homo = res.eps[nocc-1].item<double>();
        std::cout << "\n  HOMO energy:               " << std::setw(18) << homo
                  << " Ha  (" << homo * HARTREE_TO_EV << " eV)\n";
        if (nocc < n) {
            double lumo = res.eps[nocc].item<double>();
            std::cout << "  LUMO energy:               " << std::setw(18) << lumo
                      << " Ha  (" << lumo * HARTREE_TO_EV << " eV)\n";
            std::cout << "  HOMO-LUMO gap:             " << std::setw(18) << lumo - homo
                      << " Ha  (" << (lumo-homo) * HARTREE_TO_EV << " eV)\n";
        }
    }

    print_mo_energies(res.eps, "RHF", nocc);
}

void print_uhf_result(const UHFResult& res, const Molecule& mol, const BasisSet& bs) {
    std::cout << "\n╔══════════════════════════════════════════════════════════╗\n";
    std::cout << "║                    UHF Results                          ║\n";
    std::cout << "╚══════════════════════════════════════════════════════════╝\n";
    std::cout << "  Converged: " << (res.converged ? "YES" : "NO")
              << "  Iterations: " << res.n_iter << "\n\n";
    std::cout << std::fixed << std::setprecision(10);
    std::cout << "  Nuclear repulsion energy:  " << std::setw(18) << res.E_nuc  << " Ha\n";
    std::cout << "  Electronic energy:         " << std::setw(18) << res.E_elec << " Ha\n";
    std::cout << "  ─────────────────────────────────────────────────────────\n";
    std::cout << "  Total energy:              " << std::setw(18) << res.E_total << " Ha\n";
    std::cout << "                             " << std::setw(18)
              << res.E_total * HARTREE_TO_EV << " eV\n";

    std::cout << "\n  Spin: <S²> = " << std::fixed << std::setprecision(4)
              << res.S2_uhf << "  (exact = " << res.S2_exact << ")";
    double contamination = res.S2_uhf - res.S2_exact;
    if (contamination > 0.05)
        std::cout << "  ← spin contaminated!";
    std::cout << "\n";

    int na = mol.n_alpha(), nb = mol.n_beta();
    print_mo_energies(res.eps_a, "α (UHF)", na);
    print_mo_energies(res.eps_b, "β (UHF)", nb);
}

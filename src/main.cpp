#define TOML_HEADER_ONLY 1
#include <toml.hpp>
#include "molecule.hpp"
#include "basis.hpp"
#include "scf.hpp"
#include "output.hpp"
#include <iostream>
#include <filesystem>
#include <stdexcept>
#include <string>
#include <torch/torch.h>

namespace fs = std::filesystem;

// ── Locate the basis set file ─────────────────────────────────────────────
// Search: 1) relative to input file, 2) <executable_dir>/basis/, 3) ./basis/
static fs::path find_basis_file(const std::string& basis_name,
                                const fs::path& input_dir,
                                const fs::path& exe_dir) {
    // Lowercase the name for file lookup
    std::string name_lower = basis_name;
    for (auto& c : name_lower) c = (char)std::tolower((unsigned char)c);

    std::vector<fs::path> candidates = {
        input_dir / "basis" / name_lower,
        exe_dir   / "basis" / name_lower,
        fs::path("basis") / name_lower,
        input_dir / name_lower,
    };
    for (auto& c : candidates)
        if (fs::exists(c)) return c;

    throw std::runtime_error("Cannot find basis set file for \"" + basis_name +
                             "\". Searched: " + (input_dir/"basis"/name_lower).string());
}

// ── Parse TOML molecule description ──────────────────────────────────────
static Molecule parse_molecule(const toml::table& tbl) {
    Molecule mol;
    mol.charge       = tbl["system"]["charge"].value_or(0);
    mol.multiplicity = tbl["system"]["multiplicity"].value_or(1);

    auto atoms_arr = tbl["atoms"].as_array();
    if (!atoms_arr)
        throw std::runtime_error("TOML: no [[atoms]] array found.");

    for (auto& node : *atoms_arr) {
        auto* atbl = node.as_table();
        if (!atbl) throw std::runtime_error("TOML: each [[atoms]] entry must be a table.");
        Atom a;
        a.symbol    = (*atbl)["symbol"].value<std::string>().value();
        a.Z         = atomic_number(a.symbol);
        double x    = (*atbl)["x"].value<double>().value_or(0.0);
        double y    = (*atbl)["y"].value<double>().value_or(0.0);
        double z    = (*atbl)["z"].value<double>().value_or(0.0);
        a.coords    = { x * ANG2BOHR, y * ANG2BOHR, z * ANG2BOHR };
        mol.atoms.push_back(std::move(a));
    }

    // Sanity check
    int ne = mol.n_electrons();
    if (ne <= 0)
        throw std::runtime_error("Molecule has 0 or negative electrons.");
    int unpaired = mol.multiplicity - 1;
    if ((ne % 2) != (unpaired % 2))
        throw std::runtime_error("Inconsistent charge/multiplicity (parity mismatch).");
    if (mol.n_alpha() < 0 || mol.n_beta() < 0)
        throw std::runtime_error("Negative alpha/beta electron count — check charge/multiplicity.");

    return mol;
}

// ── Parse SCF settings ────────────────────────────────────────────────────
static SCFSettings parse_scf_settings(const toml::table& tbl) {
    SCFSettings s;
    if (auto* scf = tbl["scf"].as_table()) {
        s.max_iter    = (*scf)["max_iter"].value_or(100);
        s.convergence = (*scf)["convergence"].value_or(1e-8);
        s.diis_size   = (*scf)["diis_size"].value_or(8);
        s.screen_tol  = (*scf)["screen_tol"].value_or(1e-10);
    }
    return s;
}

// ── Main entry point ──────────────────────────────────────────────────────
int main(int argc, char* argv[]) {
    // Disable gradient tracking — we only need forward-mode tensor ops
    torch::NoGradGuard no_grad;

    std::string input_path = "input.toml";
    if (argc >= 2) input_path = argv[1];

    fs::path inp_path(input_path);
    if (!fs::exists(inp_path)) {
        std::cerr << "Error: input file not found: " << input_path << "\n";
        return 1;
    }

    fs::path input_dir = inp_path.parent_path();
    if (input_dir.empty()) input_dir = ".";

    fs::path exe_dir = fs::path(argv[0]).parent_path();
    if (exe_dir.empty()) exe_dir = ".";

    std::cout << "╔══════════════════════════════════════════════════════════╗\n";
    std::cout << "║           Hartree-Fock / C++23 + LibTorch                ║\n";
    std::cout << "╚══════════════════════════════════════════════════════════╝\n";
    std::cout << "  Input: " << fs::absolute(inp_path) << "\n";

    // Parse TOML
    toml::table tbl;
    try {
        tbl = toml::parse_file(inp_path.string());
    } catch (const toml::parse_error& e) {
        std::cerr << "TOML parse error: " << e.what() << "\n";
        return 1;
    }

    Molecule mol;
    SCFSettings cfg;
    try {
        mol = parse_molecule(tbl);
        cfg = parse_scf_settings(tbl);
    } catch (const std::exception& e) {
        std::cerr << "Input error: " << e.what() << "\n";
        return 1;
    }

    // Basis set
    std::string basis_name = tbl["system"]["basis"].value_or(std::string("sto-3g"));
    fs::path basis_file;
    try {
        basis_file = find_basis_file(basis_name, input_dir, exe_dir);
    } catch (const std::exception& e) {
        std::cerr << "Basis error: " << e.what() << "\n";
        return 1;
    }
    std::cout << "  Basis: " << basis_name << "  (" << basis_file << ")\n";

    BasisSet bs;
    try {
        bs = load_basis(mol, basis_file);
    } catch (const std::exception& e) {
        std::cerr << "Basis load error: " << e.what() << "\n";
        return 1;
    }

    print_molecule(mol);
    print_basis_info(bs, mol);

    // Decide RHF vs UHF
    bool use_uhf = (mol.multiplicity != 1) || (mol.n_electrons() % 2 != 0);

    try {
        if (use_uhf) {
            auto res = run_uhf(mol, bs, cfg);
            print_uhf_result(res, mol, bs);
        } else {
            auto res = run_rhf(mol, bs, cfg);
            print_rhf_result(res, mol, bs);
        }
    } catch (const std::exception& e) {
        std::cerr << "SCF error: " << e.what() << "\n";
        return 1;
    }

    std::cout << "\n  Done.\n";
    return 0;
}

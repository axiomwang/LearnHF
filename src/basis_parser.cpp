// BSE Gaussian format parser
// Format (from Basis Set Exchange v0.9):
//   Lines starting with '!' are comments
//   "ELEMENT  0"  starts an element block
//   "SHELL_TYPE  N_PRIM  SCALE" defines a contracted shell
//     SHELL_TYPE ∈ {S, P, D, F, G, H, SP}
//   Following lines have (exponent  coeff) or (exponent  coeff_s  coeff_p) for SP
//   "****" terminates an element block
//   "ECP" starts an ECP block; skip until next "****"
//
// SP shells are split into a separate S and P shell with the same primitives.
// Primitive normalization is applied here so SCF code uses pre-normalized coefficients.

#include "basis.hpp"
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <iostream>
#include <algorithm>
#include <cctype>
#include <cmath>
#include <map>

// ── Normalization of a Cartesian GTO primitive ──────────────────────────────
// N(α, nx, ny, nz) = [ (2α/π)^(3/4) * (4α)^L / sqrt(product of double_fact) ]
// where L = nx+ny+nz and double_fact = (2nx-1)!! * (2ny-1)!! * (2nz-1)!!
// For a shell we normalise assuming all exponents contribute equally
// (the "Pitzer" convention used by BSE).
// For the overall contraction normalization we follow Szabo & Ostlund.

static double double_factorial(int n) {
    // (2n-1)!! = 1*3*5*..*(2n-1);  (-1)!! = 1 by convention
    if (n <= 0) return 1.0;
    double r = 1.0;
    for (int i = n; i >= 1; i -= 2) r *= i;
    return r;
}

// Normalization of primitive GTO with exponent alpha and angular momentum l
// (using the spherically-averaged formula, suitable for all Cartesian components)
static double prim_norm(double alpha, int l) {
    // N = (2α/π)^(3/4) * sqrt( (4α)^l / (2l-1)!! )
    // where (2l-1)!! = df(2l-1)
    double df   = double_factorial(2*l - 1);
    double pow4 = std::pow(4.0 * alpha, l);
    double N    = std::pow(2.0 * alpha / M_PI, 0.75) * std::sqrt(pow4 / df);
    return N;
}

// ── String utilities ─────────────────────────────────────────────────────────
static std::string trim(const std::string& s) {
    size_t a = s.find_first_not_of(" \t\r\n");
    if (a == std::string::npos) return "";
    size_t b = s.find_last_not_of(" \t\r\n");
    return s.substr(a, b - a + 1);
}

static std::string to_upper(std::string s) {
    for (auto& c : s) c = (char)std::toupper((unsigned char)c);
    return s;
}

// Angular momentum token to integer
static int shell_type_to_l(const std::string& t) {
    if (t == "S")  return 0;
    if (t == "P")  return 1;
    if (t == "D")  return 2;
    if (t == "F")  return 3;
    if (t == "G")  return 4;
    if (t == "H")  return 5;
    if (t == "I")  return 6;
    throw std::runtime_error("Unknown shell type: " + t);
}

// ── Parser ───────────────────────────────────────────────────────────────────

// Data for one element extracted from the basis file
struct ElemBasisData {
    struct RawShell {
        char type;       // 'S','P','D','F','G','H','I', or '*' for SP
        std::vector<std::pair<double,double>> prims_s; // (exp, coeff_s)
        std::vector<std::pair<double,double>> prims_p; // (exp, coeff_p) [SP only]
    };
    std::vector<RawShell> shells;
};

static std::map<std::string, ElemBasisData>
parse_basis_file(const std::filesystem::path& path) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("Cannot open basis file: " + path.string());

    std::map<std::string, ElemBasisData> result;
    std::string line;
    std::string cur_elem;
    ElemBasisData::RawShell cur_shell;
    bool in_shell = false;
    bool in_ecp   = false;

    auto flush_shell = [&]() {
        if (in_shell && !cur_shell.prims_s.empty()) {
            result[cur_elem].shells.push_back(cur_shell);
        }
        cur_shell = {};
        in_shell  = false;
    };

    while (std::getline(in, line)) {
        std::string s = trim(line);
        if (s.empty() || s[0] == '!') continue;

        // ECP block start
        if (s.size() >= 3 && to_upper(s.substr(0,3)) == "ECP") {
            in_ecp = true;
            flush_shell();
            continue;
        }
        // "****" terminates element block (and ECP blocks)
        if (s == "****") {
            flush_shell();
            cur_elem.clear();
            in_ecp = false;
            continue;
        }
        if (in_ecp) continue;

        // Element header: "SYMBOL  0"
        // Heuristic: if there are exactly 2 tokens, the second is an integer,
        // and the first is a known/reasonable element symbol.
        {
            std::istringstream iss(s);
            std::string tok1, tok2;
            iss >> tok1 >> tok2;
            std::string rest;
            if (!(iss >> rest) && !tok2.empty()) {
                // possibly "ELEMENT  0" or "END" or shell definition
                // Shell definition lines have a letter token + integer + float
                // Element line has letter + "0"
                bool tok2_is_zero = (tok2 == "0");
                // Check the first token is a plausible element symbol
                // (1-2 chars, first uppercase letter)
                bool looks_like_elem = !tok1.empty()
                    && std::isupper((unsigned char)tok1[0])
                    && tok1.size() <= 3
                    && tok2_is_zero;
                if (looks_like_elem) {
                    flush_shell();
                    cur_elem = tok1;
                    in_ecp   = false;
                    continue;
                }
            }
        }

        if (cur_elem.empty()) continue;

        // Shell definition line: "TYPE  NPRIM  SCALE"
        // TYPE is one or more uppercase letters (S/P/D/F/G/H/SP)
        {
            std::istringstream iss(s);
            std::string type_tok;
            int nprim;
            double scale;
            if (iss >> type_tok >> nprim >> scale) {
                std::string ut = to_upper(type_tok);
                if (ut == "S" || ut == "P" || ut == "D" || ut == "F" ||
                    ut == "G" || ut == "H" || ut == "I" || ut == "SP") {
                    flush_shell();
                    in_shell       = true;
                    cur_shell.type = (ut == "SP") ? '*' : ut[0];
                    cur_shell.prims_s.reserve((size_t)nprim);
                    cur_shell.prims_p.reserve((ut == "SP") ? (size_t)nprim : 0);
                    continue;
                }
            }
        }

        // Primitive data line
        if (in_shell) {
            std::istringstream iss(s);
            // Replace D/d exponent notation with E/e
            std::string fixed_s;
            fixed_s.reserve(s.size());
            for (char c : s) {
                if ((c == 'D' || c == 'd') &&
                    !fixed_s.empty() &&
                    (std::isdigit((unsigned char)fixed_s.back()) || fixed_s.back() == '.'))
                    fixed_s += 'E';
                else
                    fixed_s += c;
            }
            std::istringstream iss2(fixed_s);
            double exp_val, cs;
            if (!(iss2 >> exp_val >> cs)) continue;
            if (cur_shell.type == '*') {
                double cp;
                if (!(iss2 >> cp)) cp = 0.0;
                cur_shell.prims_s.push_back({exp_val, cs});
                cur_shell.prims_p.push_back({exp_val, cp});
            } else {
                cur_shell.prims_s.push_back({exp_val, cs});
            }
        }
    }
    // Flush any remaining data
    flush_shell();
    return result;
}

// ── Public API ───────────────────────────────────────────────────────────────

BasisSet load_basis(const Molecule& mol, const std::filesystem::path& basis_file) {
    auto elem_data = parse_basis_file(basis_file);

    BasisSet bs;
    for (size_t ai = 0; ai < mol.atoms.size(); ++ai) {
        const Atom& atom = mol.atoms[ai];
        auto it = elem_data.find(atom.symbol);
        if (it == elem_data.end()) {
            throw std::runtime_error("Basis set \"" + basis_file.filename().string() +
                                     "\" does not contain element: " + atom.symbol);
        }
        for (auto& raw : it->second.shells) {
            if (raw.type == '*') {
                // SP shell → split into S and P
                for (int shell_l : {0, 1}) {
                    Shell sh;
                    sh.atom_idx = (int)ai;
                    sh.center   = atom.coords;
                    sh.l        = shell_l;
                    const auto& src = (shell_l == 0) ? raw.prims_s : raw.prims_p;
                    sh.prims.reserve(src.size());
                    for (auto& [exp_v, coeff] : src) {
                        double N = prim_norm(exp_v, shell_l);
                        sh.prims.push_back({exp_v, coeff * N});
                    }
                    bs.shells.push_back(std::move(sh));
                }
            } else {
                Shell sh;
                sh.atom_idx = (int)ai;
                sh.center   = atom.coords;
                sh.l        = shell_type_to_l(std::string(1, raw.type));
                sh.prims.reserve(raw.prims_s.size());
                for (auto& [exp_v, coeff] : raw.prims_s) {
                    double N = prim_norm(exp_v, sh.l);
                    sh.prims.push_back({exp_v, coeff * N});
                }
                bs.shells.push_back(std::move(sh));
            }
        }
    }
    return bs;
}

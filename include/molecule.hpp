#pragma once
#include <string>
#include <vector>
#include <array>
#include <stdexcept>
#include <unordered_map>
#include <cmath>

// Angstroms → Bohr conversion
inline constexpr double ANG2BOHR = 1.8897259886;

// Atomic number lookup (symbol → Z)
inline int atomic_number(const std::string& sym) {
    static const std::unordered_map<std::string, int> table = {
        {"H",1},{"He",2},{"Li",3},{"Be",4},{"B",5},{"C",6},{"N",7},{"O",8},
        {"F",9},{"Ne",10},{"Na",11},{"Mg",12},{"Al",13},{"Si",14},{"P",15},
        {"S",16},{"Cl",17},{"Ar",18},{"K",19},{"Ca",20},{"Sc",21},{"Ti",22},
        {"V",23},{"Cr",24},{"Mn",25},{"Fe",26},{"Co",27},{"Ni",28},{"Cu",29},
        {"Zn",30},{"Ga",31},{"Ge",32},{"As",33},{"Se",34},{"Br",35},{"Kr",36},
        {"Rb",37},{"Sr",38},{"Y",39},{"Zr",40},{"Nb",41},{"Mo",42},{"Tc",43},
        {"Ru",44},{"Rh",45},{"Pd",46},{"Ag",47},{"Cd",48},{"In",49},{"Sn",50},
        {"Sb",51},{"Te",52},{"I",53},{"Xe",54}
    };
    auto it = table.find(sym);
    if (it == table.end())
        throw std::runtime_error("Unknown element symbol: " + sym);
    return it->second;
}

inline std::string atomic_symbol(int Z) {
    static const std::vector<std::string> syms = {
        "","H","He","Li","Be","B","C","N","O","F","Ne",
        "Na","Mg","Al","Si","P","S","Cl","Ar","K","Ca",
        "Sc","Ti","V","Cr","Mn","Fe","Co","Ni","Cu","Zn",
        "Ga","Ge","As","Se","Br","Kr","Rb","Sr","Y","Zr",
        "Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd","In","Sn",
        "Sb","Te","I","Xe"
    };
    if (Z < 1 || Z >= (int)syms.size())
        throw std::runtime_error("Unknown atomic number: " + std::to_string(Z));
    return syms[Z];
}

struct Atom {
    std::string symbol;
    int Z;                        // atomic number
    std::array<double,3> coords;  // in Bohr
};

struct Molecule {
    std::vector<Atom> atoms;
    int charge       = 0;
    int multiplicity = 1;   // 2S+1; 1 → RHF, >1 → UHF

    int n_electrons() const {
        int total_Z = 0;
        for (auto& a : atoms) total_Z += a.Z;
        return total_Z - charge;
    }

    int n_alpha() const {
        int ne = n_electrons();
        int unpaired = multiplicity - 1;
        return (ne + unpaired) / 2;
    }

    int n_beta() const {
        return n_electrons() - n_alpha();
    }

    // Nuclear repulsion energy (Hartree)
    double nuclear_repulsion() const {
        double E = 0.0;
        for (size_t i = 0; i < atoms.size(); ++i)
            for (size_t j = i+1; j < atoms.size(); ++j) {
                double r2 = 0.0;
                for (int k = 0; k < 3; ++k) {
                    double d = atoms[i].coords[k] - atoms[j].coords[k];
                    r2 += d*d;
                }
                E += atoms[i].Z * atoms[j].Z / std::sqrt(r2);
            }
        return E;
    }
};

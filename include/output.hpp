#pragma once
#include "scf.hpp"
#include "molecule.hpp"

void print_molecule(const Molecule& mol);
void print_basis_info(const BasisSet& bs, const Molecule& mol);
void print_rhf_result(const RHFResult& res, const Molecule& mol, const BasisSet& bs);
void print_uhf_result(const UHFResult& res, const Molecule& mol, const BasisSet& bs);

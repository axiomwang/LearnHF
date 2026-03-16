``output.hpp`` / ``output.cpp`` — Result Printing
==================================================

.. contents:: On this page
   :depth: 2
   :local:

Overview
--------

Pretty-printing functions that format and display molecular geometry,
basis set information, and SCF results to ``stdout`` using box-drawing
Unicode characters.

**Header:** ``include/output.hpp``
**Source:** ``src/output.cpp``

Functions
---------

.. cpp:function:: void print_molecule(const Molecule& mol)

   Prints a formatted geometry table showing each atom's symbol,
   atomic number, and Cartesian coordinates in both Bohr and Ångströms.

   Example output:

   .. code-block:: none

      ╔═══════════════════════════════════════════════╗
      ║               Molecular Geometry              ║
      ╠═══════════════════════════════════════════════╣
      ║  Atom    Z        X (Bohr)    Y (Bohr)    Z (Bohr) ║
      ╠═══════════════════════════════════════════════╣
      ║   O      8      0.000000    0.000000    0.222517 ║
      ║   H      1      0.000000    1.430901   -0.890067 ║
      ║   H      1      0.000000   -1.430901   -0.890067 ║
      ╠═══════════════════════════════════════════════╣
      ║  Charge: 0   Multiplicity: 1   Electrons: 10  ║
      ╚═══════════════════════════════════════════════╝

.. cpp:function:: void print_basis_info(const BasisSet& bs, const Molecule& mol)

   Prints a summary of the loaded basis set: number of shells, total
   number of basis functions, and a per-element shell breakdown.

.. cpp:function:: void print_rhf_result(const RHFResult& res, const Molecule& mol, const BasisSet& bs)

   :param res: Converged (or attempted) RHF result.
   :param mol: Molecule (needed for :math:`E_\text{nuc}`).
   :param bs: Basis set (needed to determine occupied/virtual boundary).

   Prints:

   * Convergence status and number of iterations.
   * Electronic, nuclear-repulsion, and total energies.
   * Orbital energy table (up to first 20 MOs) with occupation.
   * HOMO energy, LUMO energy, and HOMO-LUMO gap.

   Example energy block:

   .. code-block:: none

      ╔══════════════════════════════════╗
      ║          RHF Results             ║
      ╠══════════════════════════════════╣
      ║  Converged:   YES  Iterations: 12 ║
      ║  E_elec :   -84.1978 Ha          ║
      ║  E_nuc  :     9.1806 Ha          ║
      ║  E_total:   -75.0172 Ha          ║
      ╠══════════════════════════════════╣
      ║  HOMO :  -0.3971 Ha              ║
      ║  LUMO :   0.5758 Ha              ║
      ║  Gap  :   0.9729 Ha              ║
      ╚══════════════════════════════════╝

.. cpp:function:: void print_uhf_result(const UHFResult& res, const Molecule& mol, const BasisSet& bs)

   :param res: Converged (or attempted) UHF result.
   :param mol: Molecule.
   :param bs: Basis set.

   All output of :cpp:func:`print_rhf_result` plus:

   * Separate alpha and beta orbital energy tables.
   * :math:`\langle\hat{S}^2\rangle_\text{exact}` and
     :math:`\langle\hat{S}^2\rangle_\text{UHF}`.
   * A spin-contamination warning if

     .. math::

        \langle\hat{S}^2\rangle_\text{UHF} - \langle\hat{S}^2\rangle_\text{exact}
        > 0.1 \cdot \langle\hat{S}^2\rangle_\text{exact}.

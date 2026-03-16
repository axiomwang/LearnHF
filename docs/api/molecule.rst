``molecule.hpp`` — Molecular Data Structures
=============================================

.. contents:: On this page
   :depth: 2
   :local:

Overview
--------

Defines the core chemistry data types: the atomic coordinates, the
molecular system (charge, multiplicity, atoms), and several helper
functions for element lookup and physical constants.

**Source file:** ``include/molecule.hpp``

Constants
---------

.. c:macro:: ANG2BOHR

   Conversion factor from Ångströms to Bohr radii.

   .. code-block:: cpp

      inline constexpr double ANG2BOHR = 1.8897259886;

   All atomic coordinates stored internally are in **Bohr**.  Coordinates
   read from the TOML input file (assumed Ångströms) are multiplied by
   this constant immediately on parse.

Utility Functions
-----------------

.. cpp:function:: int atomic_number(const std::string& sym)

   :param sym: Chemical symbol string (case-sensitive), e.g. ``"O"``, ``"Fe"``.
   :returns: Atomic number :math:`Z`.
   :throws std::runtime_error: If the symbol is not found in the look-up table.

   Returns the atomic number for a given element symbol.
   Covers H (1) through Xe (54).

   .. code-block:: cpp

      int z = atomic_number("C");  // z == 6

.. cpp:function:: std::string atomic_symbol(int Z)

   :param Z: Atomic number 1–54.
   :returns: Element symbol string.
   :throws std::runtime_error: If *Z* is out of the supported range.

   Inverse of :cpp:func:`atomic_number`.

   .. code-block:: cpp

      std::string sym = atomic_symbol(8);  // sym == "O"

Data Structures
---------------

.. cpp:struct:: Atom

   Represents a single nucleus.

   .. cpp:member:: std::string symbol

      Element symbol, e.g. ``"H"``.

   .. cpp:member:: int Z

      Atomic number.

   .. cpp:member:: std::array<double,3> coords

      Cartesian coordinates **(in Bohr)**.

.. cpp:struct:: Molecule

   A complete molecular system.

   .. cpp:member:: std::vector<Atom> atoms

      Ordered list of atoms.

   .. cpp:member:: int charge

      Net molecular charge (default ``0``).

   .. cpp:member:: int multiplicity

      Spin multiplicity :math:`2S+1` (default ``1``).
      ``1`` → singlet (RHF eligible); ``2`` → doublet (UHF), etc.

   .. cpp:function:: int n_electrons() const

      Returns :math:`\sum_A Z_A - \text{charge}`.

   .. cpp:function:: int n_alpha() const

      Number of :math:`\alpha`-spin electrons:
      :math:`(N_e + (M-1))/2` where :math:`M = \text{multiplicity}`.

   .. cpp:function:: int n_beta() const

      Number of :math:`\beta`-spin electrons:
      :math:`N_e - N_\alpha`.

   .. cpp:function:: double nuclear_repulsion() const

      Nuclear repulsion energy in Hartree:

      .. math::

         E_\text{nuc} = \sum_{A < B} \frac{Z_A Z_B}{R_{AB}}.

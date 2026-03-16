``basis.hpp`` / ``basis_parser.cpp`` — Basis Set Types & Loader
===============================================================

.. contents:: On this page
   :depth: 2
   :local:

Overview
--------

Defines the contracted-Gaussian basis-set data structures and a parser
for Basis Set Exchange (BSE) Gaussian-format files.

**Header:** ``include/basis.hpp``
**Source:** ``src/basis_parser.cpp``

Data Structures
---------------

.. cpp:struct:: Primitive

   One uncontracted Gaussian primitive contribution.

   .. cpp:member:: double exp

      Exponent :math:`\alpha_k`.

   .. cpp:member:: double coeff

      Contraction coefficient :math:`d_k`, **already multiplied** by the
      primitive normalisation constant

      .. math::

         N_\text{prim}(\alpha, \ell)
         = \left(\frac{2\alpha}{\pi}\right)^{3/4}
           \sqrt{\frac{(4\alpha)^\ell}{(2\ell-1)!!}}

      during basis loading.

.. cpp:struct:: Shell

   A contracted Gaussian shell: all primitives share the same angular
   momentum :math:`\ell` and nuclear centre :math:`\mathbf{A}`.

   The shell function is

   .. math::

      \chi_\mu(\mathbf{r}) = N_\mu\,
          x^{n_x} y^{n_y} z^{n_z}\,
          \sum_k d_k\, e^{-\alpha_k |\mathbf{r}-\mathbf{A}|^2}

   where :math:`(n_x,n_y,n_z)` is one of the :math:`N_\text{cart}(\ell)`
   Cartesian components and :math:`N_\mu` is a shell-level normalisation
   factor applied after contraction (see :ref:`algorithms:Shell Normalization`).

   .. cpp:member:: int atom_idx

      Index into the parent :cpp:struct:`Molecule`'s atoms array.

   .. cpp:member:: std::array<double,3> center

      Shell centre in **Bohr** (copied from the atom's coordinates).

   .. cpp:member:: int l

      Angular momentum quantum number (0 = S, 1 = P, 2 = D, 3 = F, 4 = G).

   .. cpp:member:: std::vector<Primitive> prims

      Primitives with pre-normalised coefficients.

   .. cpp:function:: int n_cart() const

      Returns :math:`(\ell+1)(\ell+2)/2` — the number of distinct
      Cartesian functions in this shell.

.. cpp:struct:: BasisSet

   The complete AO basis for one molecule.

   .. cpp:member:: std::vector<Shell> shells

      Ordered shell list — the ordering determines the AO basis-function
      index mapping.

   .. cpp:function:: int n_basis() const

      Total number of Cartesian basis functions
      :math:`N_b = \sum_i N_\text{cart}(\ell_i)`.

   .. cpp:function:: int shell_start(int i) const

      Returns the index of the first basis function in shell *i*
      (i.e. :math:`\sum_{j<i} N_\text{cart}(\ell_j)`).

Cartesian Ordering
------------------

.. cpp:function:: std::vector<std::array<int,3>> cart_components(int l)

   Returns the ordered list of :math:`(n_x, n_y, n_z)` tuples for all
   Cartesian functions with total angular momentum :math:`\ell = l`.

   **Ordering convention** — iterates ``nx`` descending from ``l`` to ``0``;
   for each ``nx`` iterates ``ny`` descending from ``l-nx`` to ``0``;
   ``nz = l - nx - ny``.  This matches the Gaussian09/Psi4 convention.

   Example for :math:`\ell = 2` (D shell):

   .. code-block:: none

      (2,0,0) d_xx
      (1,1,0) d_xy
      (1,0,1) d_xz
      (0,2,0) d_yy
      (0,1,1) d_yz
      (0,0,2) d_zz

Basis Loader
------------

.. cpp:function:: BasisSet load_basis(const Molecule& mol, const std::filesystem::path& basis_file)

   :param mol: The molecule; element symbols are used to select the
               correct sections from the file.
   :param basis_file: Filesystem path to a BSE Gaussian-format basis file.
   :returns: A fully normalised :cpp:struct:`BasisSet`.
   :throws std::runtime_error: On file-not-found or parse errors.

   **File format** — BSE Gaussian-format (``#FORMAT`` tag optional).
   Each element block begins with ``****``, followed by the element symbol
   and blocks of shell definitions such as:

   .. code-block:: none

      O   0
      S   3   1.00
           130.7093200   0.1543290
            23.8088610   0.5353281
             5.4424170   0.4446345
      P   3   1.00
             5.0331513  -0.9996723
             1.1695961   1.5591627
             0.3803890   0.9444749
      ****

   Special handling:

   * **SP shells** — lines with two coefficients (S and P on the same
     primitive) are split into two separate shells: one ``S`` and one ``P``.
   * **D-notation exponents** — ``1.23D+04`` is converted to ``1.23E+04``.
   * **ECP blocks** — detected by the ``ECP`` header keyword and skipped
     with a warning printed to ``stderr``.  Effective-core-potential
     integrals are *not* evaluated.

   After parsing, :ref:`shell normalisation <algorithms:Shell Normalization>`
   is applied to every shell.

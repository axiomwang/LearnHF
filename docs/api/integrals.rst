``integrals.hpp`` / ``integrals.cpp`` — Integral Engine
========================================================

.. contents:: On this page
   :depth: 2
   :local:

Overview
--------

Computes all one- and two-electron AO integrals using the
McMurchie-Davidson scheme.  All returned tensors are
``torch::kDouble``, two-dimensional (or rank-4 is never returned
directly — the Coulomb/exchange matrices are already contracted with a
density matrix).

**Header:** ``include/integrals.hpp``
**Source:** ``src/integrals.cpp``

For the mathematical background see :ref:`algorithms:McMurchie-Davidson Integral Engine`.

One-Electron Integrals
-----------------------

.. cpp:function:: torch::Tensor build_overlap(const BasisSet& bs)

   :returns: Symmetric :math:`N_b \times N_b` overlap matrix
             :math:`S_{\mu\nu} = \langle\mu|\nu\rangle`.

   Key property: :math:`S = X^{-T} X^{-1}` where :math:`X` is the
   canonical orthogonalization matrix (see :ref:`principles:Canonical Orthogonalization`).

.. cpp:function:: torch::Tensor build_kinetic(const BasisSet& bs)

   :returns: Symmetric :math:`N_b \times N_b` kinetic energy matrix
             :math:`T_{\mu\nu} = \langle\mu|-\tfrac{1}{2}\nabla^2|\nu\rangle`.

.. cpp:function:: torch::Tensor build_nuclear(const BasisSet& bs, const Molecule& mol)

   :param bs: Basis set.
   :param mol: Molecule; atomic charges :math:`Z_A` and positions
               :math:`\mathbf{R}_A` enter the nuclear-attraction sum.
   :returns: Symmetric :math:`N_b \times N_b` nuclear-attraction matrix

             .. math::

                V_{\mu\nu} = \sum_A \left\langle\mu\left|
                \frac{-Z_A}{|\mathbf{r}-\mathbf{R}_A|}\right|\nu\right\rangle.

.. cpp:function:: torch::Tensor build_hcore(const BasisSet& bs, const Molecule& mol)

   Convenience inline that returns
   :cpp:func:`build_kinetic` ``+`` :cpp:func:`build_nuclear`:

   .. math::

      H^{\rm core}_{\mu\nu} = T_{\mu\nu} + V_{\mu\nu}.

   Defined inline in the header.

Two-Electron Integrals
-----------------------

.. cpp:function:: torch::Tensor build_schwarz(const BasisSet& bs)

   :returns: Symmetric :math:`N_b \times N_b` Schwarz bound matrix
             :math:`Q_{\mu\nu} = \sqrt{(\mu\nu|\mu\nu)}`.

   Used to screen shell quartets and individual basis-function pairs in
   :cpp:func:`build_JK`.  Call this **once** before the SCF iteration and
   pass the result to every :cpp:func:`build_JK` call.

.. cpp:function:: std::pair<torch::Tensor, torch::Tensor> build_JK(const BasisSet& bs, const torch::Tensor& P, const torch::Tensor& screen, double tol = 1e-10)

   :param bs: Basis set.
   :param P: AO density matrix :math:`\mathbf{P}`, shape
             :math:`[N_b, N_b]`, ``torch::kDouble``.
   :param screen: Schwarz matrix from :cpp:func:`build_schwarz`.
   :param tol: Schwarz screening threshold (default ``1e-10``).
   :returns: ``{J, K}`` — Coulomb and exchange matrices,
             both :math:`N_b \times N_b`.

   Computes in a single loop over shell quartets:

   .. math::

      J_{\mu\nu} &= \sum_{\lambda\sigma} P_{\lambda\sigma}\,(\mu\nu|\lambda\sigma) \\
      K_{\mu\nu} &= \sum_{\lambda\sigma} P_{\lambda\sigma}\,(\mu\sigma|\lambda\nu)

   exploiting **8-fold permutational symmetry** and **Schwarz screening**.

   Two screening levels are applied:

   1. **Shell level**: skip entire shell quartet
      :math:`(\mathbf{a}\mathbf{b}|\mathbf{c}\mathbf{d})` if
      :math:`Q_{ab}^\text{shell} \cdot Q_{cd}^\text{shell} < \varepsilon`.
   2. **Function level**: within a surviving quartet, skip individual
      :math:`(\mu\nu|\lambda\sigma)` if
      :math:`|Q_{\mu\nu}| \cdot |Q_{\lambda\sigma}| < \varepsilon`.

   The Fock matrix is then assembled by the caller as
   :math:`F = H^{\rm core} + J - \tfrac{1}{2}K` (RHF) or the
   spin-specific form for UHF.

Internal Functions (implementation details)
--------------------------------------------

The following are in the anonymous namespace of ``integrals.cpp`` and
not part of the public API, but are documented here for reference.

``fill_E(E, alpha, beta, XAB)``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Fills the :math:`9 \times 9 \times 18` E-coefficient table for one
Cartesian direction with exponents :math:`\alpha`, :math:`\beta` and
centre separation :math:`X_{AB} = A - B` (in Bohr).

Uses horizontal recursion for :math:`j = 0` first, then vertical
recursion for :math:`j > 0`.

``fill_R(R, p, PC, nmax)``
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Fills the :math:`20 \times 20 \times 20 \times 20` R-integral table
:math:`R^n_{tuv}` up to order *nmax*, for exponent *p* and
displacement :math:`\mathbf{P}-\mathbf{C}`.

The Boys values at each base index are supplied externally via the
global Boys vector computed once per shell pair/quartet.

``prim_norm(alpha, l)``
~~~~~~~~~~~~~~~~~~~~~~~

Returns the normalisation constant for a primitive Gaussian with
exponent :math:`\alpha` and angular momentum :math:`\ell`:

.. math::

   N = \left(\frac{2\alpha}{\pi}\right)^{3/4}
       \sqrt{\frac{(4\alpha)^\ell}{(2\ell-1)!!}}

Used during basis loading in ``basis_parser.cpp``.

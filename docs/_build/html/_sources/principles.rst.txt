Theory & Principles
====================

.. contents:: On this page
   :depth: 3
   :local:

This page summarises the quantum-chemical theory underlying **hftest**.
For full derivations the reader is referred to the canonical texts
[Szabo1989]_, [HJO2000]_, and [Cramer2004]_.

The Electronic Schrödinger Equation
-------------------------------------

Within the Born-Oppenheimer approximation the nuclei are treated as
fixed classical point charges.  The goal is to solve

.. math::
   :label: schr

   \hat{H}_\text{el}\, \Psi = E_\text{el}\, \Psi

for the *N*-electron wave function :math:`\Psi`, where

.. math::

   \hat{H}_\text{el} = -\frac{1}{2}\sum_{i=1}^N \nabla_i^2
       - \sum_{i=1}^N \sum_{A=1}^M \frac{Z_A}{r_{iA}}
       + \sum_{i=1}^N \sum_{j>i}^N \frac{1}{r_{ij}}.

The total energy (observable) is

.. math::

   E_\text{total} = E_\text{el} + E_\text{nuc}\,,\qquad
   E_\text{nuc} = \sum_{A < B} \frac{Z_A Z_B}{R_{AB}}\,.

The Hartree-Fock Approximation
--------------------------------

HF approximates :math:`\Psi` by a single *Slater determinant*:

.. math::

   \Psi_\text{HF} = \frac{1}{\sqrt{N!}}\det\!\bigl[\chi_1\chi_2\cdots\chi_N\bigr]

where :math:`\chi_i` are *spin orbitals* — products of a spatial
orbital :math:`\phi_i(\mathbf{r})` and a spin function
(:math:`\alpha` or :math:`\beta`).

The variational principle applied to this determinant yields the
**Hartree-Fock equations**

.. math::
   :label: hfeq

   \hat{F}\,\phi_i = \varepsilon_i\,\phi_i\,,

where :math:`\hat{F}` is the **Fock operator**

.. math::

   \hat{F}(\mathbf{r}_1) = \hat{h}(\mathbf{r}_1)
       + \sum_{j=1}^N \Bigl[\hat{J}_j(\mathbf{r}_1) - \hat{K}_j(\mathbf{r}_1)\Bigr].

The one-electron Hamiltonian :math:`\hat{h}` contains kinetic energy
and nuclear attraction; :math:`\hat{J}_j` is the **Coulomb operator**
(mean-field repulsion) and :math:`\hat{K}_j` is the non-local
**exchange operator** (a pure quantum-mechanical effect with no
classical analogue).

LCAO Expansion and the Roothaan Equations
-------------------------------------------

Molecular orbitals are expanded in a finite *atomic-orbital* (AO)
basis :math:`\{\mu\}` of size :math:`N_b`:

.. math::

   \phi_i(\mathbf{r}) = \sum_{\mu=1}^{N_b} C_{\mu i}\,\chi_\mu(\mathbf{r}).

Substituting into :eq:`hfeq` and projecting onto the basis gives the
**Roothaan-Hall** (closed-shell) or **Pople-Nesbet** (open-shell)
matrix equations:

.. math::
   :label: roothaan

   \mathbf{F}\,\mathbf{C} = \mathbf{S}\,\mathbf{C}\,\boldsymbol{\varepsilon}

where :math:`\mathbf{S}_{\mu\nu} = \langle\mu|\nu\rangle` is the AO
overlap matrix, :math:`\mathbf{C}` is the MO coefficient matrix, and
:math:`\boldsymbol{\varepsilon}` is diagonal with orbital energies.

The Fock matrix in the AO basis is

.. math::

   F_{\mu\nu} = H^{\rm core}_{\mu\nu}
       + \sum_{\lambda\sigma} P_{\lambda\sigma}
         \Bigl[(\mu\nu|\lambda\sigma) - \tfrac{1}{2}(\mu\sigma|\lambda\nu)\Bigr]

where :math:`(\mu\nu|\lambda\sigma) = \iint
\chi_\mu(\mathbf{r}_1)\chi_\nu(\mathbf{r}_1)\,\frac{1}{r_{12}}\,
\chi_\lambda(\mathbf{r}_2)\chi_\sigma(\mathbf{r}_2)\,
d\mathbf{r}_1\,d\mathbf{r}_2`
are the **two-electron repulsion integrals** (ERIs), and the **density
matrix** is

.. math::

   P_{\mu\nu} = 2\sum_{i=1}^{N/2} C_{\mu i} C_{\nu i}
   \quad \text{(closed-shell, RHF)}.

Electronic Energy
~~~~~~~~~~~~~~~~~~

For RHF the converged electronic energy is

.. math::

   E_\text{el} = \tfrac{1}{2}\operatorname{Tr}\!\bigl[
       \mathbf{P}\bigl(\mathbf{H}^{\rm core} + \mathbf{F}\bigr)\bigr].

Restricted vs Unrestricted HF
-------------------------------

Restricted Hartree-Fock (RHF)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

All electrons occupy *paired* spatial orbitals:
:math:`\alpha`- and :math:`\beta`-spin electrons of the same pair share
the same spatial function.  RHF is only applicable to closed-shell
(singlet) systems.  It leads to the single Fock matrix above.

Unrestricted Hartree-Fock (UHF)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

:math:`\alpha`- and :math:`\beta`-spin electrons use **independent**
sets of spatial orbitals
:math:`\{\phi_i^\alpha\}` and :math:`\{\phi_j^\beta\}`.

The separate density matrices are

.. math::

   P^\alpha_{\mu\nu} = \sum_{i=1}^{N_\alpha} C^\alpha_{\mu i} C^\alpha_{\nu i}\,,\qquad
   P^\beta_{\mu\nu}  = \sum_{i=1}^{N_\beta}  C^\beta_{\mu i}  C^\beta_{\nu i}\,,

and the two Fock matrices are

.. math::

   F^\alpha_{\mu\nu} &= H^{\rm core}_{\mu\nu}
       + \sum_{\lambda\sigma}
         \bigl(P^\alpha_{\lambda\sigma}+P^\beta_{\lambda\sigma}\bigr)(\mu\nu|\lambda\sigma)
         - P^\alpha_{\lambda\sigma}(\mu\sigma|\lambda\nu) \\[4pt]
   F^\beta_{\mu\nu}  &= H^{\rm core}_{\mu\nu}
       + \sum_{\lambda\sigma}
         \bigl(P^\alpha_{\lambda\sigma}+P^\beta_{\lambda\sigma}\bigr)(\mu\nu|\lambda\sigma)
         - P^\beta_{\lambda\sigma}(\mu\sigma|\lambda\nu).

The electronic energy is

.. math::

   E_\text{el} = \tfrac{1}{2}\operatorname{Tr}\!\bigl[
       (\mathbf{P}^\alpha+\mathbf{P}^\beta)\mathbf{H}^{\rm core}
       + \mathbf{P}^\alpha\mathbf{F}^\alpha
       + \mathbf{P}^\beta \mathbf{F}^\beta\bigr].

Spin Contamination
~~~~~~~~~~~~~~~~~~~

UHF wave functions are not eigenfunctions of :math:`\hat{S}^2` in
general.  The expectation value

.. math::

   \langle\hat{S}^2\rangle_\text{UHF}
   = S(S+1) + N_\beta - \sum_{i,j} \bigl|\langle\phi^\alpha_i|\phi^\beta_j\rangle\bigl|^2

deviates from the exact value :math:`S(S+1)` by the *spin-contamination*
term, evaluated from the overlap of occupied alpha and beta MO blocks.

Gaussian Basis Sets
--------------------

The AO basis functions are *contracted Gaussian-type orbitals* (cGTOs):

.. math::

   \chi_\mu(\mathbf{r}) \propto
       \underbrace{x^{n_x} y^{n_y} z^{n_z}}_{\text{angular part}}\,
       \sum_k d_k\, e^{-\alpha_k |\mathbf{r}-\mathbf{A}|^2}

where :math:`(n_x,n_y,n_z)` specifies the Cartesian angular momentum,
:math:`\mathbf{A}` is the nuclear centre, the :math:`\alpha_k` are
primitive exponents, and the :math:`d_k` are contraction coefficients
pre-normalised during basis loading.

The angular quantum number :math:`\ell = n_x+n_y+n_z` labels the shell
type: :math:`\ell=0` (S), :math:`1` (P), :math:`2` (D), :math:`3` (F),
:math:`4` (G).  A shell of angular momentum :math:`\ell` contains

.. math::

   N_\text{cart}(\ell) = \tfrac{(\ell+1)(\ell+2)}{2}

Cartesian basis functions.

Canonical Orthogonalization
-----------------------------

Equation :eq:`roothaan` is a *generalised* eigenvalue problem because
:math:`\mathbf{S}` is not the identity when basis functions overlap.
**Canonical orthogonalization** transforms to an orthonormal basis:

1. Diagonalise  :math:`\mathbf{S} = \mathbf{U}\boldsymbol{\Lambda}\mathbf{U}^\top`.
2. Define  :math:`\mathbf{X} = \mathbf{U}\boldsymbol{\Lambda}^{-1/2}`, keeping
   only columns with :math:`\lambda_i \ge 10^{-8}` (discarding linearly
   dependent combinations).
3. Transform  :math:`\mathbf{F}' = \mathbf{X}^\top\mathbf{F}\mathbf{X}`.
4. Solve the *ordinary* eigenvalue problem  :math:`\mathbf{F}'\mathbf{C}' = \mathbf{C}'\boldsymbol{\varepsilon}`.
5. Back-transform  :math:`\mathbf{C} = \mathbf{X}\mathbf{C}'`.

References
----------

.. [Szabo1989]  A. Szabo and N. S. Ostlund, *Modern Quantum Chemistry*,
   Dover (1989).
.. [HJO2000]    T. Helgaker, P. Jørgensen, and J. Olsen,
   *Molecular Electronic-Structure Theory*, Wiley (2000).
.. [Cramer2004] C. J. Cramer, *Essentials of Computational Chemistry*,
   2nd ed., Wiley (2004).

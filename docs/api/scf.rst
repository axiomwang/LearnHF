``scf.hpp`` / ``scf.cpp`` — SCF Solvers
========================================

.. contents:: On this page
   :depth: 2
   :local:

Overview
--------

Self-consistent field (SCF) solvers for both closed-shell (RHF) and
open-shell (UHF) wave functions, with Pulay DIIS convergence
acceleration and canonical orthogonalization.

**Header:** ``include/scf.hpp``
**Source:** ``src/scf.cpp``

For the theoretical background see
:ref:`principles:Restricted vs Unrestricted HF` and
:ref:`algorithms:Pulay DIIS Convergence Acceleration`.

Settings
--------

.. cpp:struct:: SCFSettings

   All SCF parameters, parsed from the ``[scf]`` section of the TOML
   input.

   .. cpp:member:: int max_iter

      Maximum number of SCF iterations before stopping.
      Default ``100``.

   .. cpp:member:: double convergence

      Convergence threshold :math:`\varepsilon` applied to
      **both** :math:`|\Delta E|` and :math:`\text{RMS}(\Delta\mathbf{P})`.
      Default ``1e-8``.

   .. cpp:member:: int diis_size

      Number of (Fock, error) vector pairs retained in DIIS memory.
      When the queue is full the oldest pair is discarded.
      Default ``8``.

   .. cpp:member:: double screen_tol

      Schwarz screening threshold forwarded to
      :cpp:func:`build_JK`.  Default ``1e-10``.

Result Types
------------

.. cpp:struct:: RHFResult

   Converged RHF solution.

   .. cpp:member:: bool converged

      ``true`` if convergence was reached within ``max_iter``.

   .. cpp:member:: int n_iter

      Number of SCF iterations performed.

   .. cpp:member:: double E_elec

      Electronic energy :math:`E_\text{el}` in Hartree.

   .. cpp:member:: double E_nuc

      Nuclear repulsion energy in Hartree.

   .. cpp:member:: double E_total

      Total energy :math:`E_\text{total} = E_\text{el} + E_\text{nuc}`.

   .. cpp:member:: torch::Tensor C

      MO coefficient matrix, shape :math:`[N_b, N_b]`.
      Columns are MOs ordered by energy (lowest first).

   .. cpp:member:: torch::Tensor eps

      Orbital energies, shape :math:`[N_b]`, ascending order.

   .. cpp:member:: torch::Tensor P

      Final AO density matrix, shape :math:`[N_b, N_b]`.

   .. cpp:member:: torch::Tensor F

      Final Fock matrix (possibly DIIS-extrapolated), shape :math:`[N_b, N_b]`.

.. cpp:struct:: UHFResult

   Converged UHF solution.

   Follows the same energy fields as :cpp:struct:`RHFResult`, plus:

   .. cpp:member:: torch::Tensor Ca

      Alpha MO coefficient matrix.

   .. cpp:member:: torch::Tensor Cb

      Beta MO coefficient matrix.

   .. cpp:member:: torch::Tensor eps_a

      Alpha orbital energies (ascending).

   .. cpp:member:: torch::Tensor eps_b

      Beta orbital energies (ascending).

   .. cpp:member:: torch::Tensor Pa

      Alpha density matrix.

   .. cpp:member:: torch::Tensor Pb

      Beta density matrix.

   .. cpp:member:: double S2_exact

      Exact expectation value :math:`S(S+1)` where :math:`S = (N_\alpha - N_\beta)/2`.

   .. cpp:member:: double S2_uhf

      UHF expectation value of :math:`\hat{S}^2` including spin contamination:

      .. math::

         \langle\hat{S}^2\rangle_\text{UHF}
         = S(S+1) + N_\beta
           - \sum_{i,j}\bigl|\langle\phi^\alpha_i|\phi^\beta_j\rangle\bigr|^2.

Solver Functions
----------------

.. cpp:function:: RHFResult run_rhf(const Molecule& mol, const BasisSet& bs, const SCFSettings& settings)

   :param mol: Molecular system.
   :param bs: Fully loaded and normalised basis set.
   :param settings: SCF tuning parameters.
   :returns: :cpp:struct:`RHFResult` with converged (or best-so-far) data.

   **Algorithm:**

   1. Build :math:`\mathbf{S}`,  :math:`\mathbf{H}^{\rm core}`,
      Schwarz matrix :math:`\mathbf{Q}`.
   2. Canonical orthogonalization → :math:`\mathbf{X}` (threshold :math:`10^{-8}`).
   3. Core-Hamiltonian initial guess: diagonalise :math:`\mathbf{H}^{\rm core}`,
      form :math:`\mathbf{P}^{(0)}`.
   4. DIIS-accelerated iteration:

      a. Form :math:`\mathbf{J},\mathbf{K}` via :cpp:func:`build_JK` with density
         :math:`\mathbf{P}`.
      b. Assemble :math:`\mathbf{F} = \mathbf{H}^{\rm core} + \mathbf{J} - \tfrac{1}{2}\mathbf{K}`.
      c. Compute DIIS error :math:`\mathbf{e} = \mathbf{F}\mathbf{P}\mathbf{S} - \mathbf{S}\mathbf{P}\mathbf{F}`.
      d. Push :math:`(\mathbf{F}, \mathbf{e})` to DIIS; extrapolate :math:`\mathbf{F}_\text{DIIS}`.
      e. Compute energy :math:`E = \tfrac{1}{2}\operatorname{Tr}[\mathbf{P}(\mathbf{H}^{\rm core}+\mathbf{F})]`.
      f. Diagonalise :math:`\mathbf{X}^\top\mathbf{F}_\text{DIIS}\mathbf{X}`; form new :math:`\mathbf{P}`.
      g. Check :math:`|\Delta E| < \varepsilon` and :math:`\text{RMS}(\Delta\mathbf{P}) < \varepsilon`.

   Prints a convergence table to ``stdout`` during iteration.

.. cpp:function:: UHFResult run_uhf(const Molecule& mol, const BasisSet& bs, const SCFSettings& settings)

   :param mol: Molecular system.
   :param bs: Basis set.
   :param settings: SCF settings.
   :returns: :cpp:struct:`UHFResult`.

   Same overall algorithm as :cpp:func:`run_rhf` but maintains **separate**
   DIIS queues, density matrices, and Fock matrices for :math:`\alpha` and
   :math:`\beta` spins.

   The UHF energy functional is:

   .. math::

      E_\text{el} = \tfrac{1}{2}\operatorname{Tr}\bigl[
          (\mathbf{P}^\alpha+\mathbf{P}^\beta)\mathbf{H}^{\rm core}
          + \mathbf{P}^\alpha\mathbf{F}^\alpha
          + \mathbf{P}^\beta\mathbf{F}^\beta\bigr].

Internal: ``canonical_orthogonalization``
-----------------------------------------

.. cpp:function:: static std::pair<torch::Tensor, int> canonical_orthogonalization(const torch::Tensor& S, double thresh = 1e-8)

   Returns the transformation matrix :math:`\mathbf{X}` and the number of
   linearly independent columns retained.  Eigenvalues of :math:`\mathbf{S}`
   below *thresh* are discarded.

   Uses ``at::linalg_eigh(S, "L")`` which returns eigenvalues in ascending
   order; the last *n_ind* are kept.

Internal: ``DIIS``
-------------------

.. cpp:struct:: DIIS

   Pulay DIIS extrapolation state machine.

   .. cpp:function:: explicit DIIS(int sz)

      Construct with maximum subspace size *sz*.

   .. cpp:function:: void push(const torch::Tensor& F, const torch::Tensor& e)

      Append a (Fock matrix, error vector) pair.  If the queue exceeds
      *max_size*, the oldest pair is evicted.

   .. cpp:function:: torch::Tensor extrapolate()

      Solve the DIIS linear system (see :ref:`algorithms:DIIS Linear System`)
      and return the extrapolated Fock matrix.  Falls back to the most
      recent Fock matrix if the :math:`\mathbf{B}` matrix is singular
      (``at::linalg_solve`` fails).

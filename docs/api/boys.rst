``boys.hpp`` / ``boys.cpp`` — Boys Function
===========================================

.. contents:: On this page
   :depth: 2
   :local:

Overview
--------

Efficient, numerically stable computation of the Boys function
:math:`F_n(x)` required by the McMurchie-Davidson integral engine.

**Header:** ``include/boys.hpp``
**Source:** ``src/boys.cpp``

For the mathematical definition and derivation of the three evaluation
regimes see :ref:`algorithms:Boys Function`.

Functions
---------

.. cpp:function:: double boys_F(int n, double x)

   :param n: Order :math:`n \ge 0`.
   :param x: Argument :math:`x \ge 0`.
   :returns: :math:`F_n(x) = \int_0^1 t^{2n} e^{-xt^2}\,dt`.

   Convenience wrapper that allocates a temporary vector of size
   :math:`n+1`, calls :cpp:func:`boys_Fvec`, and returns the last element.
   Prefer :cpp:func:`boys_Fvec` in hot loops.

.. cpp:function:: void boys_Fvec(int nmax, double x, double* F)

   :param nmax: Maximum order; computes :math:`F_0, F_1, \ldots, F_{n_\text{max}}`.
   :param x: Argument :math:`x \ge 0`.
   :param F: **Output** array of size at least :math:`n_\text{max}+1`
             (caller-allocated).  ``F[n]`` receives :math:`F_n(x)`.

   Fills all Boys values in a single pass.  This is the primary entry
   point used by the integral engine, which requests
   :math:`n_\text{max} = \ell_1+\ell_2+\ell_3+\ell_4`
   per primitive quartet.

   The three evaluation strategies used internally:

   .. list-table::
      :header-rows: 1
      :widths: 20 80

      * - Range
        - Method
      * - :math:`x \le 10^{-14}`
        - Trivial limit: :math:`F_n = 1/(2n+1)`.
      * - :math:`0 < x < 1`
        - Taylor series
          :math:`\sum_{k\ge 0} (-x)^k / [k!\,(2n+2k+1)]`
          (converges in < 30 terms).
      * - :math:`x \ge 1`
        - :math:`F_0 = \tfrac{\sqrt\pi}{2\sqrt{x}}\operatorname{erf}(\sqrt{x})`,
          then upward recurrence
          :math:`F_{n+1} = [(2n+1)F_n - e^{-x}]/(2x)`;
          asymptotic fallback if result is non-finite.

Example
-------

.. code-block:: cpp

   #include "boys.hpp"

   // Single value
   double f2 = boys_F(2, 3.5);   // F_2(3.5)

   // Full vector (preferred in integral loops)
   int nmax = 4;
   double F[5];
   boys_Fvec(nmax, 1.23, F);
   // F[0] = F_0(1.23), ..., F[4] = F_4(1.23)

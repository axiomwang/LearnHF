#pragma once
#include <cmath>

// Compute the Boys function F_n(x) = ∫₀¹ t^{2n} exp(-x t²) dt
// Valid for n ≥ 0 and x ≥ 0.
// Uses incomplete Gamma: F_n(x) = (1/2) x^{-(n+1/2)} γ(n+1/2, x)
// but implemented via:
//   - Taylor series upward from F_0 for small x
//   - Downward recursion from asymptotic for large x
double boys_F(int n, double x);

// Compute all Boys values F_0(x) … F_nmax(x) in a single pass.
// Result stored in F[0..nmax] (caller allocates, size nmax+1).
void boys_Fvec(int nmax, double x, double* F);

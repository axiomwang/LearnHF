// McMurchie-Davidson integral engine
// References:
//   Helgaker, Jørgensen, Olsen: "Molecular Electronic-Structure Theory" (2000)
//   McMurchie & Davidson, J. Comput. Phys. 26 (1978) 218.
//   Obara & Saika, J. Chem. Phys. 84 (1986) 3963.

#include "integrals.hpp"
#include "boys.hpp"
#include <cmath>
#include <vector>
#include <array>
#include <cassert>
#include <stdexcept>

// ═══════════════════════════════════════════════════════════════════════════
//  McMurchie-Davidson E-coefficients
//  E_t^{ij}(α, β, X_AB) — Hermite Gaussian expansion
//
//  Recursion (HJO eq. 9.5.6):
//    E_t^{i+1,j} = 1/(2p) E_{t-1}^{ij} + X_PA E_t^{ij} + (t+1) E_{t+1}^{ij}
//    E_t^{i,j+1} = 1/(2p) E_{t-1}^{ij} + X_PB E_t^{ij} + (t+1) E_{t+1}^{ij}
//  Base: E_0^{00} = exp(-mu * X_AB^2),  mu = alpha*beta / (alpha+beta)
//  Convention: E_t^{ij} = 0 if t < 0 or t > i+j.
//
//  We store coefficients in a 3D array E[i][j][t] with i,j ≤ MAX_L, t ≤ i+j.
//  MAX_L = 8 covers up to H-type AOs (l=5 → 2*5=10 OK, but G=4 → 8 sufficient).
// ═══════════════════════════════════════════════════════════════════════════

static constexpr int MAX_LI = 9;  // max angular momentum index (0..8)

using ETable = double[MAX_LI][MAX_LI][2*MAX_LI];

static void fill_E(ETable E, double alpha, double beta, double XAB) {
    double p  = alpha + beta;
    double mu = alpha * beta / p;
    double XPA = beta / p * (-XAB);  // P - A with PA = beta/p * BA = -beta/p * XAB
    double XPB = alpha / p * XAB;    // P - B = alpha/p * AB = alpha/p * XAB
    // But convention: X_AB = A - B, so P = (α·A + β·B)/p
    // X_PA = P - A = -β/p * (A-B) = β/p * X_BA = -β/p * X_AB
    double base = std::exp(-mu * XAB * XAB);

    // Zero-fill
    for (int i = 0; i < MAX_LI; ++i)
        for (int j = 0; j < MAX_LI; ++j)
            for (int t = 0; t < 2*MAX_LI; ++t)
                E[i][j][t] = 0.0;

    E[0][0][0] = base;

    // Horizontal recursion to build E[i][0][t] for i = 1..MAX_LI-1
    for (int i = 0; i < MAX_LI-1; ++i) {
        for (int t = 0; t <= i+1; ++t) {
            double v = 0.0;
            if (t > 0 && i >= 0) v += 0.5/p * E[i][0][t-1];
            v += XPA * E[i][0][t];
            if (t+1 < 2*MAX_LI) v += (t+1) * E[i][0][t+1];
            E[i+1][0][t] = v;
        }
    }

    // Vertical recursion for j > 0
    for (int j = 0; j < MAX_LI-1; ++j) {
        for (int i = 0; i < MAX_LI; ++i) {
            if (i + j + 1 >= 2*MAX_LI) continue;
            for (int t = 0; t <= i+j+1; ++t) {
                double v = 0.0;
                if (t > 0) v += 0.5/p * E[i][j][t-1];
                v += XPB * E[i][j][t];
                if (t+1 < 2*MAX_LI) v += (t+1) * E[i][j][t+1];
                E[i][j+1][t] = v;
            }
        }
    }
}

// ═══════════════════════════════════════════════════════════════════════════
//  Auxiliary Hermite integrals R_{tuv}^n(p, RPC)
//  R_{tuv}^n = (-2p)^n ∂^t/∂P_x^t ∂^u/∂P_y^u ∂^v/∂P_z^v (1/r) evaluated
//  via Boys: R_{000}^n = F_n(p·RPC²) * (-2p)^n (for charge at C)
//
//  Recursion (HJO eq. 9.9.13):
//    R_{t+1,u,v}^n = t·R_{t-1,u,v}^{n+1} + X_PC·R_{t,u,v}^{n+1}
//  similarly for u, v.
//  Base: R_{000}^n = (-2p)^n F_n(p·RPC²)  [or with appropriate sign]
//
//  We use the standard convention from Helgaker et al., equations (9.9.12-14):
//    R_{t+1,u,v}^n = t R_{t-1,u,v}^{n+1} + (P_x - C_x) R_{tuv}^{n+1}
//  Base (n=0 auxiliary output):
//    R_{000}^n = (-2p)^n * F_n(p * |PC|²)
// ═══════════════════════════════════════════════════════════════════════════

// Max total degree t+u+v for R tensor; for ERI with 2 shells of max l=4 each:
// t+u+v ≤ 2*(4+4) = 16 → need order up to 16.
static constexpr int MAX_ORDER = 20;
// Dimensions: [order_n][t][u][v]
using RTable = double[MAX_ORDER][MAX_ORDER][MAX_ORDER][MAX_ORDER];

static void fill_R(double R[][MAX_ORDER][MAX_ORDER][MAX_ORDER],
                   double p, double Xpc, double Ypc, double Zpc, int nmax)
{
    double rpc2 = Xpc*Xpc + Ypc*Ypc + Zpc*Zpc;

    // Zero fill relevant portion
    for (int n = 0; n <= nmax; ++n)
        for (int t = 0; t <= nmax; ++t)
            for (int u = 0; u <= nmax-t; ++u)
                for (int v = 0; v <= nmax-t-u; ++v)
                    R[n][t][u][v] = 0.0;

    // Base: R[n][0][0][0] = (-2p)^n * F_n(p*rpc2)
    std::vector<double> Fv((size_t)(nmax+1));
    boys_Fvec(nmax, p * rpc2, Fv.data());
    double pow_m2p = 1.0;
    for (int n = 0; n <= nmax; ++n) {
        R[n][0][0][0] = pow_m2p * Fv[(size_t)n];
        pow_m2p *= -2.0 * p;
    }

    // Build up t+u+v for each n (decreasing n as we go up in t,u,v)
    // Recursion (Helgaker 9.9.13-15):
    //   R^n_{t+1,u,v} = t * R^{n+1}_{t-1,u,v} + Xpc * R^{n+1}_{t,u,v}
    //   R^n_{t,u+1,v} = u * R^{n+1}_{t,u-1,v} + Ypc * R^{n+1}_{t,u,v}
    //   R^n_{t,u,v+1} = v * R^{n+1}_{t,u,v-1} + Zpc * R^{n+1}_{t,u,v}
    for (int n = nmax; n >= 0; --n) {
        int nl = nmax - n;
        // Recurse in t
        for (int t = 0; t <= nl; ++t) {
            if (t > 0 && n+1 <= nmax) {
                for (int u = 0; u <= nl-t; ++u)
                    for (int v = 0; v <= nl-t-u; ++v) {
                        double val = Xpc * R[n+1][t-1][u][v];
                        if (t-1 > 0) val += (t-1) * R[n+1][t-2][u][v];
                        R[n][t][u][v] = val;
                    }
            }
        }
        // Recurse in u
        for (int t = 0; t <= nl; ++t)
            for (int u = 1; u <= nl-t; ++u)
                for (int v = 0; v <= nl-t-u; ++v) {
                    if (n+1 > nmax) continue;
                    double val = Ypc * R[n+1][t][u-1][v];
                    if (u-1 > 0) val += (u-1) * R[n+1][t][u-2][v];
                    R[n][t][u][v] = val;
                }
        // Recurse in v
        for (int t = 0; t <= nl; ++t)
            for (int u = 0; u <= nl-t; ++u)
                for (int v = 1; v <= nl-t-u; ++v) {
                    if (n+1 > nmax) continue;
                    double val = Zpc * R[n+1][t][u][v-1];
                    if (v-1 > 0) val += (v-1) * R[n+1][t][u][v-2];
                    R[n][t][u][v] = val;
                }
    }
}

// ═══════════════════════════════════════════════════════════════════════════
//  Primitive integrals (single pair of primitives)
// ═══════════════════════════════════════════════════════════════════════════

// 3D overlap of primitive Gaussians:
//   S(a,b,XAB) = Σ_t E_t^{ab} sqrt(π/p)
//   where a=(ax,ay,az), b=(bx,by,bz), XAB = (XA-XB component)
// Returns S_ab for all Cartesian components a,b.

// Compute sqrt(pi/p) * E_0^{a,b} for a 1D pair
static double overlap_1d(int ai, int bi, double alpha, double beta,
                          double XAB, const ETable& E) {
    double p = alpha + beta;
    return E[ai][bi][0] * std::sqrt(M_PI / p);
}

// Compute 3D primitive overlap integral
// ra=(ax,ay,az), rb=(bx,by,bz)
static double overlap_prim(double alpha, double beta,
                            const std::array<double,3>& A,
                            const std::array<double,3>& B,
                            const std::array<int,3>& a,
                            const std::array<int,3>& b)
{
    double p = alpha + beta;
    double S = 1.0;
    for (int k = 0; k < 3; ++k) {
        ETable E;
        fill_E(E, alpha, beta, A[k] - B[k]);
        S *= std::sqrt(M_PI / p) * E[a[k]][b[k]][0];
    }
    return S;
}

// Kinetic energy primitive:
// T_ab = alpha*(2*La+3)*S_ab - 2*alpha^2*S(a,b+2) - 0.5*Lb*(Lb-1)*S(a,b-2)
// where b+2 means increase one component of b by 2, summed over x,y,z.
// We use standard relation from Obara & Saika.
static double kinetic_prim(double alpha, double beta,
                            const std::array<double,3>& A,
                            const std::array<double,3>& B,
                            const std::array<int,3>& a,
                            const std::array<int,3>& b)
{
    double La = a[0]+a[1]+a[2];
    double Lb = b[0]+b[1]+b[2];
    // T_x + T_y + T_z where:
    // T_i = -0.5*[ Lb(Lb-1) * S(a, b-2e_i) - 2*(2Lb+1)*beta*S(a,b)
    //              + 4*beta^2 * S(a, b+2e_i) ]
    // This is the standard Obara-Saika result.
    // Equivalent form used here (from Szabo & Ostlund appendix + OSA):
    // T_ab = beta*(2Lb+3)*S(a,b) - 2*beta^2*S(a,b+2_xyz) - 0.5*Lb*(Lb-1)*...
    // Actually use the Helgaker et al. formula:
    // T = Σ_i [ β(β_i - 1/2)S(a, b - 2e_i) + (2β_i + 1)β S(a, b) - 2β² S(a, b+2e_i) ]
    // where β_i is b[i] for component i.
    // The correct formula (HJO 9.2.15):
    //  T_ab = β(2b_x+1)S_{a,b} - 2β² S_{a,b+2x} - (1/2)b_x(b_x-1)S_{a,b-2x} + [y] + [z]
    double T = 0.0;
    for (int k = 0; k < 3; ++k) {
        // Term 1: β*(2*b[k]+1) * S(a,b) [same b as input but scaled differently]
        // We need S for three different b-variations:
        // Compute full 3D overlap for each case

        auto bk = b;

        // S(a, b): already can compute
        double S0 = overlap_prim(alpha, beta, A, B, a, b);

        // S(a, b + 2*e_k)
        auto bp2 = b; bp2[k] += 2;
        double Sp2 = overlap_prim(alpha, beta, A, B, a, bp2);

        // S(a, b - 2*e_k) → 0 if b[k] < 2
        double Sm2 = 0.0;
        if (b[k] >= 2) {
            auto bm2 = b; bm2[k] -= 2;
            Sm2 = overlap_prim(alpha, beta, A, B, a, bm2);
        }

        T += beta * (2.0*b[k] + 1.0) * S0
           - 2.0 * beta * beta * Sp2
           - 0.5 * b[k] * (b[k]-1) * Sm2;
    }
    return T;
}

// Nuclear attraction primitive for charge Z at C:
// V_ab = -Z * 2π/p * Σ_{tuv} E_t^{ax ax} E_u^{ay ay} E_v^{az az} * R_{tuv}^0(p, PC)
static double nuclear_prim(double alpha, double beta,
                            const std::array<double,3>& A,
                            const std::array<double,3>& B,
                            const std::array<int,3>& a,
                            const std::array<int,3>& b,
                            double Z,
                            const std::array<double,3>& C)
{
    double p = alpha + beta;
    std::array<double,3> P = {
        (alpha*A[0] + beta*B[0]) / p,
        (alpha*A[1] + beta*B[1]) / p,
        (alpha*A[2] + beta*B[2]) / p
    };
    std::array<double,3> PC = { P[0]-C[0], P[1]-C[1], P[2]-C[2] };

    int lA = a[0]+a[1]+a[2];
    int lB = b[0]+b[1]+b[2];
    int nmax = lA + lB;

    // Fill E tables per dimension
    ETable Ex, Ey, Ez;
    fill_E(Ex, alpha, beta, A[0]-B[0]);
    fill_E(Ey, alpha, beta, A[1]-B[1]);
    fill_E(Ez, alpha, beta, A[2]-B[2]);

    // Allocate R table on stack (MAX_ORDER fixed size)
    static thread_local double R[MAX_ORDER][MAX_ORDER][MAX_ORDER][MAX_ORDER];
    fill_R(R, p, PC[0], PC[1], PC[2], nmax);

    double sum = 0.0;
    for (int t = 0; t <= a[0]+b[0]; ++t)
        for (int u = 0; u <= a[1]+b[1]; ++u)
            for (int v = 0; v <= a[2]+b[2]; ++v)
                sum += Ex[a[0]][b[0]][t] * Ey[a[1]][b[1]][u] * Ez[a[2]][b[2]][v]
                       * R[0][t][u][v];

    return -Z * 2.0 * M_PI / p * sum;
}

// ═══════════════════════════════════════════════════════════════════════════
//  2e- integral - ERI (μν|λσ) -- one shell quartet
//
//  (ab|cd) = (2π^{5/2}) / (p·q·sqrt(p+q)) * Σ Eabtuv * (-1)^{τυφ} Ecdτυφ * R(tuv+τυφ)
//  where p = α+β (for a,b), q = γ+δ (for c,d)
//  R is built at P-Q (center P from (a,b) pair, Q from (c,d) pair).
//
//  Returns a vector of integrals for all (a_cart × b_cart × c_cart × d_cart)
//  combinations in standard Cartesian order.
// ═══════════════════════════════════════════════════════════════════════════

// Compute all ERI values for a set of primitive pairs
static double eri_prim(
    double alpha, const std::array<double,3>& A, const std::array<int,3>& a,
    double beta,  const std::array<double,3>& B, const std::array<int,3>& b,
    double gamma, const std::array<double,3>& C, const std::array<int,3>& c,
    double delta, const std::array<double,3>& D, const std::array<int,3>& d)
{
    double p = alpha + beta;
    double q = gamma + delta;
    double rho = p * q / (p + q);

    std::array<double,3> P = {
        (alpha*A[0] + beta*B[0]) / p,
        (alpha*A[1] + beta*B[1]) / p,
        (alpha*A[2] + beta*B[2]) / p
    };
    std::array<double,3> Q_ = {
        (gamma*C[0] + delta*D[0]) / q,
        (gamma*C[1] + delta*D[1]) / q,
        (gamma*C[2] + delta*D[2]) / q
    };
    std::array<double,3> PQ = { P[0]-Q_[0], P[1]-Q_[1], P[2]-Q_[2] };

    int lA = a[0]+a[1]+a[2];
    int lB = b[0]+b[1]+b[2];
    int lC = c[0]+c[1]+c[2];
    int lD = d[0]+d[1]+d[2];
    int nmax = lA + lB + lC + lD;

    ETable Eabx, Eaby, Eabz;
    ETable Ecdx, Ecdy, Ecdz;
    fill_E(Eabx, alpha, beta, A[0]-B[0]);
    fill_E(Eaby, alpha, beta, A[1]-B[1]);
    fill_E(Eabz, alpha, beta, A[2]-B[2]);
    fill_E(Ecdx, gamma, delta, C[0]-D[0]);
    fill_E(Ecdy, gamma, delta, C[1]-D[1]);
    fill_E(Ecdz, gamma, delta, C[2]-D[2]);

    static thread_local double R[MAX_ORDER][MAX_ORDER][MAX_ORDER][MAX_ORDER];
    fill_R(R, rho, PQ[0], PQ[1], PQ[2], nmax);

    double sum = 0.0;
    for (int t = 0; t <= lA+lB; ++t)
      for (int u = 0; u <= lA+lB-t; ++u)
        for (int v = 0; v <= lA+lB-t-u; ++v) {
            double Eab = Eabx[a[0]][b[0]][t]
                       * Eaby[a[1]][b[1]][u]
                       * Eabz[a[2]][b[2]][v];
            if (std::abs(Eab) < 1e-20) continue;
            for (int tau = 0; tau <= lC+lD; ++tau)
              for (int ups = 0; ups <= lC+lD-tau; ++ups)
                for (int phi = 0; phi <= lC+lD-tau-ups; ++phi) {
                    double Ecd = Ecdx[c[0]][d[0]][tau]
                               * Ecdy[c[1]][d[1]][ups]
                               * Ecdz[c[2]][d[2]][phi];
                    if (std::abs(Ecd) < 1e-20) continue;
                    int sign = ((tau+ups+phi) % 2 == 0) ? 1 : -1;
                    sum += Eab * sign * Ecd * R[0][t+tau][u+ups][v+phi];
                }
        }

    double prefac = 2.0 * std::pow(M_PI, 2.5) / (p * q * std::sqrt(p + q));
    return prefac * sum;
}

// ═══════════════════════════════════════════════════════════════════════════
//  Public integral builders
// ═══════════════════════════════════════════════════════════════════════════

// Contraction normalization: the BSE stores coefficients such that
//   φ_μ(r) = Σ_k c_k * N_k * g(α_k, l, r-A)
// where N_k is the primitive norm already folded into prims.coeff.
// But the contracted shell still needs an overall contraction normalization
// N_contract such that <μ|μ> = 1.
// We include this normalization post-hoc by dividing by sqrt(S_μμ).

static void normalize_shell(Shell& sh) {
    // Compute the diagonal overlap <μ|μ> for the first Cartesian component
    // (all components of the same shell have the same norm by symmetry).
    int l = sh.l;
    std::array<int,3> a = {l, 0, 0};  // use first component
    const auto& A = sh.center;

    double S = 0.0;
    for (auto& p1 : sh.prims)
        for (auto& p2 : sh.prims)
            S += p1.coeff * p2.coeff
                 * overlap_prim(p1.exp, p2.exp, A, A, a, a);

    double N = 1.0 / std::sqrt(S);
    for (auto& p : sh.prims) p.coeff *= N;
}

torch::Tensor build_overlap(const BasisSet& bs_in) {
    BasisSet bs = bs_in;
    for (auto& sh : bs.shells) normalize_shell(sh);

    int n = bs.n_basis();
    auto S = torch::zeros({n, n}, torch::kDouble);
    auto Sa = S.accessor<double,2>();

    auto comps = [](int l) { return cart_components(l); };

    for (size_t i = 0; i < bs.shells.size(); ++i) {
        const Shell& si = bs.shells[i];
        int si0 = bs.shell_start((int)i);
        auto ci = comps(si.l);
        for (size_t j = 0; j <= i; ++j) {
            const Shell& sj = bs.shells[j];
            int sj0 = bs.shell_start((int)j);
            auto cj = comps(sj.l);
            for (int mi = 0; mi < si.n_cart(); ++mi) {
                for (int mj = 0; mj < sj.n_cart(); ++mj) {
                    double val = 0.0;
                    for (auto& pi : si.prims)
                        for (auto& pj : sj.prims)
                            val += pi.coeff * pj.coeff
                                   * overlap_prim(pi.exp, pj.exp,
                                                  si.center, sj.center,
                                                  ci[mi], cj[mj]);
                    Sa[si0+mi][sj0+mj] = val;
                    Sa[sj0+mj][si0+mi] = val;
                }
            }
        }
    }
    return S;
}

torch::Tensor build_kinetic(const BasisSet& bs_in) {
    BasisSet bs = bs_in;
    for (auto& sh : bs.shells) normalize_shell(sh);

    int n = bs.n_basis();
    auto T = torch::zeros({n, n}, torch::kDouble);
    auto Ta = T.accessor<double,2>();

    for (size_t i = 0; i < bs.shells.size(); ++i) {
        const Shell& si = bs.shells[i];
        int si0 = bs.shell_start((int)i);
        auto ci = cart_components(si.l);
        for (size_t j = 0; j <= i; ++j) {
            const Shell& sj = bs.shells[j];
            int sj0 = bs.shell_start((int)j);
            auto cj = cart_components(sj.l);
            for (int mi = 0; mi < si.n_cart(); ++mi) {
                for (int mj = 0; mj < sj.n_cart(); ++mj) {
                    double val = 0.0;
                    for (auto& pi : si.prims)
                        for (auto& pj : sj.prims)
                            val += pi.coeff * pj.coeff
                                   * kinetic_prim(pi.exp, pj.exp,
                                                  si.center, sj.center,
                                                  ci[mi], cj[mj]);
                    Ta[si0+mi][sj0+mj] = val;
                    Ta[sj0+mj][si0+mi] = val;
                }
            }
        }
    }
    return T;
}

torch::Tensor build_nuclear(const BasisSet& bs_in, const Molecule& mol) {
    BasisSet bs = bs_in;
    for (auto& sh : bs.shells) normalize_shell(sh);

    int n = bs.n_basis();
    auto V = torch::zeros({n, n}, torch::kDouble);
    auto Va = V.accessor<double,2>();

    for (size_t i = 0; i < bs.shells.size(); ++i) {
        const Shell& si = bs.shells[i];
        int si0 = bs.shell_start((int)i);
        auto ci = cart_components(si.l);
        for (size_t j = 0; j <= i; ++j) {
            const Shell& sj = bs.shells[j];
            int sj0 = bs.shell_start((int)j);
            auto cj = cart_components(sj.l);
            for (int mi = 0; mi < si.n_cart(); ++mi) {
                for (int mj = 0; mj < sj.n_cart(); ++mj) {
                    double val = 0.0;
                    for (const Atom& atm : mol.atoms) {
                        for (auto& pi : si.prims)
                            for (auto& pj : sj.prims)
                                val += pi.coeff * pj.coeff
                                       * nuclear_prim(pi.exp, pj.exp,
                                                      si.center, sj.center,
                                                      ci[mi], cj[mj],
                                                      (double)atm.Z, atm.coords);
                    }
                    Va[si0+mi][sj0+mj] = val;
                    Va[sj0+mj][si0+mi] = val;
                }
            }
        }
    }
    return V;
}

// ── Schwarz screening matrix ─────────────────────────────────────────────────

torch::Tensor build_schwarz(const BasisSet& bs_in) {
    BasisSet bs = bs_in;
    for (auto& sh : bs.shells) normalize_shell(sh);

    int n = bs.n_basis();
    auto Q = torch::zeros({n, n}, torch::kDouble);
    auto Qa = Q.accessor<double,2>();

    for (size_t i = 0; i < bs.shells.size(); ++i) {
        const Shell& si = bs.shells[i];
        int si0 = bs.shell_start((int)i);
        auto ci = cart_components(si.l);
        for (size_t j = 0; j <= i; ++j) {
            const Shell& sj = bs.shells[j];
            int sj0 = bs.shell_start((int)j);
            auto cj = cart_components(sj.l);
            for (int mi = 0; mi < si.n_cart(); ++mi) {
                for (int mj = 0; mj < sj.n_cart(); ++mj) {
                    double val = 0.0;
                    for (auto& pi : si.prims) for (auto& pj : sj.prims)
                        for (auto& pk : si.prims) for (auto& pl : sj.prims)
                            val += pi.coeff * pj.coeff * pk.coeff * pl.coeff
                                   * eri_prim(pi.exp, si.center, ci[mi],
                                              pj.exp, sj.center, cj[mj],
                                              pk.exp, si.center, ci[mi],
                                              pl.exp, sj.center, cj[mj]);
                    double q = (val > 0.0) ? std::sqrt(val) : 0.0;
                    Qa[si0+mi][sj0+mj] = q;
                    Qa[sj0+mj][si0+mi] = q;
                }
            }
        }
    }
    return Q;
}

// ── Build J and K matrices (Direct SCF) ─────────────────────────────────────

std::pair<torch::Tensor, torch::Tensor>
build_JK(const BasisSet& bs_in, const torch::Tensor& P_in,
         const torch::Tensor& screen, double tol)
{
    BasisSet bs = bs_in;
    for (auto& sh : bs.shells) normalize_shell(sh);

    int n = bs.n_basis();
    auto J = torch::zeros({n, n}, torch::kDouble);
    auto K = torch::zeros({n, n}, torch::kDouble);
    auto Ja = J.accessor<double,2>();
    auto Ka = K.accessor<double,2>();
    auto Pa = P_in.accessor<double,2>();
    auto Qa = screen.accessor<double,2>();

    int nsh = (int)bs.shells.size();

    for (int s1 = 0; s1 < nsh; ++s1) {
        const Shell& sh1 = bs.shells[s1];
        int b1 = bs.shell_start(s1);
        auto c1 = cart_components(sh1.l);
        for (int s2 = 0; s2 <= s1; ++s2) {
            const Shell& sh2 = bs.shells[s2];
            int b2 = bs.shell_start(s2);
            auto c2 = cart_components(sh2.l);

            // Screen: need max Q[mu,nu] over shell pair (s1,s2)
            double Qmn_max = 0.0;
            for (int m1 = 0; m1 < sh1.n_cart(); ++m1)
                for (int m2 = 0; m2 < sh2.n_cart(); ++m2)
                    Qmn_max = std::max(Qmn_max, std::abs(Qa[b1+m1][b2+m2]));

            for (int s3 = 0; s3 < nsh; ++s3) {
                const Shell& sh3 = bs.shells[s3];
                int b3 = bs.shell_start(s3);
                auto c3 = cart_components(sh3.l);
                for (int s4 = 0; s4 <= s3; ++s4) {
                    const Shell& sh4 = bs.shells[s4];
                    int b4 = bs.shell_start(s4);
                    auto c4 = cart_components(sh4.l);

                    double Qls_max = 0.0;
                    for (int m3 = 0; m3 < sh3.n_cart(); ++m3)
                        for (int m4 = 0; m4 < sh4.n_cart(); ++m4)
                            Qls_max = std::max(Qls_max, std::abs(Qa[b3+m3][b4+m4]));

                    if (Qmn_max * Qls_max < tol) continue;

                    // Compute ERI for this shell quartet
                    for (int m1 = 0; m1 < sh1.n_cart(); ++m1)
                      for (int m2 = 0; m2 < sh2.n_cart(); ++m2)
                        for (int m3 = 0; m3 < sh3.n_cart(); ++m3)
                          for (int m4 = 0; m4 < sh4.n_cart(); ++m4) {
                            // Cauchy-Schwarz per basis function pair
                            if (Qa[b1+m1][b2+m2] * Qa[b3+m3][b4+m4] < tol) continue;

                            double eri = 0.0;
                            for (auto& p1 : sh1.prims)
                              for (auto& p2 : sh2.prims)
                                for (auto& p3 : sh3.prims)
                                  for (auto& p4 : sh4.prims)
                                    eri += p1.coeff * p2.coeff * p3.coeff * p4.coeff
                                           * eri_prim(p1.exp, sh1.center, c1[m1],
                                                      p2.exp, sh2.center, c2[m2],
                                                      p3.exp, sh3.center, c3[m3],
                                                      p4.exp, sh4.center, c4[m4]);

                            int mu = b1+m1, nu = b2+m2, lm = b3+m3, sg = b4+m4;
                            // s2 ≤ s1 factor (shell degeneracy)
                            double d12 = (s1 == s2) ? 1.0 : 2.0;
                            double d34 = (s3 == s4) ? 1.0 : 2.0;

                            // Coulomb: J[μν] += P[λσ] * (μν|λσ)
                            double Pls  = Pa[lm][sg] * d34;
                            Ja[mu][nu] += Pls * eri;
                            if (mu != nu) Ja[nu][mu] += Pls * eri;

                            // Exchange: K[μλ] += P[νσ] * (μν|λσ)
                            double Pns  = Pa[nu][sg] * d34;
                            Ka[mu][lm] += 0.5 * Pns * eri * d12;
                            if (mu != lm) Ka[lm][mu] += 0.5 * Pns * eri * d12;

                            double Pnl  = Pa[nu][lm] * d34;
                            Ka[mu][sg] += 0.5 * Pnl * eri * d12;
                            if (mu != sg) Ka[sg][mu] += 0.5 * Pnl * eri * d12;
                        }
                }
            }
        }
    }
    return {J, K};
}

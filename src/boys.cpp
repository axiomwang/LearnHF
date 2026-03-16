#include "boys.hpp"

#include <cmath>
#include <vector>

namespace {

double boys_series_single(int n, double x) {
    double sum = 0.0;
    double term_factor = 1.0;  // (-x)^k / k!
    for (int k = 0; k < 120; ++k) {
        double term = term_factor / (2.0 * n + 2.0 * k + 1.0);
        sum += term;
        if (std::abs(term) < 1e-16) {
            break;
        }
        term_factor *= (-x) / (k + 1.0);
    }
    return sum;
}

}  // namespace

void boys_Fvec(int nmax, double x, double* F) {
    if (x <= 1e-14) {
        for (int n = 0; n <= nmax; ++n) {
            F[n] = 1.0 / (2.0 * n + 1.0);
        }
        return;
    }

    if (x < 1.0) {
        for (int n = 0; n <= nmax; ++n) {
            F[n] = boys_series_single(n, x);
        }
        return;
    }

    double sqrtx = std::sqrt(x);
    F[0] = 0.5 * std::sqrt(M_PI / x) * std::erf(sqrtx);

    double expx = std::exp(-x);
    for (int n = 0; n < nmax; ++n) {
        F[n + 1] = ((2.0 * n + 1.0) * F[n] - expx) / (2.0 * x);
        if (!std::isfinite(F[n + 1])) {
            double df = 1.0;
            for (int k = 1; k <= n + 1; ++k) {
                df *= (2.0 * k - 1.0);
            }
            F[n + 1] = df * std::sqrt(M_PI) /
                       (std::pow(2.0, n + 2.0) * std::pow(x, n + 1.5));
        }
    }
}

double boys_F(int n, double x) {
    std::vector<double> values(static_cast<size_t>(n + 1));
    boys_Fvec(n, x, values.data());
    return values[static_cast<size_t>(n)];
}

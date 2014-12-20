#include "numerics.h"
#include <algorithm>
#include <iostream>
#include <cassert>

namespace tbem {

/* Evaluate the chebyshev polynomials up through degree n = n_max
 * at the point x_hat in the range [-1, 1].
 */
std::vector<double> cheb_polys(double x_hat, int n_max) {
    std::vector<double> res(n_max + 1);
    res[0] = 1.0; 
    if (n_max == 0) {
        return res;
    }
    res[1] = x_hat; 
    if (n_max == 1) {
        return res;
    }
    for (int i = 2; i <= n_max; i++) {
        // The Chebyshev polynomial recurrence relation.
        res[i] = 2 * x_hat * res[i - 1] - res[i - 2];
    }
    return res;
}

/* Compute the Chebyshev anterpolation and interpolation operator.
 * See the Black-Box Fast Multipole Method paper by Fong and Darve 2009
 * for details.
 */
double s_n(double x_hat, double y_hat, unsigned int n) {
    assert(n >= 0);
    auto x_cheb = cheb_polys(x_hat, n - 1);
    auto y_cheb = cheb_polys(y_hat, n - 1);
    double result = 0.0;
    for (unsigned int i = 1; i < n; i++) {
        result += x_cheb[i] * y_cheb[i];
    }
    result *= (2.0 / n); 
    result += (1.0 / n);
    return result;
}

/* Same as s_n above but optimized. Mainly, the cheb_polys function is folded
 * in so there isn't a std::vector creation and deletion for every call*/
double s_n_fast(double x_hat, double y_hat, unsigned int n) {
    if (n == 1) {
        return 1.0;
    }
    if (n == 2) {
        return 0.5 + (x_hat * y_hat);
    }
    double result = x_hat * y_hat;
    double prev_x_hat = x_hat;
    double prev_prev_x_hat = 1.0;
    double prev_y_hat = y_hat;
    double prev_prev_y_hat = 1.0;
    for (unsigned int i = 2; i < n; i++) {
        // The Chebyshev polynomial recurrence relation.
        double cur_x_hat = 2 * x_hat * prev_x_hat - prev_prev_x_hat;
        double cur_y_hat = 2 * y_hat * prev_y_hat - prev_prev_y_hat;
        result += cur_x_hat * cur_y_hat;

        prev_prev_y_hat = prev_y_hat;
        prev_y_hat = cur_y_hat;
        prev_prev_x_hat = prev_x_hat;
        prev_x_hat = cur_x_hat;
    }
    result *= 2.0 / n;
    result += 1.0 / n;
    return result;
}

/* Compute the Chebyshev points of the first kind. These are the roots
 * of the Chebyshev polynomial of degree n.
 */
std::vector<double> cheb_pts_first_kind(unsigned int n_pts) {
    std::vector<double> x(n_pts);
    for (unsigned int i = 1; i <= n_pts; i++) {
        x[i - 1] = std::cos((2 * i - 1) * M_PI / (2 * n_pts));
    }
    return x;
}

} //END NAMESPACE tbem

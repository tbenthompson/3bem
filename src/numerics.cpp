#include "numerics.h"
#include <algorithm>
#include <iostream>
#include <cassert>


/* equivalent to range(min, max) in python */
std::vector<int> naturals(int min, int max) {
    std::vector<int> indices(max);
    std::iota(indices.begin(), indices.end(), min);
    return indices;
}

/* equivalent to range(0, max) in python */
std::vector<int> naturals(int max) {
    return naturals(0, max);
}

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

/* Compute the Double exponential (also called Tanh-Sinh) 
 * quadrature rule with n points.
 */
QuadratureRule double_exp(int n, double h) {
    QuadratureRule retval;
    for (int i = -n; i <= n; i++) {
        const double sinhterm = 0.5 * M_PI * std::sinh(i * h);
        const double coshsinh = std::cosh(sinhterm);
        const double x = std::tanh(sinhterm);
        const double w = h * 0.5 * M_PI * std::cosh(i * h) / (coshsinh * coshsinh);
        retval.push_back(std::make_pair(x, w));
    }
    return retval;
}

/* Evaluate legendre polynomials P_{n}(x) and P_{n - 1}(x).
 * This is a helper function for the Gaussian quadrature algorithm.
 */
std::pair<double, double> legendre_and_n_minus_1(unsigned int n, 
                                                 double x) {
    double p_cur = 1.;
    double p_last = 0.;
    double p_temp;
    for (unsigned int j=0; j<n; ++j)
    {
        p_temp = p_last;
        p_last = p_cur;
        p_cur = ((2.*j+1.)*x*p_last-j*p_temp)/(j+1);
    }
    return std::make_pair(p_cur, p_last);
}

/* Compute the gaussian quadrature rule with n points.
 * The algorithm used first computes an initial analytic approximation
 * to the roots of the Legendre polynomial of degree n.
 * Then, the analytic approximation is refined using Newton's method.
 * This should work for rules up to about n = 1000.
 */
QuadratureRule gauss(unsigned int n) {
    QuadratureRule retval(n);
    const double tolerance = 1e-14;
    //Because gaussian quadrature rules are symmetric, I only 
    const unsigned int m = (n+1)/2;
    for (unsigned int i=0;i < m;i++)
    {
        // Initial guess.
        double x = std::cos(M_PI * (i + (3.0/4.0)) / (n + 0.5));

        double dp = 0;
        double dx = 10;

        // Perform newton iterations until the quadrature points
        // have converged.
        while (std::fabs(dx) > tolerance)
        {
            std::pair<double, double> p_n_and_nm1 =
                legendre_and_n_minus_1(n, x);
            double p_n = p_n_and_nm1.first;
            double p_nm1 = p_n_and_nm1.second;
            dp = (n + 1) * (x * p_n - p_nm1) / (x * x - 1);
            dx = p_n / dp;
            x = x - dx;
        }

        double w = 2 * (n + 1) * (n + 1) / (n * n * (1 - x * x) * dp * dp);
        retval[i] = std::make_pair(-x, w);
        retval[n - i - 1] = std::make_pair(x, w);
    }

    return retval;
}

/* A nonlinear mapping that smooths out a singular or near-singular integral.
 *
 * n is the number of points of the underlying Gauss quadrature rule.
 * x0 is the location in reference space of the singularity or most singular point.
 * q specifies how smooth to make the integrand. It must be an odd positive integer.
 *
 * This quadrature rule is presented in section 4.1 of the paper:
 * Aimi, a., and M. Diligenti. “Numerical Integration in 3D Galerkin BEM Solution of HBIEs.” Computational Mechanics 28, no. 3–4 (April 01, 2002): 233–49. doi:10.1007/s00466-001-0284-9.
 * To derive the version used below:
 */
QuadratureRule diligenti_mapping(unsigned int n, double x0, int q) {
    double x0_mapped = from_11_to_01(x0);
    double inv_q = 1.0 / q;
    double qth_root_x0 = pow(x0_mapped, inv_q);
    double qth_root_1mx0 = pow(1 - x0_mapped, inv_q);
    double zs = qth_root_x0 / (qth_root_x0 + qth_root_1mx0);
    double delta = (1 - x0_mapped) / pow((1 - zs), q);

    auto underlying = gauss(n);
    QuadratureRule retval(n);
    for (unsigned int i = 0; i < n; i++) {
        double z_hat = from_11_to_01(underlying[i].first) - zs;
        retval[i].first = from_01_to_11(x0_mapped + delta * pow(z_hat, q));
        retval[i].second = underlying[i].second * q * delta * pow(z_hat, q - 1);
    }
    return retval;
}

/* A helper function for testing the quadrature rules. Accepts a function
 * and integrates it according to the specified quadrature rule.
 */
double integrate(QuadratureRule& qr, std::function<double (double)> fnc) {
    double integral_val = 0;
    for (auto xw: qr) {
        integral_val += xw.second * fnc(xw.first);
    }
    return integral_val;
}

/* Another helper function, but for 2D integration. 
 * TODO: Make 1D and 2D quadrature more similar.
 */
double integrate(QuadratureRule2D& qr, std::function<double (double,double)> fnc) {
    double integral_val = 0;
    for (unsigned int i = 0; i < qr.x_hat.size(); i++) {
        integral_val += qr.weights[i] * fnc(qr.x_hat[i], qr.y_hat[i]);
    }
    return integral_val;
}


/* Produce a 2D tensor product quadrature rule from the product of two
 * one quadrature rule. 
 */
QuadratureRule2D tensor_product(QuadratureRule xq, QuadratureRule yq) {
    unsigned int xn = xq.size();
    unsigned int yn = yq.size();
    QuadratureRule2D retval(xn * yn);
    for(unsigned int i = 0; i < xn; i++) {
        for(unsigned int j = 0; j < yn; j++) {
            int idx_2d = i * yn + j;
            retval.x_hat[idx_2d] = xq[i].first;
            retval.y_hat[idx_2d] = yq[j].first;
            retval.weights[idx_2d] = xq[i].second * yq[j].second;
        }
    }
    return retval;
}

/* Converts the square [-1,1]x[-1,1] to the triangle (0,0)-(1,0)-(0,1)
 */
QuadratureRule2D square_to_tri(QuadratureRule2D square_quad) {
    // assert(square_quad.x_hat.size() == square_quad.y_hat.size())
    // assert(square_quad.y_hat.size() == square_quad.weights.size())

    unsigned int nq = square_quad.x_hat.size();
    QuadratureRule2D retval(nq);
    for (unsigned int i = 0; i < nq; i++) {
        double x_01 = from_11_to_01(square_quad.x_hat[i]);
        double y_01 = from_11_to_01(square_quad.y_hat[i]);
        retval.x_hat[i] = x_01 * (1 - y_01);
        retval.y_hat[i] = y_01;
        retval.weights[i] = (square_quad.weights[i] / 4.0) * (1 - y_01);
    }
    return retval;
}

/* Produces a 2D tensor product gaussian quadrature rule. The number of gauss
 * points in each dimension is the same.
 */
QuadratureRule2D tensor_gauss(int n_pts) {
    auto g1d = gauss(n_pts);
    return tensor_product(g1d, g1d);
}

/* Produces a 2D tensor product gaussian quadrature rule mapped into the unit
 * triangle (0,0)-(1,0)-(0,1).
 */
QuadratureRule2D tri_gauss(int n_pts) {
    return square_to_tri(tensor_gauss(n_pts));
}

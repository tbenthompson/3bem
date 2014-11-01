#include "quadrature.h"
#include "numerics.h"
#include <cmath>

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
    for (unsigned int j=0; j<n; ++j)
    {
        double p_temp = p_last;
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

/* Produces a 2D tensor product double exponential rule. */
QuadratureRule2D tensor_double_exp(int n_pts, double h) {
    auto de1d = double_exp(n_pts, h);
    return tensor_product(de1d, de1d);
}

/* Produces a 2D tensor product double exponetial rule mapped into the unit
 * triangle (0,0)-(1,0)-(0,1).
 */
QuadratureRule2D tri_double_exp(int n_pts, double h) {
    return square_to_tri(tensor_double_exp(n_pts, h));
}

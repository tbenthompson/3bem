#include "quadrature.h"
#include "numerics.h"
#include <cmath>

namespace tbem {

/* Evaluate legendre polynomials P_{n}(x) and P_{n - 1}(x).
 * This is a helper function for the Gaussian quadrature algorithm.
 */
std::pair<double, double> legendre_and_n_minus_1(size_t n, 
                                                 double x) {
    double p_cur = 1.;
    double p_last = 0.;
    for (size_t j=0; j<n; ++j)
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
QuadRule<1> gauss(size_t n) {
    assert(n > 0);
    std::vector<Vec2<double>> points(n);
    const double tolerance = 1e-14;
    //Because gaussian quadrature rules are symmetric, I only compute half of
    //the points and then mirror across x = 0.
    const size_t m = (n+1)/2;
    for (size_t i = 0; i < m; i++)
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
        points[i] = {-x, w};
        points[n - i - 1] = {x, w};
    }

    QuadRule<1> retval;
    for (std::size_t i = 0; i < n; i++) {
        retval.push_back(QuadPt<1>{points[i][0], points[i][1]});
    }

    return retval;
}

QuadRule<1> sinh_transform(const QuadRule<1>& gauss_rule, double a,
    double b, bool iterated_sinh) 
{
    auto mu_0 = 0.5 * (std::asinh((1.0 + a) / b) + std::asinh((1.0 - a) / b));
    auto eta_0 = 0.5 * (std::asinh((1.0 + a) / b) - std::asinh((1.0 - a) / b));
    auto start_q = gauss_rule;
    if (iterated_sinh) {
        start_q.clear();
        double a_1 = eta_0 / mu_0;
        double b_1 = M_PI / (2 * mu_0);
        double mu_1 = 0.5 * 
            (std::asinh((1.0 + a_1) / b_1) + std::asinh((1.0 - a_1) / b_1));
        double eta_1 = 0.5 * 
            (std::asinh((1.0 + a_1) / b_1) - std::asinh((1.0 - a_1) / b_1));
        for (size_t i = 0; i < gauss_rule.size(); i++) {
            double u = gauss_rule[i].x_hat[0];
            double u_w = gauss_rule[i].w;
            double s = a_1 + b_1 * std::sinh(mu_1 * u - eta_1);
            double jacobian = b_1 * mu_1 * std::cosh(mu_1 * u - eta_1);
            double s_w = u_w * jacobian;
            start_q.push_back({{s}, s_w});
        }
    }
    std::vector<QuadPt<1>> q_pts;
    for (size_t i = 0; i < start_q.size(); i++) {
        auto s = start_q[i].x_hat[0];
        auto x = a + b * std::sinh(mu_0 * s - eta_0);
        auto jacobian = b * mu_0 * std::cosh(mu_0 * s - eta_0);
        auto w = start_q[i].w * jacobian; 
        q_pts.push_back({x, w});
    }
    return q_pts;
}

/* Produce a 2D tensor product quadrature rule from the product of two
 * one quadrature rule. 
 */
QuadRule<2> tensor_product(QuadRule<1> xq, QuadRule<1> yq) {
    size_t xn = xq.size();
    size_t yn = yq.size();
    QuadRule<2> retval;
    for(size_t i = 0; i < xn; i++) {
        for(size_t j = 0; j < yn; j++) {
            retval.push_back(QuadPt<2>{
                {xq[i].x_hat[0], yq[j].x_hat[0]},
                xq[i].w * yq[j].w
            });
        }
    }
    return retval;
}

/* Converts the square [-1,1]x[-1,1] to the triangle (0,0)-(1,0)-(0,1)
 */
QuadRule<2> square_to_tri(QuadRule<2> square_quad) {
    size_t nq = square_quad.size();
    QuadRule<2> retval;
    for (size_t i = 0; i < nq; i++) {
        double x_01 = from_11_to_01(square_quad[i].x_hat[0]);
        double y_01 = from_11_to_01(square_quad[i].x_hat[1]);
        double w = square_quad[i].w;
        retval.push_back(QuadPt<2>{
            {x_01 * (1 - y_01), y_01},
            (w / 4.0) * (1 - y_01)
        });
    }
    return retval;
}

/* Produces a 2D tensor product gaussian quadrature rule. The number of gauss
 * points in each dimension is the same.
 */
QuadRule<2> tensor_gauss(size_t n_pts) {
    auto g1d = gauss(n_pts);
    return tensor_product(g1d, g1d);
}

/* Produces a 2D tensor product gaussian quadrature rule mapped into the unit
 * triangle (0,0)-(1,0)-(0,1).
 */
QuadRule<2> tri_gauss(size_t n_pts) {
    return square_to_tri(tensor_gauss(n_pts));
}

} //END NAMESPACE tbem

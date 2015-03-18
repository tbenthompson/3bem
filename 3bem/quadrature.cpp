#include "quadrature.h"
#include "numerics.h"
#include <cmath>
#include <cassert>
#include <map>

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
    // assert(square_quad.x_hat.size() == square_quad.y_hat.size())
    // assert(square_quad.y_hat.size() == square_quad.weights.size())

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
QuadRule<2> tensor_gauss(int n_pts) {
    auto g1d = gauss(n_pts);
    return tensor_product(g1d, g1d);
}

/* Produces a 2D tensor product gaussian quadrature rule mapped into the unit
 * triangle (0,0)-(1,0)-(0,1).
 */
QuadRule<2> tri_gauss(int n_pts) {
    return square_to_tri(tensor_gauss(n_pts));
}


std::vector<double> get_singular_steps(int n_steps) {
    static constexpr double initial_dist = 1.0;
    std::vector<double> dist(n_steps);
    for (int nf = 0; nf < n_steps; nf++) {
        dist[nf] = initial_dist / (std::pow(2, nf));
    }
    return dist;
}

template <size_t dim>
QuadStrategy<dim>::QuadStrategy(int obs_order):
    QuadStrategy(obs_order, obs_order, 8, 3.0, 1e-4)
{}

template <>
QuadStrategy<2>::QuadStrategy(int obs_order, int src_far_order, 
                           int n_singular_steps, double far_threshold,
                           double near_tol):
    obs_quad(gauss(obs_order)),
    src_far_quad(gauss(src_far_order)),
    far_threshold(far_threshold),
    n_singular_steps(n_singular_steps),
    singular_steps(get_singular_steps(n_singular_steps)),
    near_tol(near_tol)
{}

template <>
QuadStrategy<3>::QuadStrategy(int obs_order, int src_far_order, 
                           int n_singular_steps, double far_threshold,
                           double near_tol):
    obs_quad(tri_gauss(obs_order)),
    src_far_quad(tri_gauss(src_far_order)),
    far_threshold(far_threshold),
    n_singular_steps(n_singular_steps),
    singular_steps(get_singular_steps(n_singular_steps)),
    near_tol(near_tol)
{}

template class QuadStrategy<2>;
template class QuadStrategy<3>;

} //END NAMESPACE tbem

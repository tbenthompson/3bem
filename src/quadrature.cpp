#include "quadrature.h"
#include "numerics.h"
#include <cmath>
#include <cassert>

/* A helper function for integrating a given function using a quadrature rule.
 * Via templating, can be used with 1D, 2D, double, Vec3<double> quadrature.
 */
template <typename T, int dim>
T integrate(const std::vector<QuadPt<dim>>& qr, 
            const std::function<T(std::array<double,dim>)>& fnc) {
    T integral_val = qr[0].w * fnc(qr[0].x_hat);;
    for (unsigned int i = 1; i < qr.size(); i++) {
        integral_val += qr[i].w * fnc(qr[i].x_hat);
    }
    return integral_val;
}

//Explicitly instantiate the only reasonable options for the templated
//integrate function.
template double integrate(const QuadRule<1>&, 
        const std::function<double(std::array<double,1>)>&);
template double integrate(const QuadRule<2>&,
        const std::function<double(std::array<double,2>)>&);
template Vec3<double> integrate(const QuadRule<1>&, 
        const std::function<Vec3<double>(std::array<double,1>)>&);
template Vec3<double> integrate(const QuadRule<2>&,
        const std::function<Vec3<double>(std::array<double,2>)>&);

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
QuadRule<1> gauss(unsigned int n) {
    assert(n > 0);
    QuadRule<1> retval(n);
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
        retval[i] = {{-x},w};
        retval[n - i - 1] = {{x},w};
    }

    return retval;
}

/* Produce a 2D tensor product quadrature rule from the product of two
 * one quadrature rule. 
 */
QuadRule<2> tensor_product(QuadRule<1> xq, QuadRule<1> yq) {
    unsigned int xn = xq.size();
    unsigned int yn = yq.size();
    QuadRule<2> retval(xn * yn);
    for(unsigned int i = 0; i < xn; i++) {
        for(unsigned int j = 0; j < yn; j++) {
            int idx_2d = i * yn + j;
            retval[idx_2d] = {
                {xq[i].x_hat[0], yq[j].x_hat[0]},
                xq[i].w * yq[j].w
            };
        }
    }
    return retval;
}

/* Converts the square [-1,1]x[-1,1] to the triangle (0,0)-(1,0)-(0,1)
 */
QuadRule<2> square_to_tri(QuadRule<2> square_quad) {
    // assert(square_quad.x_hat.size() == square_quad.y_hat.size())
    // assert(square_quad.y_hat.size() == square_quad.weights.size())

    unsigned int nq = square_quad.size();
    QuadRule<2> retval(nq);
    for (unsigned int i = 0; i < nq; i++) {
        double x_01 = from_11_to_01(square_quad[i].x_hat[0]);
        double y_01 = from_11_to_01(square_quad[i].x_hat[1]);
        double w = square_quad[i].w;
        retval[i] = {
            {x_01 * (1 - y_01), y_01},
            (w / 4.0) * (1 - y_01)
        };
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

template <int dim>
QuadStrategy<dim>::QuadStrategy(int obs_order):
    QuadStrategy(obs_order, obs_order, obs_order * 2, 8, 3.0, 1e-4)
{}

template <>
QuadStrategy<2>::QuadStrategy(int obs_order, int src_far_order, int src_near_order,
                           int n_singular_steps, double far_threshold,
                           double singular_tol):
    obs_quad(gauss(obs_order)),
    src_far_quad(gauss(src_far_order)),
    src_near_quad(gauss(src_near_order)),
    far_threshold(far_threshold),
    n_singular_steps(n_singular_steps),
    singular_steps(get_singular_steps(n_singular_steps)),
    singular_tol(singular_tol)
{}

template <>
QuadStrategy<3>::QuadStrategy(int obs_order, int src_far_order, int src_near_order,
                           int n_singular_steps, double far_threshold,
                           double singular_tol):
    obs_quad(tri_gauss(obs_order)),
    src_far_quad(tri_gauss(src_far_order)),
    src_near_quad(tri_gauss(src_near_order)),
    far_threshold(far_threshold),
    n_singular_steps(n_singular_steps),
    singular_steps(get_singular_steps(n_singular_steps)),
    singular_tol(singular_tol)
{}

template class QuadStrategy<2>;
template class QuadStrategy<3>;

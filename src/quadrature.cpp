#include "quadrature.h"
#include "numerics.h"
#include <cmath>

/* A helper function for integrating a given function using a quadrature rule.
 * Via templating, can be used with 1D, 2D, double, Vec3<double> quadrature.
 */
template <typename T, int dim>
T integrate(const std::vector<QuadPt<dim>>& qr, 
            std::function<T(std::array<double,dim>)> fnc) {
    T integral_val = qr[0].w * fnc(qr[0].x_hat);;
    for (unsigned int i = 1; i < qr.size(); i++) {
        integral_val += qr[i].w * fnc(qr[i].x_hat);
    }
    return integral_val;
}

//Explicitly instantiate the only reasonable options for the templated
//integrate function.
template double integrate(const QuadRule1d&, 
        std::function<double(std::array<double,1>)>);
template double integrate(const QuadRule2d&,
        std::function<double(std::array<double,2>)>);
template Vec3<double> integrate(const QuadRule1d&, 
        std::function<Vec3<double>(std::array<double,1>)>);
template Vec3<double> integrate(const QuadRule2d&,
        std::function<Vec3<double>(std::array<double,2>)>);


/* Compute the Double exponential (also called Tanh-Sinh) 
 * quadrature rule with 2n + 1 points.
 */
QuadRule1d double_exp(int n, double h) {
    QuadRule1d retval;
    for (int i = -n; i <= n; i++) {
        const double sinhterm = 0.5 * M_PI * std::sinh(i * h);
        const double coshsinh = std::cosh(sinhterm);
        const double x = std::tanh(sinhterm);
        const double w = h * 0.5 * M_PI * std::cosh(i * h) / (coshsinh * coshsinh);
        retval.push_back({{x},w});
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
QuadRule1d gauss(unsigned int n) {
    QuadRule1d retval(n);
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
QuadRule1d diligenti_mapping(unsigned int n, double x0, int q) {
    double x0_mapped = from_11_to_01(x0);
    double inv_q = 1.0 / q;
    double qth_root_x0 = pow(x0_mapped, inv_q);
    double qth_root_1mx0 = pow(1 - x0_mapped, inv_q);
    double zs = qth_root_x0 / (qth_root_x0 + qth_root_1mx0);
    double delta = (1 - x0_mapped) / pow((1 - zs), q);

    auto underlying = gauss(n);
    QuadRule1d retval(n);
    for (unsigned int i = 0; i < n; i++) {
        double z_hat = from_11_to_01(underlying[i].x_hat[0]) - zs;
        double x = from_01_to_11(x0_mapped + delta * pow(z_hat, q));
        double w = underlying[i].w * q * delta * pow(z_hat, q - 1);
        retval[i] = {{x},w};
    }
    return retval;
}


/* Produce a 2D tensor product quadrature rule from the product of two
 * one quadrature rule. 
 */
QuadRule2d tensor_product(QuadRule1d xq, QuadRule1d yq) {
    unsigned int xn = xq.size();
    unsigned int yn = yq.size();
    QuadRule2d retval(xn * yn);
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
QuadRule2d square_to_tri(QuadRule2d square_quad) {
    // assert(square_quad.x_hat.size() == square_quad.y_hat.size())
    // assert(square_quad.y_hat.size() == square_quad.weights.size())

    unsigned int nq = square_quad.size();
    QuadRule2d retval(nq);
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
QuadRule2d tensor_gauss(int n_pts) {
    auto g1d = gauss(n_pts);
    return tensor_product(g1d, g1d);
}

/* Produces a 2D tensor product gaussian quadrature rule mapped into the unit
 * triangle (0,0)-(1,0)-(0,1).
 */
QuadRule2d tri_gauss(int n_pts) {
    return square_to_tri(tensor_gauss(n_pts));
}

/* Produces a 2D tensor product double exponential rule. */
QuadRule2d tensor_double_exp(int n_pts, double h) {
    auto de1d = double_exp(n_pts, h);
    return tensor_product(de1d, de1d);
}

/* Produces a 2D tensor product double exponetial rule mapped into the unit
 * triangle (0,0)-(1,0)-(0,1).
 */
QuadRule2d tri_double_exp(int n_pts, double h) {
    return square_to_tri(tensor_double_exp(n_pts, h));
}

QuadRule2d tri_double_exp(int n_pts) {
    double h = 0.6 / std::log(n_pts);
    return tri_double_exp(n_pts, h);
}

std::vector<double> get_singular_steps(int n_steps) {
    static constexpr double initial_dist = 1.0;
    std::vector<double> dist(n_steps);
    for (int nf = 0; nf < n_steps; nf++) {
        dist[nf] = initial_dist / (std::pow(2, nf));
    }
    return dist;
}

QuadStrategy::QuadStrategy(int obs_order):
    QuadStrategy(obs_order, obs_order, obs_order * 2, 4, 3.0, 1e-2)
{}

QuadStrategy::QuadStrategy(int obs_order, int src_far_order, int src_near_order,
                           int n_singular_steps, double far_threshold,
                           double singular_tol):
    obs_quad(tri_gauss(obs_order)),
    src_far_quad(tri_gauss(src_far_order)),
    src_near_quad(tri_double_exp(src_near_order)),
    far_threshold(far_threshold),
    n_singular_steps(n_singular_steps),
    singular_steps(get_singular_steps(n_singular_steps)),
    singular_tol(singular_tol)
{}

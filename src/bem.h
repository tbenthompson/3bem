#ifndef __AAAAAAAAA_BEM_H
#define __AAAAAAAAA_BEM_H

#include <functional>
#include <array>
#include <vector>
#include "vec.h"
#include "numerics.h"
#include "quadrature.h"

class Mesh;
template <typename T, int M>
class Taylor;

const int taylor_degree = 6;

template <typename T>
using GenericKernel = std::function<T
    (const T&, const Vec3<T>&, const Vec3<double>&, const Vec3<double>&)>;

typedef GenericKernel<double> Kernel;
typedef GenericKernel<Taylor<double,taylor_degree>> TaylorKernel;

/* Temporary data store for the evaluation of nearfield integrals. */
class NearEval {
public:
    NearEval(int n_steps);

    static constexpr double initial_dist = 1.0;

    const int n_steps;
    //TODO: const (remove state!)
    std::vector<QuadratureRule2D> quad;
    std::vector<double> dist;
};


class FaceInfo {
public:
    FaceInfo(const Mesh& mesh, int face_index);
    
    const std::array<int,3>& face;
    const std::array<Vec3<double>,3> corners;
    const Vec3<double> unscaled_normal;
    const double area;
    const double jacobian;
    const Vec3<double> normal;
};

template <typename T>
struct SrcPointInfo {
    SrcPointInfo(const QuadratureRule2D& quad_rule,
              const GenericKernel<T>& kernel,
              const FaceInfo& face,
              const Vec3<T>& obs_loc,
              const Vec3<double>& obs_n,
              int q_index) {
        const double x_hat = quad_rule.x_hat[q_index];
        const double y_hat = quad_rule.y_hat[q_index];
        const double q_wt = quad_rule.weights[q_index];
        basis = linear_basis(x_hat, y_hat);

        const auto src_pt = ref_to_real(x_hat, y_hat, face.corners);
        const Vec3<T> d = {
            src_pt[0] - obs_loc[0], 
            src_pt[1] - obs_loc[1],
            src_pt[2] - obs_loc[2]
        };
        const auto r2 = hypot2(d);
        const auto kernel_val = kernel(r2, d, face.normal, obs_n);

        weighted_kernel = kernel_val * q_wt * face.jacobian;
    }
    
    Vec3<double> basis;
    T weighted_kernel;
};

/* Perform the richardson extrapolation for the nearfield quadrature. 
 */
//TODO: Use the diligenti mapping quadrature for the nearly 
//singular quadratures required.
template <typename T>
T richardson_step(const std::vector<T>& values) {
    assert(values.size() > 1);
    int n_steps = values.size();
    std::vector<T> last_level = values;
    std::vector<T> this_level;
    int error_order = 1;
    for (int m = 1; m < n_steps; m++) {
        this_level.resize(n_steps - m);
        for (int i = 0; i < n_steps - m; i++) {
            const double mult = pow(2, error_order);
            const double factor = 1.0 / (mult - 1.0);
            T low = last_level[i];
            T high = last_level[i + 1];
            T moreacc = factor * (mult * high - low);
            this_level[i] = moreacc;
        }
        error_order++;
        last_level = this_level;
    }
    // std::cout << this_level[0] << std::endl;
    return this_level[0];
}

template <typename T>
Vec3<T> basis_integrals(const QuadratureRule2D& quad_rule,
                        const GenericKernel<T>& kernel,
                        const FaceInfo& face,
                        const Vec3<T>& obs_loc,
                        const Vec3<double>& obs_n) {

    Vec3<T> result = {0,0,0};
    for (unsigned int src_q = 0; src_q < quad_rule.x_hat.size(); src_q++) {
        SrcPointInfo<T> pt(quad_rule, kernel, face, obs_loc, obs_n, src_q);
        result += pt.weighted_kernel * pt.basis;
    }
    return result;
}

template <typename T>
T integral(const QuadratureRule2D& quad_rule,
           const GenericKernel<T>& kernel,
           const FaceInfo& face,
           const Vec3<double>& src_vals,
           const Vec3<T>& obs_loc,
           const Vec3<double>& obs_n) {

    T result = 0.0;
    for (unsigned int src_q = 0; src_q < quad_rule.x_hat.size(); src_q++) {
        SrcPointInfo<T> pt(quad_rule, kernel, face, obs_loc, obs_n, src_q);
        result += pt.weighted_kernel * dot(pt.basis, src_vals);
    }
    return result;
}

std::vector<double> integral_equation_vector(const Mesh& src_mesh,
                              const QuadratureRule2D& src_quad,
                              const Kernel& kernel,
                              const TaylorKernel& t_kernel,
                              const NearEval& near_eval, 
                              const Vec3<double>& obs_pt,
                              const Vec3<double>& obs_normal,
                              double obs_len_scale,
                              const double far_threshold = 3.0);

double eval_integral_equation(const Mesh& src_mesh,
                              const QuadratureRule2D& src_quad,
                              const Kernel& kernel,
                              const TaylorKernel& t_kernel,
                              const NearEval& near_eval, 
                              const Vec3<double>& obs_pt,
                              const Vec3<double>& obs_normal,
                              double obs_len_scale,
                              const std::vector<double>& src_strength,
                              const double far_threshold = 3.0);

std::vector<std::vector<double>> interact_matrix(const Mesh& src_mesh,
                                    const Mesh& obs_mesh,
                                    const QuadratureRule2D& src_quad,
                                    const QuadratureRule2D& obs_quad,
                                    const Kernel& kernel,
                                    const TaylorKernel& t_kernel,
                                    int n_steps,
                                    double far_threshold = 3.0);

std::vector<double> direct_interact(const Mesh& src_mesh,
                                    const Mesh& obs_mesh,
                                    const QuadratureRule2D& src_quad,
                                    const QuadratureRule2D& obs_quad,
                                    const Kernel& kernel,
                                    const TaylorKernel& t_kernel,
                                    const std::vector<double>& src_strength,
                                    int n_steps, 
                                    const double far_threshold = 3.0);

std::vector<double> mass_term(const Mesh& obs_mesh,
                              const QuadratureRule2D& obs_quad,
                              const std::vector<double>& strengths);


double get_len_scale(Mesh& mesh, int which_face, int q);
#endif

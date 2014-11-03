#ifndef __AAAAAAAAA_BEM_H
#define __AAAAAAAAA_BEM_H

#include <functional>
#include <array>
#include <vector>
#include "vec.h"

class QuadratureRule2D;
class Mesh;

/* Temporary data store for the evaluation of nearfield integrals. */
class NearEval {
public:
    NearEval(int n_steps);

    static constexpr double initial_dist = 1.0;

    const int n_steps;
    std::vector<QuadratureRule2D> quad;
    std::vector<double> dist;
};

typedef std::function<double (double,
                              const Vec3<double>&,
                              const Vec3<double>&,
                              const Vec3<double>&)> KernelFnc;

class FaceInfo {
public:
    FaceInfo(const Mesh& mesh, int face_index);
    
    const std::array<int,3>& face;
    double area;
    double jacobian;
    Vec3<double> normal;
    const std::array<Vec3<double>,3> corners;
};

/* Perform the richardson extrapolation for the nearfield quadrature. 
 */
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
        //TODO: Consider steps of two in error_order as an optional feature
        //TODO: Consider allowing setting the maximum error and then 
        // adaptively building
        //TODO: Integrate this into the main loop?
        //TODO: Use the diligenti mapping quadrature for the nearly 
        //singular quadratures required.
        last_level = this_level;
    }
    // std::cout << this_level[0] << std::endl;
    return this_level[0];
}

Vec3<double> basis_integrals(const QuadratureRule2D& quad_rule,
                                     const KernelFnc& kernel,
                                     const FaceInfo& face,
                                     const Vec3<double>& obs_loc,
                                     const Vec3<double>& obs_n);

double integral(const QuadratureRule2D& quad_rule,
                const KernelFnc& kernel,
                const FaceInfo& face,
                const Vec3<double>& src_vals,
                const Vec3<double>& obs_loc,
                const Vec3<double>& obs_n);

std::vector<double> integral_equation_vector(const Mesh& src_mesh,
                              const QuadratureRule2D& src_quad,
                              const KernelFnc& kernel,
                              const NearEval& near_eval, 
                              const Vec3<double>& obs_pt,
                              const Vec3<double>& obs_normal,
                              double obs_len_scale,
                              const double far_threshold = 3.0);

double eval_integral_equation(const Mesh& src_mesh,
                              const QuadratureRule2D& src_quad,
                              const KernelFnc& kernel,
                              const NearEval& near_eval, 
                              const Vec3<double>& obs_pt,
                              const Vec3<double>& obs_normal,
                              double obs_len_scale,
                              const std::vector<double>& src_strength,
                              const double far_threshold = 3.0);

std::vector<std::vector<double>> interact_matrix(Mesh& src_mesh,
                                    Mesh& obs_mesh,
                                    QuadratureRule2D& src_quad,
                                    QuadratureRule2D& obs_quad,
                                    const KernelFnc& kernel,
                                    int n_steps,
                                    double far_threshold = 3.0);

std::vector<double> direct_interact(Mesh& src_mesh,
                                    Mesh& obs_mesh,
                                    QuadratureRule2D& src_quad,
                                    QuadratureRule2D& obs_quad,
                                    const KernelFnc& kernel,
                                    std::vector<double>& src_strength,
                                    int n_steps, 
                                    const double far_threshold = 3.0);

std::vector<double> mass_term(const Mesh& obs_mesh,
                              const QuadratureRule2D& obs_quad,
                              const std::vector<double>& strengths);


double get_len_scale(Mesh& mesh, int which_face, int q);
#endif

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
                              Vec3<double>,
                              Vec3<double>)> KernelFnc;

double integral(const QuadratureRule2D& quad_rule,
                const KernelFnc& kernel,
                const std::array<Vec3<double>,3>& src_locs,
                const Vec3<double>& src_vals,
                const Vec3<double>& obs_loc);

double eval_integral_equation(const Mesh& src_mesh,
                              const QuadratureRule2D& src_quad,
                              const KernelFnc& kernel,
                              const NearEval& near_eval, 
                              std::array<double,3> obs_pt,
                              std::array<double,3> obs_normal,
                              double obs_len_scale,
                              std::vector<double>& src_strength,
                              const double far_threshold = 3.0);

double richardson_step(const std::vector<double>& values);

//TODO: TEst me!
std::vector<double> direct_interact(Mesh& src_mesh,
                                    Mesh& obs_mesh,
                                    QuadratureRule2D& src_quad,
                                    QuadratureRule2D& obs_quad,
                                    const KernelFnc& kernel,
                                    std::vector<double>& src_strength,
                                    int n_steps, 
                                    const double far_threshold = 3.0);

//TODO: TEst me!
std::vector<double> mass_term(const Mesh& obs_mesh,
                              const QuadratureRule2D& obs_quad,
                              const std::vector<double>& strengths);

double get_len_scale(Mesh& mesh, int which_face, int q);
#endif

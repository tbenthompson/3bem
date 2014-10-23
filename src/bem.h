#ifndef __AAAAAAAAA_BEM_H
#define __AAAAAAAAA_BEM_H

#include <functional>
#include <array>
#include <vector>

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
                              std::array<double,3>,
                              std::array<double,3>)> KernelFnc;

double integral(const QuadratureRule2D& quad_rule,
                const KernelFnc& kernel,
                const std::array<std::array<double,3>,3>& src_locs,
                const std::array<double,3>& src_vals,
                const std::array<double,3>& obs_loc);

double eval_integral_equation(const Mesh& src_mesh,
                              const QuadratureRule2D& src_quad,
                              const KernelFnc& kernel,
                              const NearEval& near_eval, 
                              std::array<double,3> obs_pt,
                              std::array<double,3> obs_normal,
                              std::vector<double>& src_strength);

double richardson_step(const std::vector<double>& values);
#endif

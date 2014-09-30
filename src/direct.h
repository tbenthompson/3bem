#ifndef __DIRECT_EVALUATION_H
#define __DIRECT_EVALUATION_H
#include <array>
#include <vector>
#include <functional>
#include "numerics.h"

typedef std::function<double (std::array<double,3>, std::array<double,3>)> Kernel;

inline double one_kernel(std::array<double,3>, std::array<double,3>) {
    return 1.0;
}

const double eps = 1e-15;
inline double laplace_single(std::array<double,3> x0, std::array<double,3> x1) {
    double r2 = dist2<3>(x0, x1);
    return 1.0 / (4 * M_PI * sqrt(r2 + eps));
}

inline std::vector<double> direct_n_body(std::vector<std::array<double,3>>& src_locs,
                                  std::vector<std::array<double,3>>& obs_locs,
                                  Kernel kernel,
                                  std::vector<double>& values) 
{
    std::vector<double> out_vals(obs_locs.size(), 0.0);
#pragma omp parallel for
    for (unsigned int i = 0; i < obs_locs.size(); ++i) {
        for (unsigned int j = 0; j < src_locs.size(); ++j) {
            out_vals[i] += kernel(obs_locs[i], src_locs[j]);
        }
    }
    return out_vals;
}
#endif

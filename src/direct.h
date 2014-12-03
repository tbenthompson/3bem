#ifndef __DIRECT_EVALUATION_H
#define __DIRECT_EVALUATION_H
#include <array>
#include <vector>
#include <functional>
#include "numerics.h"
#include <immintrin.h>

typedef std::function<double (double, double, double, double, double, double)> Kernel;

inline double one_kernel(double,double,double,double,double,double) {
    return 1.0;
}

const double eps = 1e-15;
inline double laplace_single(double x0, double y0, double z0,
                             double x1, double y1, double z1) {
    double dx = x0 - x1;
    double dy = y0 - y1;
    double dz = z0 - z1;
    double r2 = dx * dx + dy * dy + dz * dz;
    if (r2 < eps) {
        return 0.0;
    }
    return 1.0 / (4 * M_PI * std::sqrt(r2));
}

//Naive double for loop version of direct_n_body
inline std::vector<double> direct_n_body(std::array<std::vector<double>,3>& src_locs,
                                         std::array<std::vector<double>,3>& obs_locs,
                                         Kernel kernel,
                                         std::vector<double>& values) 
{
    std::vector<double> out_vals(obs_locs[0].size(), 0.0);
#pragma omp parallel for 
    for (unsigned int i = 0; i < obs_locs[0].size(); ++i) {
        for (unsigned int j = 0; j < src_locs[0].size(); ++j) {
            double kernel_val = kernel(obs_locs[0][i], obs_locs[1][i], obs_locs[2][i],
                                       src_locs[0][j], src_locs[1][j], src_locs[2][j]);
            out_vals[i] += values[j] * kernel_val;
        }
    }
    return out_vals;
}
#endif

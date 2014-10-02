#ifndef __DIRECT_EVALUATION_H
#define __DIRECT_EVALUATION_H
#include <array>
#include <vector>
#include <functional>
#include "numerics.h"
#include <immintrin.h>

typedef std::function<double (std::array<double,3>, std::array<double,3>)> Kernel;

inline double one_kernel(std::array<double,3>, std::array<double,3>) {
    return 1.0;
}

const double eps = 1e-15;
inline double laplace_single(std::array<double,3> x0, std::array<double,3> x1) {
    double r2 = dist2<3>(x0, x1);
    if (r2 < eps) {
        return 0.0;
    }
    return 1.0 / (4 * M_PI * sqrt(r2));
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
            double kernel_val = kernel(obs_locs[i], src_locs[j]);
            out_vals[i] += values[j] * kernel_val;
        }
    }
    return out_vals;
}

//TODO: Try this with aligned float arrays.
inline std::vector<float> vec_direct_n_body(std::array<std::vector<float>,3>& src_locs,
                                  std::array<std::vector<float>,3>& obs_locs,
                                  std::vector<float>& values) 
{
    const float factor = 1.0 / (4 * M_PI);
    const __m256 factor_rep = _mm256_broadcast_ss(&factor);

    std::vector<float> out_vals(obs_locs[0].size(), 0.0);
#pragma omp parallel for
    for (unsigned int i = 0; i < obs_locs[0].size(); i += 8) {
        __m256 temp_out = _mm256_setzero_ps();
        __m256 obs_loc_x = _mm256_loadu_ps(&obs_locs[0][i]);
        __m256 obs_loc_y = _mm256_loadu_ps(&obs_locs[1][i]);
        __m256 obs_loc_z = _mm256_loadu_ps(&obs_locs[2][i]);
        for (unsigned int j = 0; j < src_locs[0].size(); j++) {
            __m256 src_loc_x = _mm256_broadcast_ss(&src_locs[0][j]);
            __m256 src_loc_y = _mm256_broadcast_ss(&src_locs[1][j]);
            __m256 src_loc_z = _mm256_broadcast_ss(&src_locs[2][j]);
            __m256 src_values = _mm256_broadcast_ss(&values[j]);
            __m256 dx = _mm256_sub_ps(obs_loc_x, src_loc_x);
            __m256 dy = _mm256_sub_ps(obs_loc_y, src_loc_y);
            __m256 dz = _mm256_sub_ps(obs_loc_z, src_loc_z);
            __m256 dx2 = _mm256_mul_ps(dx, dx);
            __m256 dy2 = _mm256_mul_ps(dy, dy);
            __m256 dz2 = _mm256_mul_ps(dz, dz);
            __m256 r2 = _mm256_add_ps(_mm256_add_ps(dx2, dy2), dz2);
            __m256 inv_r = _mm256_rsqrt_ps(r2);
            __m256 kernel_eval = _mm256_mul_ps(inv_r, factor_rep);
            __m256 o = _mm256_mul_ps(kernel_eval, src_values);
            temp_out = _mm256_add_ps(o, temp_out);
        }
        _mm256_storeu_ps(&out_vals[i], temp_out);
    }
    return out_vals;
}
#endif

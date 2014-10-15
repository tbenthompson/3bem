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
    return 1.0 / (4 * M_PI * sqrt(r2));
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

//This version uses OpenMP/AVX on std::vectors.
inline 
std::vector<float> vec_direct_n_body(const std::array<std::vector<float>,3>& src_locs,
                                     const std::array<std::vector<float>,3>& obs_locs,
                                     unsigned int src_start, unsigned int src_end,
                                     unsigned int obs_start, unsigned int obs_end,
                                     const std::vector<float>& values) 
{
    const float factor = 1.0 / (4 * M_PI);
    const __m256 factor_rep = _mm256_broadcast_ss(&factor);

    std::vector<float> out_vals(obs_end - obs_start);
    for (unsigned int i = obs_start; i < obs_end - 8; i += 8) {
        __m256 temp_out = _mm256_setzero_ps();
        __m256 obs_loc_x = _mm256_loadu_ps(&obs_locs[0][i]);
        __m256 obs_loc_y = _mm256_loadu_ps(&obs_locs[1][i]);
        __m256 obs_loc_z = _mm256_loadu_ps(&obs_locs[2][i]);
        for (unsigned int j = src_start; j < src_end; j++) {
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
        _mm256_storeu_ps(&out_vals[i - obs_start], temp_out);
    }
    return out_vals;
}

inline 
std::vector<float> vec_direct_n_body(const std::array<std::vector<float>,3>& src_locs,
                                     const std::array<std::vector<float>,3>& obs_locs,
                                     std::vector<float>& values) {
    return vec_direct_n_body(src_locs, obs_locs, 
                             0, src_locs.size(),
                             0, obs_locs.size(), values);
}


//This version uses OpenMP/AVX on aligned float arrays. It is not or is barely faster
//than the unaligned version. Experimentation!
// Further: does blocking help? Maybe not.
// Further: fused multiply add -> a bit faster (~15 GFlop/s)
// Further: loop unrolling? --> YES! (210 GFlop/s -> 255 GFlop/s
// Best result so far: ~ 89 instr/cycle --> 302 Gigaflops
inline float* really_fast_vec_direct_n_body(float* srcfx, float* srcfy, float* srcfz,
                                            float* obsfx, float* obsfy, float* obsfz,
                                            float* values, int n_src, int n_obs)
{
    const float factor = 1.0 / (4 * M_PI);
    const __m256 factor_rep = _mm256_broadcast_ss(&factor);

    float* out_vals = (float*)_mm_malloc(n_obs * sizeof(float), 32);
    for (int i = 0; i < n_obs; i++) {
        out_vals[i] = 0.0;
    }
    int block_width = 800;
    int obs_blocks = n_obs / block_width;
    int src_blocks = n_src / block_width;
    int blocks = obs_blocks * src_blocks;
#pragma omp parallel for
    for (int b = 0; b < blocks; b++) {
        int i_block = b % obs_blocks;
        int i_start = block_width * i_block;
        int i_end = block_width * (i_block + 1);

        int j_block = (b - i_block) / obs_blocks;
        int j_start = block_width * j_block;
        int j_end = block_width * (j_block + 1);
        for (int i = i_start; i < i_end; i += 8) {
            auto obs_loc_x = _mm256_load_ps(&obsfx[i]);
            auto obs_loc_y = _mm256_load_ps(&obsfy[i]);
            auto obs_loc_z = _mm256_load_ps(&obsfz[i]);
            __m256 temp0 = _mm256_setzero_ps();
            __m256 temp1 = _mm256_setzero_ps();
            __m256 temp2 = _mm256_setzero_ps();
            __m256 temp3 = _mm256_setzero_ps();
            for (int j = j_start; j < j_end; j += 4) {
                auto src_loc_x0 = _mm256_broadcast_ss(&srcfx[j]);
                auto src_loc_x1 = _mm256_broadcast_ss(&srcfx[j + 1]);
                auto src_loc_x2 = _mm256_broadcast_ss(&srcfx[j + 2]);
                auto src_loc_x3 = _mm256_broadcast_ss(&srcfx[j + 3]);
                auto src_loc_y0 = _mm256_broadcast_ss(&srcfy[j]);
                auto src_loc_y1 = _mm256_broadcast_ss(&srcfy[j + 1]);
                auto src_loc_y2 = _mm256_broadcast_ss(&srcfy[j + 2]);
                auto src_loc_y3 = _mm256_broadcast_ss(&srcfy[j + 3]);
                auto src_loc_z0 = _mm256_broadcast_ss(&srcfz[j]);
                auto src_loc_z1 = _mm256_broadcast_ss(&srcfz[j + 1]);
                auto src_loc_z2 = _mm256_broadcast_ss(&srcfz[j + 2]);
                auto src_loc_z3 = _mm256_broadcast_ss(&srcfz[j + 3]);
                auto src_values0 = _mm256_broadcast_ss(&values[j]);
                auto src_values1 = _mm256_broadcast_ss(&values[j + 1]);
                auto src_values2 = _mm256_broadcast_ss(&values[j + 2]);
                auto src_values3 = _mm256_broadcast_ss(&values[j + 3]);
                auto dx0 = obs_loc_x - src_loc_x0;
                auto dx1 = obs_loc_x - src_loc_x1;
                auto dx2 = obs_loc_x - src_loc_x2;
                auto dx3 = obs_loc_x - src_loc_x3;
                auto dy0 = obs_loc_y - src_loc_y0;
                auto dy1 = obs_loc_y - src_loc_y1;
                auto dy2 = obs_loc_y - src_loc_y2;
                auto dy3 = obs_loc_y - src_loc_y3;
                auto dz0 = obs_loc_z - src_loc_z0;
                auto dz1 = obs_loc_z - src_loc_z1;
                auto dz2 = obs_loc_z - src_loc_z2;
                auto dz3 = obs_loc_z - src_loc_z3;
                auto r20 = dx0 * dx0;
                auto r21 = dx1 * dx1;
                auto r22 = dx2 * dx2;
                auto r23 = dx3 * dx3;
                r20 = dy0 * dy0 + r20;
                r21 = dy1 * dy1 + r21;
                r22 = dy2 * dy2 + r22;
                r23 = dy3 * dy3 + r23;
                r20 = dz0 * dz0 + r20;
                r21 = dz1 * dz1 + r21;
                r22 = dz2 * dz2 + r22;
                r23 = dz3 * dz3 + r23;
                auto src_times_factor0 = _mm256_mul_ps(src_values0, factor_rep);
                auto src_times_factor1 = _mm256_mul_ps(src_values1, factor_rep);
                auto src_times_factor2 = _mm256_mul_ps(src_values2, factor_rep);
                auto src_times_factor3 = _mm256_mul_ps(src_values3, factor_rep);
                auto inv_r0 = _mm256_rsqrt_ps(r20);
                auto inv_r1 = _mm256_rsqrt_ps(r21);
                auto inv_r2 = _mm256_rsqrt_ps(r22);
                auto inv_r3 = _mm256_rsqrt_ps(r23);
                temp0 = inv_r0 * src_times_factor0 + temp0;
                temp1 = inv_r1 * src_times_factor1 + temp1;
                temp2 = inv_r2 * src_times_factor2 + temp2;
                temp3 = inv_r3 * src_times_factor3 + temp3;
            }
            auto cur_value = _mm256_load_ps(&out_vals[i]);
            cur_value = cur_value + temp0 + temp1 + temp2 + temp3;
            _mm256_store_ps(&out_vals[i], cur_value);
        }
    }
    return out_vals;
}
#endif

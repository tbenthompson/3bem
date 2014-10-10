#include "UnitTest++.h"
#include <iostream>
#include <math.h>
#include "test_shared.h"
#include "direct.h"

int main() {
    int n = 12 * 800 * 10;
    bool check = n < 2000;
    std::array<std::vector<double>,3> src =
        {random_list(n), random_list(n), random_list(n)};
    std::array<std::vector<double>,3> obs =
        {random_list(n), random_list(n), random_list(n)};
    float *obsfx = (float*)_mm_malloc(n * sizeof(float), 32);
    float *obsfy = (float*)_mm_malloc(n * sizeof(float), 32);
    float *obsfz = (float*)_mm_malloc(n * sizeof(float), 32);
    float *srcfx = new float[n];
    float *srcfy = new float[n];
    float *srcfz = new float[n];
    float *fast_values = new float[n];
    std::vector<float> values(n, 1.0f);
    std::array<std::vector<float>,3> srcf;
    std::array<std::vector<float>,3> obsf;
    std::array<std::vector<double>,3> slow_src;
    std::array<std::vector<double>,3> slow_obs;
    std::vector<double> slow_values(n, 1.0);
    for (int d = 0; d < 3; d++) {
        srcf[d].resize(n);
        obsf[d].resize(n);
        slow_src[d].resize(n);
        slow_obs[d].resize(n);
        for (int i = 0; i < n; i++) {
            srcfx[i] = src[0][i];
            obsfx[i] = obs[0][i];
            srcfy[i] = src[1][i];
            obsfy[i] = obs[1][i];
            srcfz[i] = src[2][i];
            obsfz[i] = obs[2][i];
            slow_src[i][d] = src[d][i];
            slow_obs[i][d] = obs[d][i];
            fast_values[i] = 1.0f;
            srcf[d][i] = src[d][i];
            obsf[d][i] = obs[d][i];
        }
    }
    TIC;
    auto result1 = vec_direct_n_body(srcf, obsf, values);
    TOC("Fast Direct N Body");
    TIC2;
    auto result2 = really_fast_vec_direct_n_body(srcfx, srcfy, srcfz, obsfx, obsfy, obsfz, fast_values, n, n);
    TOC("Not Really Much Faster Direct N Body");
    long interacts = ((long)n) * ((long)n);
    // In order: 3 Subtracts, 1 multiply, 2 FMAs (count 2), 1 multiply, 1 inv_sqrt (count for 4), 1 FMA (count 2)
    double ops_per_interact = 15;
    double proc = 3.4e9;
    double cycles = (((double)time_ms) / 1000.0) * proc;
    double instr_per_cycle = (interacts * ops_per_interact) / cycles;
    double gigaflops = instr_per_cycle * proc / 1e9;
    // Conclusion: aligned memory access is not actually faster.
    // approximate instructions per cycle = 72 --> very close to optimal
    std::cout << "Cycles: " << cycles << std::endl;
    std::cout << "Instructions/cycle: " << instr_per_cycle << std::endl;
    std::cout << "GFlop/s: " << gigaflops << std::endl;

    if (check) {
        Kernel k = laplace_single;
        auto exact = direct_n_body(slow_src, slow_obs, k, slow_values);
        // CHECK_ARRAY_CLOSE(exact, result1, n, 1e-1);
        CHECK_ARRAY_CLOSE(exact, result2, n, 1e-1);
    }
    // for (auto r: result) {
    //     std::cout << r << std::endl;
    // }
}

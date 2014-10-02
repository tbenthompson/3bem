#include "UnitTest++.h"
#include <iostream>
#include <math.h>
#include "test_shared.h"
#include "direct.h"

const double PI = 4.0 * std::atan(1.0);

TEST(Direct) {
    CHECK_CLOSE(one_kernel({0, 0, 0}, {0, 0, 0}), 1.0, 1e-14);

    int n = 8 * 250;
    auto src = random_pts(n);
    auto obs = random_pts(n);
    std::vector<double> values(n, 1.0);

    Kernel k = laplace_single;
    TIC
    auto result = direct_n_body(src, obs, k, values);
    TOC("Direct N Body");
    // std::cout << "number of interaction: " << n * n << std::endl;
    


    // std::vector<double> correct(n);
    // for (int i = 0; i < n; ++i) correct[i] = n;
    // CHECK_ARRAY_CLOSE(correct, result, n - 1, 1e-8);
}

TEST(FastDirect) {
    int n = 8 * 12500;
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
    std::vector<std::array<double, 3>> slow_src(n);
    std::vector<std::array<double, 3>> slow_obs(n);
    std::vector<double> slow_values(n, 1.0);
    for (int d = 0; d < 3; d++) {
        srcf[d].resize(n);
        obsf[d].resize(n);
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
    TOC("Not Really Faster Direct N Body");
    // Conclusion: aligned memory access is not actually faster.
    // std::cout << "number of interaction: " << ((long)n) * ((long)n) << std::endl;

    if (check) {
        Kernel k = laplace_single;
        auto exact = direct_n_body(slow_src, slow_obs, k, slow_values);
        CHECK_ARRAY_CLOSE(exact, result1, n, 1e-1);
        CHECK_ARRAY_CLOSE(exact, result2, n, 1e-1);
    }
    // for (auto r: result) {
    //     std::cout << r << std::endl;
    // }
}

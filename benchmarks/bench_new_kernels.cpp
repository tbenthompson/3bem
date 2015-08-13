#include "benchmark/benchmark.h"
#include "new_laplace_kernels.h"
#include "laplace_kernels.h"
#include "util.h"

using namespace tbem;

static void old_laplace_single_kernel(benchmark::State& state)
{
    size_t n = state.range_x();
    NEWLaplaceSingle<2> K;
    auto obs_pts = random_pts<2>(n);
    auto src_pts = random_pts<2>(n);
    auto obs_normals = random_pts<2>(n);
    auto src_normals = random_pts<2>(n);
    while (state.KeepRunning()) {
        std::vector<double> data(n * n);
        for (size_t i = 0; i < n; i++) {
            for (size_t j = 0; j < n; j++) {
                auto correct = LaplaceSingle<2>()(
                    obs_pts[i], src_pts[j],
                    obs_normals[i], src_normals[j]
                );
                data[i * n + j] = correct[0][0];
            }
        }
    }
}

static void laplace_single_kernel(benchmark::State& state)
{
    size_t n = state.range_x();
    NEWLaplaceSingle<2> K;
    auto obs_pts = random_pts<2>(n);
    auto src_pts = random_pts<2>(n);
    auto obs_normals = random_pts<2>(n);
    auto src_normals = random_pts<2>(n);
    while (state.KeepRunning()) {
        auto result = K(obs_pts, src_pts, obs_normals, src_normals);
    }
}

BENCHMARK(old_laplace_single_kernel)->Range(1, 2000);
BENCHMARK(laplace_single_kernel)->Range(1, 2000);

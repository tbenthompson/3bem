#include "catch.hpp"
#include "new_laplace_kernels.h"
#include "laplace_kernels.h"
#include "util.h"

using namespace tbem;

TEST_CASE("Laplace", "[new_kernels]")
{
    size_t n = 3;
    NEWLaplaceSingle<2> K;
    auto obs_pts = random_pts<2>(n);
    auto src_pts = random_pts<2>(n);
    auto obs_normals = random_pts<2>(n);
    auto src_normals = random_pts<2>(n);
    auto result = K(obs_pts, src_pts, obs_normals, src_normals);
    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < n; j++) {
            auto correct = LaplaceSingle<2>()(
                obs_pts[i], src_pts[j],
                obs_normals[i], src_normals[j]
            );
            REQUIRE_CLOSE(correct[0][0], result[i * n + j], 1e-10);
        }
    }
}

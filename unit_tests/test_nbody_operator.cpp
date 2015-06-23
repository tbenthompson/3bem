#include "catch.hpp"
#include "nbody_operator.h"
#include "laplace_kernels.h"
#include "util.h"
#include "mesh_gen.h"
#include "gauss_quad.h"

using namespace tbem;

TEST_CASE("nbody from bem", "[nbody_operator]")
{
    auto m = circle_mesh({0, 0}, 1.0, 1);
    auto q_src = gauss(2);
    auto q_obs = gauss(2);
    auto data = nbody_data_from_bem(m, m, q_obs, q_src);
    REQUIRE(data.obs_locs.size() == 16);
    REQUIRE_ARRAY_CLOSE(
        data.obs_locs[0], ref_to_real(q_obs[0].x_hat, m.facets[0]), 2, 1e-10
    );
    REQUIRE_ARRAY_CLOSE(
        data.src_locs[0], ref_to_real(q_src[0].x_hat, m.facets[0]), 2, 1e-10
    );
    auto edge_length = dist(m.facets[0][0], m.facets[0][1]);
    REQUIRE_CLOSE(data.src_weights[0], 0.5 * edge_length * q_src[0].w, 1e-10);
}

TEST_CASE("simple nbody operator", "[nbody_operator]")
{
    LaplaceSingle<3> K;
    NBodyData<3> data{
        {{0, 0, 0}, {0, 2, 0}, {0, 0, -1}},
        {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}},
        {{1, 0, 0}, {0, 1, 0}},
        {{0, 0, 0}, {0, 0, 0}},
        {1.0}
    };
    auto op = make_direct_nbody_operator(data, K);
    auto out = op.apply({4.0, 0.0});
    REQUIRE(out.size() == 3);
    REQUIRE_ARRAY_CLOSE(out, std::vector<double>{0.3183, 0.1424, 0.2251}, 3, 1e-3);
}

TEST_CASE("compare with eval", "[nbody_operator]") 
{
    size_t n = 20;
    NBodyData<2> data{
        random_pts<2>(n), random_pts<2>(n), 
        random_pts<2>(n), random_pts<2>(n), random_list(n)
    };
    auto input = random_list(n);
    LaplaceHypersingular<2> K;
    auto from_op = make_direct_nbody_operator(data, K).apply(input);
    auto eval = nbody_eval(K, data, input.data());
    REQUIRE_ARRAY_CLOSE(from_op, eval, n, 1e-12);
}

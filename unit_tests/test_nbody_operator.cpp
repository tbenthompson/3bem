#include "catch.hpp"
#include "nbody_operator.h"
#include "laplace_kernels.h"
#include "util.h"

using namespace tbem;

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
    // auto op = BlockDirectNBodyOperator<3,1,1>{data, K};
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

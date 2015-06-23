#include "catch.hpp"
#include "galerkin_operator.h"
#include "gauss_quad.h"
#include "mesh.h"
#include "sparse_operator.h"

using namespace tbem;

TEST_CASE("Galerkin one facet", "[galerkin_operator]") 
{
    Mesh<3> m{{
        {{
            {{0.0, 0.0, 0.0}},
            {{1.0, 0.0, 0.0}},
            {{0.0, 1.0, 0.0}}
        }}
    }};

    auto g = make_galerkin_operator(1, m, tri_gauss(3));

    SECTION("Shape") {
        REQUIRE(g.n_cols() == 9);
        REQUIRE(g.n_rows() == 3);
    }

    SECTION("Apply") {
        auto result = g.apply(std::vector<double>(g.n_cols(), 1.0));
        REQUIRE(result.size() == 3);
        // The result should be equal to:
        // \int_0^1 \int_0^{1-x} (1-x-y) dy dx = 1 / 6
        // Although that is only the integral of the first basis function, all
        // the other basis functions are symmetric
        REQUIRE_ARRAY_CLOSE(result, std::vector<double>(3, 1.0 / 6.0), 3, 1e-10);
    }
}

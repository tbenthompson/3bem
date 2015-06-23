#include "catch.hpp"
#include "interpolation_operator.h"
#include "gauss_quad.h"
#include "mesh.h"
#include "sparse_operator.h"

using namespace tbem;

TEST_CASE("Interpolation one facet", "[interpolation_operator]")
{
    Mesh<3> m{{
        {{
            {{0.0, 0.0, 0.0}},
            {{1.0, 0.0, 0.0}},
            {{0.0, 1.0, 0.0}}
        }}
    }};

    auto interp = make_interpolation_operator(1, m, tri_gauss(3));

    SECTION("Shape") {
        REQUIRE(interp.n_cols() == 3);
        REQUIRE(interp.n_rows() == 9);
    }

    SECTION("Apply") {
        auto result = interp.apply(std::vector<double>(3, 1.0));
        REQUIRE(result.size() == 9);
        REQUIRE_ARRAY_CLOSE(result, std::vector<double>(9, 1.0), 9, 1e-10);
    }
}

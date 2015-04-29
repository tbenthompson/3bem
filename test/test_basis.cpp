#include "catch.hpp"
#include "basis.h"
#include "mesh_gen.h"
#include "mesh.h"
#include "vectorx.h"

using namespace tbem;

TEST_CASE("Interpolate", "[basis]") {
    auto m = sphere_mesh({0,0,0}, 1, 2);
    auto res = interpolate<3>(m, [](const Vec<double,3>& x) {return x[0];});

    for (unsigned int i = 0; i < m.n_facets(); i++) {
        for (int d = 0; d < 3; d++) {
            REQUIRE(res[3 * i + d] == Approx(m.facets[i][d][0]).epsilon(1e-14));
        }
    }
}

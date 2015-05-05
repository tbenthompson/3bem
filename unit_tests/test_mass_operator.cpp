#include "catch.hpp"
#include "mesh_gen.h"
#include "util.h"
#include "identity_kernels.h"
#include "dense_builder.h"
#include "mass_operator.h"
#include "laplace_kernels.h"

using namespace tbem;

TEST_CASE("MassTerm", "[mass_operator]") {
    auto sphere = sphere_mesh({0,0,0}, 1.0, 0);
    std::vector<double> str(sphere.n_dofs(), 1.0);
    for (std::size_t i = 0; i < sphere.n_facets(); i++) {
        for (int d = 0; d < 3; d++) {
            if (sphere.facets[i][d][0] > 0.5) {
                str[3 * i + d] = 0.0;
            }
        }
    }
    auto mass_op = mass_operator<3,1,1>(sphere, 2);
    auto res = mass_op.apply(str);
    REQUIRE(res.size() == sphere.n_dofs());
    double true_area = 0.0;
    for (auto f: sphere.facets) {
        true_area += tri_area(f);
    }
    double mass_area = 0.0;
    for (auto r: res) {
        mass_area += r;
    }   
    REQUIRE_CLOSE(mass_area, (5.0 / 6.0) * true_area, 1e-12);
}

TEST_CASE("TensorMassTerm", "[mass_operator]") {
    auto sphere = sphere_mesh({0,0,0}, 1.0, 1);
    std::vector<double> str(3 * sphere.n_dofs(), 1.0);
    auto mass_op = mass_operator<3,3,3>(sphere, 2);
    auto res = mass_op.apply(str);
    REQUIRE((mass_op.n_block_rows() * mass_op.n_block_cols() == 9));
    REQUIRE(mass_op.n_total_rows() == 3 * sphere.n_dofs());
    REQUIRE(mass_op.galerkin.n_total_cols() == 3 * sphere.n_facets() * 4);
    REQUIRE(res.size() == mass_op.n_total_rows());

    double true_area = 0.0;
    for (auto f: sphere.facets) {
        true_area += tri_area(f);
    }

    double mass_area = 0.0;
    for (size_t i = 0; i < sphere.n_dofs(); i++) {
        mass_area += res[i];    
    }
    REQUIRE_CLOSE(mass_area, true_area, 1e-12);
}

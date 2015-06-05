#include "catch.hpp"
#include "nearfield_operator.h"
#include "continuity_builder.h"
#include "quadrature.h"
#include "laplace_kernels.h"
#include "mesh_gen.h"
#include "integral_operator.h"
#include "boost_geometry_wrapper.h"

using namespace tbem;

TEST_CASE("richardson points always inside", "[nearfield_operator]")
{
    std::vector<Vec<double,2>> polygon{
        {1, 0.1}, {0, 0}, {1, -0.1}
    };
    polygon.push_back(polygon[0]);

    auto R = 3;
    std::vector<Mesh<2>> mesh_pieces;
    for (size_t i = 0; i < polygon.size() - 1; i++) {
        mesh_pieces.push_back(line_mesh(polygon[i], polygon[i + 1]));
    }
    auto m = Mesh<2>::create_union(mesh_pieces).refine_repeatedly(R);

    auto pts = galerkin_obs_pts(m, gauss_facet<2>(2), m);

    for(auto p: pts) {
        REQUIRE(in_polygon(polygon, p.loc + p.richardson_dir));
    }
}

TEST_CASE("interior obs pts", "[nearfield_operator]")
{
    auto pts = interior_obs_pts<2>(
        {{0,1}, {1,0}},
        {{1,0}, {0,1}},
        Mesh<2>{{
            {{
                {0, 0}, {1, 0}
            }}
        }}
    );
    
    SECTION("size") {
        REQUIRE(pts.size() == 2);
    }

    SECTION("richardson length") {
        auto dir = pts[0].richardson_dir;
        REQUIRE(hypot(pts[0].richardson_dir) == 0.0);
        auto length = hypot(pts[1].richardson_dir);
        REQUIRE_CLOSE(length, std::sqrt(1.25) / 2, 1e-12);
    }
}

TEST_CASE("nearfield obs pts", "[nearfield_operator]")
{
    auto m = circle_mesh({0, 0}, 1.0, 1);
    REQUIRE(m.facets.size() == 8);
    auto pts = galerkin_obs_pts(m, gauss(2), m);
    REQUIRE(pts.size() == 16);
}

TEST_CASE("Constrained nearfield matrix", "[nearfield_operator]")
{
    auto m = circle_mesh({0, 0}, 1.0, 2);
    // Use a very large nearfield range (300000) to force everything to
    // be nearfield.
    QuadStrategy<2> qs(3, 8, 300000, 1e-4);
    LaplaceDouble<2> k;
    auto mthd = make_adaptive_integration_mthd(qs, k);
    auto galerkin = make_galerkin_operator(1, m, mthd.get_obs_quad());
    auto obs_pts = galerkin_obs_pts(m, mthd.get_obs_quad(), m);
    auto nearfield = make_nearfield_operator(obs_pts, m, mthd);
    auto matrix = galerkin.right_multiply(nearfield);
    auto dense_matrix = matrix.to_dense();

    SECTION("Dense operator equals galerkin nearfield") {
        auto dense_matrix2 = *dense_integral_operator(m, m, mthd, {m}).storage;
        REQUIRE_ARRAY_CLOSE(dense_matrix2, dense_matrix, dense_matrix2.size(), 1e-12);
    }

    SECTION("Dense condensation and sparse condensation are the same") {
        auto cm = from_constraints(
            convert_to_constraints(mesh_continuity(m.begin()))
        );
        auto test = condense_matrix(cm, cm, matrix).to_dense().data();
        auto correct = condense_matrix(cm, cm, dense_matrix).data();
        REQUIRE(test.size() == correct.size());
        REQUIRE_ARRAY_CLOSE(test, correct, correct.size(), 1e-12);
    }
}

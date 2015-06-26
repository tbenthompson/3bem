#include "catch.hpp"
#include "nearfield_operator.h"
#include "continuity_builder.h"
#include "gauss_quad.h"
#include "laplace_kernels.h"
#include "mesh_gen.h"
#include "integral_operator.h"
#include "gte_wrapper.h"
#include "util.h"

//TODO: remove
#include <iomanip>

using namespace tbem;

void test_polygon(const std::vector<Vec<double,2>>& polygon, size_t refine)
{
    // turn the points into a 3bem mesh
    std::vector<Mesh<2>> mesh_pieces;
    for (size_t i = 0; i < polygon.size() - 1; i++) {
        mesh_pieces.push_back(line_mesh(polygon[i], polygon[i + 1]));
    }
    auto m = Mesh<2>::create_union(mesh_pieces).refine_repeatedly(refine);

    // generate the observation points for a galerkin evaluation
    auto pts = galerkin_obs_pts(m, gauss_facet<2>(2), m);
    
    // and check that they are all with the polygon
    for(auto p: pts) {
        auto furthest_limit_pt = p.loc + p.richardson_dir;
        auto success = in_polygon(polygon, furthest_limit_pt);
        // std::cout << "polygon:" << std::endl;
        // for (size_t i = 0; i < polygon.size(); i++) {
        //     std::cout << polygon[i] << std::endl;
        // }
        // std::cout << "loc: " << std::setprecision(17) << p.loc << std::endl;
        // std::cout << "dir: " << p.richardson_dir << std::endl;
        // std::cout << "limitpt: " << furthest_limit_pt << std::endl;
        REQUIRE(success);
    }
}

TEST_CASE("richardson points always inside", "[nearfield_operator]")
{
    //TODO: This test should probably be in the test_limit_direction file
    for (size_t i = 0; i < 30; i++) {
        // make a random triangle
        auto polygon = random_pts<2>(3); 

        // then ensure that it is oriented counterclockwise by flipping the 
        // orientation if it isn't
        auto normal = unscaled_normal({
            Vec<double,3>{polygon[0][0], polygon[0][1], 0.0},
            Vec<double,3>{polygon[1][0], polygon[1][1], 0.0},
            Vec<double,3>{polygon[2][0], polygon[2][1], 0.0}
        });

        if (normal[2] < 0) {
            std::swap(polygon[0], polygon[1]);
        }

        // the polygon must be closed (note that this cannot be done before 
        // the clockwise check because vertices are swapped
        polygon.push_back(polygon[0]);

        test_polygon(polygon, 5);
    }
}

TEST_CASE("richardson points inside wedge -- regression test", "[nearfield_operator]")
{
    test_polygon({{0, 0}, {150e3, -25e3}, {200e3, 0}, {150e3, 5e3}, {0, 0}}, 6);
}

TEST_CASE("interior obs pts", "[nearfield_operator]")
{
    auto pts = interior_obs_pts<2>({{0, 1}, {1, 0}}, {{0, 1}, {0, 1}},
        Mesh<2>{{
            {{
                {0, 0}, {1, 0}
            }}
        }}
    );
    REQUIRE(pts.size() == 2);

    REQUIRE(hypot(pts[0].richardson_dir) == 0.0);
    auto correct = std::sqrt(1.25) * 0.4;
    REQUIRE_CLOSE(hypot(pts[1].richardson_dir), correct, 1e-12);
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
    LaplaceDouble<2> k;
    auto mthd = make_adaptive_integrator(1e-4, 3, 8, 300000, k);
    auto galerkin = make_galerkin_operator(1, m, mthd.obs_far_quad);
    auto obs_pts = galerkin_obs_pts(m, mthd.obs_far_quad, m);
    auto nearfield = make_nearfield_operator(obs_pts, m, mthd);
    auto matrix = galerkin.right_multiply(nearfield);
    auto dense_matrix = matrix.to_dense();

    SECTION("Dense operator equals galerkin nearfield") {
        auto dense_matrix2 = *dense_boundary_operator(m, m, mthd, {m}).storage;
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

#include "catch.hpp"
#include "nearfield_operator.h"
#include "continuity_builder.h"
#include "quadrature.h"
#include "laplace_kernels.h"
#include "mesh_gen.h"
#include "integral_operator.h"

#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point_xy.hpp>
#include <boost/geometry/geometries/polygon.hpp>

using namespace tbem;

TEST_CASE("richardson direction far", "[nearfield_operator]") 
{
    Facet<2> f{{{0, -1}, {1, -1}}};
    Vec<double,2> p{0, 1};

    auto dir = decide_richardson_dir(p, {{f}, f, p, 0.0});

    REQUIRE(dir == (Vec<double,2>{0, 1}));
}

TEST_CASE("richardson direction touching", "[nearfield_operator]") 
{
    Facet<2> f{{{0, 0}, {0, 1}}};
    Vec<double,2> p{0, 0};

    auto dir = decide_richardson_dir(p, {{f}, f, p, 0.0});

    REQUIRE(dir == (Vec<double,2>{-1, 0.5}));
}

TEST_CASE("richardson direction backside", "[nearfield_operator]") 
{
    Facet<2> f{{{0, -1}, {1, -1}}};
    Vec<double,2> p{0, -2};

    auto dir = decide_richardson_dir(p, {{f}, f, p, 0.0});

    std::cout << dir << std::endl;
    REQUIRE(dir == (Vec<double,2>{0, -1}));
}

TEST_CASE("richardson direction intersection", "[nearfield_operator]") 
{
    Facet<2> f{{{0, 0}, {1, 0}}};
    Facet<2> f2{{{0, 1}, {0, 0}}};
    Vec<double,2> p{0, 0};
    
    auto dir = decide_richardson_dir(p, {{f, f2}, f, p, 0.0});
    
    REQUIRE(dir == (Vec<double,2>{0.5, 0.5}));
}

TEST_CASE("richardson direction normals dont cancel", "[nearfield_operator]") 
{
    Facet<2> f{{{1, -1}, {1, 1}}};
    Facet<2> f2{{{-1, 1}, {-1, -1}}};
    Vec<double,2> p{0, 0};

    auto dir = decide_richardson_dir(p, {{f, f2}, f, p, 0.0});

    bool normals_not_canceled = (dir[0] != 0.0) || (dir[1] != 0.0);
    REQUIRE(normals_not_canceled);
}

TEST_CASE("richardson direction for reflex angle", "[nearfield_operator]")
{
    Facet<2> f{{{1, 0}, {0, 0}}};
    Facet<2> f2{{{0, 0}, {-1, 1}}};

    Vec<double,2> p{0.001, 0};
    auto dir = decide_richardson_dir(p, {{f, f2}, f, p, 0});
    REQUIRE_ARRAY_CLOSE(dir, Vec<double,2>{-0.5, -1.0}, 2, 1e-12);
    std::cout << dir << std::endl;
}

TEST_CASE("richardson points always inside", "[nearfield_operator]")
{
    std::vector<Vec<double,2>> polygon{
        {1, 0.1}, {0, 0}, {1, -0.1}
    };
    polygon.push_back(polygon[0]);

    auto R = 5;
    std::vector<Mesh<2>> mesh_pieces;
    for (size_t i = 0; i < polygon.size() - 1; i++) {
        mesh_pieces.push_back(line_mesh(polygon[i], polygon[i + 1]));
    }
    auto m = Mesh<2>::create_union(mesh_pieces).refine_repeatedly(R);

    auto pts = galerkin_obs_pts(m, gauss_facet<2>(2), m);
    // auto pts = interior_obs_pts({{0.149646, 0.0149646}}, {{0, 1}}, m);

    typedef boost::geometry::model::d2::point_xy<double> point_type;
    typedef boost::geometry::model::polygon<point_type> polygon_type;
    polygon_type poly;
    for (size_t i = 0; i < polygon.size(); i++) {
        poly.outer().push_back({polygon[i][0], polygon[i][1]});
    }

    for(auto p: pts) {
        auto interior_pt = p.loc + p.len_scale * p.richardson_dir;
        point_type bg_p(interior_pt[0], interior_pt[1]);
        REQUIRE(boost::geometry::within(bg_p, poly));
    }
}

TEST_CASE("interior obs pts", "[nearfield_operator]")
{
    auto pts = interior_obs_pts<2>(
        {{0,1}, {2,0}},
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

    SECTION("normalized richardson direction") {
        auto dir = pts[0].richardson_dir;
        REQUIRE(hypot(dir) == 1.0);
    }

    SECTION("richardson length divided by 5") {
        auto len = pts[0].len_scale;
        REQUIRE(len == 0.2);
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

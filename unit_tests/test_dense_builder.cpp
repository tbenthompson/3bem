#include "catch.hpp"
#include "dense_builder.h"
#include "integral_operator.h"
#include "numerics.h"
#include "quadrature.h"
#include "mesh.h"
#include "mesh_gen.h"
#include "util.h"
#include "identity_kernels.h"
#include "laplace_kernels.h"
#include "elastic_kernels.h"

//TODO: This file should be split up

using namespace tbem;
struct EvalProb {
    EvalProb(int refine_level, int near_eval, int gauss_order,
             Vec3<double> center = Vec3<double>{0,0,0},
             double r = 3.0):
        sphere(sphere_mesh(center, r, refine_level)),
        qs(gauss_order, near_eval, 2.0, 1e-2),
        obs_pt(random_pt<3>()),
        obs_n(random_pt<3>()),
        src_strength(std::vector<double>(sphere.n_dofs(), 1.0))
    {}

    double go(const Kernel<3,1,1>& k) 
    {
        auto mthd = make_adaptive_integration_mthd(qs, k);
        auto op = mesh_to_points_operator({obs_pt}, {obs_n}, sphere, mthd, {sphere});
        return op.apply(src_strength)[0];
    }

    Mesh<3> sphere;
    QuadStrategy<3> qs;
    Vec3<double> obs_pt;
    Vec3<double> obs_n;
    std::vector<double> src_strength;
}; 

TEST_CASE("EvalIntegralEquationSphereSurfaceArea", "[dense_builder]") 
{
    EvalProb ep(5, 3, 2);
    IdentityScalar<3> identity;
    double result = ep.go(identity);
    double exact_surf_area = 4*M_PI*9;
    REQUIRE_CLOSE(result, exact_surf_area, 1e-1);
}

TEST_CASE("ConstantLaplace", "[dense_builder]") 
{
    EvalProb ep(5, 3, 2);
    LaplaceDouble<3> double_kernel;
    double result = ep.go(double_kernel);
    REQUIRE_CLOSE(result, -1.0, 1e-3);
}

TEST_CASE("MatrixRowVsEval", "[dense_builder]") 
{
    EvalProb ep(4, 3, 2);
    double result = ep.go(LaplaceDouble<3>());
    REQUIRE_CLOSE(result, -1.0, 1e-3);
}

TEST_CASE("FacetInfo2D", "[dense_builder]") 
{
    Facet<2> f{Vec2<double>{0.0, 0.0}, Vec2<double>{3.0, 0.0}};
    auto face_info = FacetInfo<2>::build(f);
    REQUIRE(face_info.area_scale == 9);
    REQUIRE(face_info.length_scale == 3);
    REQUIRE(face_info.jacobian == 1.5);
    REQUIRE(face_info.normal == (Vec2<double>{0.0, 1.0}));
}

TEST_CASE("FacetInfo3D", "[dense_builder]") 
{
    Facet<3> f{
        Vec3<double>{0.0, 0.0, 0.0},
        Vec3<double>{2.0, 0.0, 0.0},
        Vec3<double>{0.0, 2.0, 0.0}
    };
    auto face_info = FacetInfo<3>::build(f);
    REQUIRE(face_info.area_scale == 2.0);
    REQUIRE(face_info.length_scale == std::sqrt(2.0));
    REQUIRE(face_info.jacobian == 4.0);
    REQUIRE(face_info.normal == (Vec3<double>{0.0, 0.0, 1.0}));
}

TEST_CASE("get_facet_info", "[dense_builder]") 
{
    auto m = line_mesh({0, 0}, {1, 0}).refine_repeatedly(4);
    auto f = get_facet_info(m);
    REQUIRE(f.size() == m.n_facets());
}

TEST_CASE("interior obs pts", "[dense builder]")
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

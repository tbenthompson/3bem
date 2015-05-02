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
#include "vectorx.h"

using namespace tbem;
struct EvalProb {
    EvalProb(int refine_level, int near_eval, int gauss_order,
             Vec3<double> center = Vec3<double>{0,0,0},
             double r = 3.0):
        sphere(sphere_mesh(center, r, refine_level)),
        qs(gauss_order, near_eval, 2.0, 1e-2),
        obs_pt(random_pt<3>()),
        obs_n(random_pt<3>()),
        obs_length_scale(std::sqrt(tri_area(sphere.facets[0]))),
        src_strength(std::vector<double>(sphere.n_dofs(), 1.0))
    {}

    double go(const Kernel<3,1,1>& k) {
        ObsPt<3> obs{obs_length_scale, obs_pt, obs_n, obs_n};
        auto mthd = make_adaptive_integration_mthd(qs, k);
        auto op = mesh_to_points_operator({obs}, sphere, mthd);
        return op.apply({src_strength})[0][0];
    }

    Mesh<3> sphere;
    QuadStrategy<3> qs;
    Vec3<double> obs_pt;
    Vec3<double> obs_n;
    double obs_length_scale;
    std::vector<double> src_strength;
}; 

TEST_CASE("EvalIntegralEquationSphereSurfaceArea", "[dense_builder]") {
    EvalProb ep(5, 3, 2);
    IdentityScalar<3> identity;
    double result = ep.go(identity);
    double exact_surf_area = 4*M_PI*9;
    REQUIRE_CLOSE(result, exact_surf_area, 1e-1);
}

TEST_CASE("ConstantLaplace", "[dense_builder]") {
    EvalProb ep(5, 3, 2);
    LaplaceDouble<3> double_kernel;
    double result = ep.go(double_kernel);
    REQUIRE_CLOSE(result, -1.0, 1e-3);
}

TEST_CASE("MatrixRowVsEval", "[dense_builder]") {
    EvalProb ep(4, 3, 2);
    double result = ep.go(LaplaceDouble<3>());
    REQUIRE_CLOSE(result, -1.0, 1e-3);
}

TEST_CASE("FacetInfo2D", "[dense_builder]") {
    Facet<2> f{Vec2<double>{0.0, 0.0}, Vec2<double>{3.0, 0.0}};
    auto face_info = FacetInfo<2>::build(f);
    REQUIRE(face_info.area_scale == 9);
    REQUIRE(face_info.length_scale == 3);
    REQUIRE(face_info.jacobian == 1.5);
    REQUIRE(face_info.normal == (Vec2<double>{0.0, 1.0}));
}

TEST_CASE("FacetInfo3D", "[dense_builder]") {
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

TEST_CASE("ObsPtFromFace", "[dense_builder]") {
    Facet<2> f{Vec2<double>{0.0, 0.0}, Vec2<double>{1.0, 1.0}};
    auto face_info = FacetInfo<2>::build(f);
    auto obs = ObsPt<2>::from_face({0}, face_info);
    REQUIRE(obs.len_scale == std::sqrt(2));
    REQUIRE(obs.loc == (Vec2<double>{0.5, 0.5}));
    REQUIRE(obs.normal == (Vec2<double>{-1.0 / std::sqrt(2), 1.0 / std::sqrt(2)}));
    REQUIRE(obs.richardson_dir == obs.normal);
}

void test_integral_operator(Mesh<2> m1, Mesh<2> m2) 
{
    QuadStrategy<2> qs(3);
    LaplaceDouble<2> k;
    auto mthd = make_adaptive_integration_mthd(qs, k);
    VectorX v(random_list(m2.n_dofs()));
    auto correct = dense_integral_operator(m1, m2, mthd).apply_scalar(v);
    auto other_op = integral_operator(m1, m2, mthd);
    auto other = other_op.apply_scalar(v);
    REQUIRE(other_op.n_total_rows() == m1.n_dofs());
    REQUIRE(other_op.n_total_cols() == m2.n_dofs());
    REQUIRE_ARRAY_CLOSE(correct, other, m1.n_dofs(), 1e-12);
}

TEST_CASE("IntegralOperatorSameMesh", "[dense_builder]") {
    auto m = circle_mesh({0, 0}, 1.0, 5);
    test_integral_operator(m, m);
}

TEST_CASE("IntegralOperatorDifferentMesh", "[dense_builder]") {
    auto m1 = circle_mesh({0, 0}, 1.0, 5);
    auto m2 = circle_mesh({1, 0}, 1.0, 4);
    test_integral_operator(m1, m2);
}

TEST_CASE("IntegralOperatorTensor", "[dense_builder]") {
    auto m1 = circle_mesh({0, 0}, 1.0, 5);
    auto m2 = m1;
    QuadStrategy<2> qs(3);
    ElasticHypersingular<2> k(1.0, 0.25);
    auto mthd = make_adaptive_integration_mthd(qs, k);
    BlockVectorX v{
        VectorX(random_list(m2.n_dofs())),
        VectorX(random_list(m2.n_dofs()))
    };
    auto correct = dense_integral_operator(m1, m2, mthd).apply(v);
    auto other_op = integral_operator(m1, m2, mthd);
    auto other = other_op.apply(v);
    REQUIRE(other_op.n_total_rows() == 2 * m1.n_dofs());
    REQUIRE(other_op.n_total_cols() == 2 * m2.n_dofs());
    REQUIRE_ARRAY_CLOSE(correct[0], other[0], m1.n_dofs(), 1e-12);
    REQUIRE_ARRAY_CLOSE(correct[1], other[1], m1.n_dofs(), 1e-12);
}

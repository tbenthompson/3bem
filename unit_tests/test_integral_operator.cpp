#include "catch.hpp"
#include "integral_operator.h"
#include "laplace_kernels.h"
#include "util.h"
#include "elastic_kernels.h"
#include "identity_kernels.h"
#include "mesh_gen.h"

using namespace tbem;

TEST_CASE("interior operator", "[dense_builder]") 
{
    auto sphere = sphere_mesh({0, 0, 0}, 3.0, 5);
    auto obs_pt = random_pt<3>();
    auto obs_n = random_pt<3>();
    std::vector<double> src_strength(sphere.n_dofs(), 1.0);

    SECTION("sphere surface area") {
        auto mthd = make_adaptive_integrator(1e-2, 2, 3, 2.0, IdentityScalar<3>());
        auto op = dense_interior_operator({obs_pt}, {obs_n}, sphere, mthd, sphere);
        auto result = op.apply(src_strength)[0];
        REQUIRE_CLOSE(result, 4 * M_PI * 9, 1e-1);
    }

    SECTION("constant laplace") {
        auto mthd = make_adaptive_integrator(1e-2, 2, 3, 2.0, LaplaceDouble<3>());
        auto op = dense_interior_operator({obs_pt}, {obs_n}, sphere, mthd, sphere);
        auto result = op.apply(src_strength)[0];
        REQUIRE_CLOSE(result, -1.0, 1e-3);
    }
}

void test_boundary_operator(Mesh<2> m1, Mesh<2> m2) 
{
    LaplaceDouble<2> k;
    auto mthd = make_adaptive_integrator(1e-4, 3, 8, 3.0, k);
    FMMConfig fmm_config{0.3, 30, 10000, 0.1, true};
    std::vector<double> v(random_list(m2.n_dofs()));
    auto correct = dense_boundary_operator(m1, m2, mthd, {m2}).apply(v);
    auto other_op = boundary_operator(m1, m2, mthd, fmm_config, {m2});
    auto other = other_op.apply(v);
    REQUIRE(other_op.n_rows() == m1.n_dofs());
    REQUIRE(other_op.n_cols() == m2.n_dofs());
    REQUIRE_ARRAY_CLOSE(correct, other, m1.n_dofs(), 1e-3);
}

TEST_CASE("IntegralOperatorSameMesh", "[boundary_operator]") 
{
    auto m = circle_mesh({0, 0}, 1.0, 2);
    test_boundary_operator(m, m);
}

TEST_CASE("IntegralOperatorDifferentMesh", "[boundary_operator]") 
{
    auto m1 = circle_mesh({0, 0}, 1.0, 5);
    auto m2 = circle_mesh({1, 0}, 1.0, 4);
    test_boundary_operator(m1, m2);
}

TEST_CASE("IntegralOperatorTensor", "[boundary_operator]") 
{
    auto m1 = circle_mesh({0, 0}, 1.0, 3);
    auto m2 = m1;
    ElasticHypersingular<2> k(1.0, 0.25);
    auto mthd = make_adaptive_integrator(1e-4, 3, 8, 3.0, k);
    FMMConfig fmm_config{0.3, 30, 10000, 0.1, true};
    auto v = random_list(2 * m2.n_dofs());
    auto correct = dense_boundary_operator(m1, m2, mthd, {m2}).apply(v);
    auto other_op = boundary_operator(m1, m2, mthd, fmm_config, {m2});
    auto other = other_op.apply(v);
    REQUIRE(other_op.n_rows() == 2 * m1.n_dofs());
    REQUIRE(other_op.n_cols() == 2 * m2.n_dofs());
    REQUIRE_ARRAY_CLOSE(correct, other, m1.n_dofs(), 1e-12);
}

TEST_CASE("nearfield matrix", "[boundary_operator]")
{
    auto m = circle_mesh({0, 0}, 1.0, 3);
    //far threshold of 1e10 forces nearfield matrix to be equivalent to the
    //full matrix
    auto mthd = make_adaptive_integrator(1e-4, 3, 8, 1e10, LaplaceDouble<2>());
    FMMConfig fmm_config{0.3, 30, 10000, 0.1, true};
    auto op = boundary_operator(m, m, mthd, fmm_config, m);
    auto v = random_list(m.n_dofs());

    auto nearfield_apply = op.nearfield.apply(v);
    auto full_apply = op.apply(v);

    REQUIRE_ARRAY_CLOSE(nearfield_apply, full_apply, nearfield_apply.size(), 1e-12);
}

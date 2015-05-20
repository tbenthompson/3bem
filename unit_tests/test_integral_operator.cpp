#include "catch.hpp"
#include "integral_operator.h"
#include "laplace_kernels.h"
#include "util.h"
#include "dense_builder.h"
#include "elastic_kernels.h"
#include "mesh_gen.h"

using namespace tbem;

void test_integral_operator(Mesh<2> m1, Mesh<2> m2) 
{
    QuadStrategy<2> qs(3);
    LaplaceDouble<2> k;
    auto mthd = make_adaptive_integration_mthd(qs, k);
    std::vector<double> v(random_list(m2.n_dofs()));
    auto correct = dense_integral_operator(m1, m2, mthd, {m2}).apply(v);
    auto other_op = integral_operator(m1, m2, mthd, {m2});
    auto other = other_op.apply(v);
    REQUIRE(other_op.n_rows() == m1.n_dofs());
    REQUIRE(other_op.n_cols() == m2.n_dofs());
    REQUIRE_ARRAY_CLOSE(correct, other, m1.n_dofs(), 1e-12);
}

TEST_CASE("IntegralOperatorSameMesh", "[integral_operator]") 
{
    auto m = circle_mesh({0, 0}, 1.0, 2);
    test_integral_operator(m, m);
}

TEST_CASE("IntegralOperatorDifferentMesh", "[integral_operator]") 
{
    auto m1 = circle_mesh({0, 0}, 1.0, 5);
    auto m2 = circle_mesh({1, 0}, 1.0, 4);
    test_integral_operator(m1, m2);
}

TEST_CASE("IntegralOperatorTensor", "[integral_operator]") 
{
    auto m1 = circle_mesh({0, 0}, 1.0, 5);
    auto m2 = m1;
    QuadStrategy<2> qs(3);
    ElasticHypersingular<2> k(1.0, 0.25);
    auto mthd = make_adaptive_integration_mthd(qs, k);
    auto v = random_list(2 * m2.n_dofs());
    auto correct = dense_integral_operator(m1, m2, mthd, {m2}).apply(v);
    auto other_op = integral_operator(m1, m2, mthd, {m2});
    auto other = other_op.apply(v);
    REQUIRE(other_op.n_rows() == 2 * m1.n_dofs());
    REQUIRE(other_op.n_cols() == 2 * m2.n_dofs());
    REQUIRE_ARRAY_CLOSE(correct, other, m1.n_dofs(), 1e-12);
}

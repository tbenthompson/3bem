#include "UnitTest++.h"
#include "integral_term.h"
#include "laplace_kernels.h"
#include "elastic_kernels.h"
#include "identity_kernels.h"
#include "dense_builder.h"
#include "mesh_gen.h"

using namespace tbem;

TEST(IdentityTensor) {
    IdentityTensor<3,3,3> K;
    (void)K;
}

TEST(IntegralOne) {
    QuadStrategy<2> quad_strategy(2);
    IdentityScalar<2> identity;
    ObsPt<2> obs{0.01, {0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}};
    auto facet_info = FacetInfo<2>::build({{{0,0},{1,0}}});
    auto term = make_integral_term(quad_strategy, identity, obs, facet_info);
    auto result = compute_term(term);
    CHECK_CLOSE(result[0][0][0], 0.5, 1e-6);
    CHECK_CLOSE(result[1][0][0], 0.5, 1e-6);
}

template <typename KT>
void integral_term_test(const KT& K, double claimed_distance, double exact) {
    QuadStrategy<3> quad_strategy(2);
    ObsPt<3> obs{0.01, {2.0, 2.0, 2.0}, {1.0, 0.0, 0.0}, {0.0, 0.0}};
    auto facet_info = FacetInfo<3>::build({{{0,0,0},{2,0,0},{0,1,0}}});
    auto term = make_integral_term(quad_strategy, K, obs, facet_info);
    auto result = compute_term(term);
    auto est = sum(result); 
    CHECK_CLOSE(est[0][0], exact, 1e-3);
}

TEST(IntegralLaplaceSingleFar) {
    LaplaceSingle<3> single_kernel;
    integral_term_test(single_kernel, 2.0, 0.0269063);
}

TEST(IntegralLaplaceSingleFakeNear) {
    LaplaceSingle<3> single_kernel;
    double exact = 0.0269063;
    integral_term_test(single_kernel, 1e-1, exact);
    integral_term_test(single_kernel, 1e-12, exact);
}

TEST(IntegralLaplaceDoubleFakeNear) {
    LaplaceDouble<3> double_kernel;
    double exact = -0.00621003;
    integral_term_test(double_kernel, 1.0, exact);
    integral_term_test(double_kernel, 1e-1, exact);
    integral_term_test(double_kernel, 1e-12, exact);
}

TEST(RichardsonExtrapolate) {
    // For n values in a sequence, Richardson extrapolation should be 
    // exact for polynomial sequences with degree n - 1.
    size_t n = 4;
    double error = 1e-6;
    std::vector<double> x = {1.0, 1.0 / 2, 1.0 / 4, 1.0 / 8};
    for (size_t i = 1; i < n; i++) {
        std::vector<double> input;
        for (size_t j = 0; j < n; j++) {
            input.push_back(std::pow(x[j], i) - 1.0);
        }
        double result = richardson_limit(input);
        CHECK_CLOSE(result, -1.0, error);
    }
}

TEST(RichardsonZeros) {
    std::vector<double> input = {0.0, 0.0, 0.0, 0.0};
    double result = richardson_limit(input);
    CHECK_CLOSE(result, 0.0, 1e-12);
}

TEST(TensorKernel) {
    ElasticDisplacement<2> k(30e9, 0.25);
    auto facet_info = FacetInfo<2>::build({{{-1.0, 0.0}, {1.0, 0.0}}});
    auto result = eval_point_influence({0.0}, k, facet_info, {0.0, 1.0}, {0.0, 1.0});
    CHECK_CLOSE(result[0][1][1], 8.84194e-13, 1e-17);
    CHECK_CLOSE(result[1][1][1], 8.84194e-13, 1e-17);
    CHECK_EQUAL(result[0][0][0], 0.0); CHECK_EQUAL(result[0][0][1], 0.0);
    CHECK_EQUAL(result[0][1][0], 0.0); CHECK_EQUAL(result[1][0][0], 0.0);
    CHECK_EQUAL(result[1][0][1], 0.0); CHECK_EQUAL(result[1][1][0], 0.0);
}

int main(int, char const *[])
{
    return UnitTest::RunAllTests();
}

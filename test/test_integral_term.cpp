#include "UnitTest++.h"
#include "integral_term.h"
#include "laplace_kernels.h"
#include "elastic_kernels.h"
#include "identity_kernels.h"
#include "test_shared.h"
#include "util.h"
#include "richardson.h"

using namespace tbem;

TEST(IdentityTensor) {
    IdentityTensor<3,3,3> K;
    (void)K;
}

TEST(IntegralOne) {
    QuadStrategy<2> qs(2);
    AdaptiveIntegrationMethod<2,1,1> mthd(qs);
    IdentityScalar<2> identity;
    ObsPt<2> obs{0.01, {0.0, 0.0}, {0.0, 0.0}, {0.0, 1.0}};
    auto facet_info = FacetInfo<2>::build({{{0,0},{1,0}}});
    IntegralTerm<2,1,1> term{identity, obs, facet_info};
    NearestPoint<2> nearest_pt{{0.0}, {0.0, 0.0}, 0.0, FarNearType::Farfield};
    auto result = mthd.compute_term(term, nearest_pt);
    CHECK_CLOSE(result[0][0][0], 0.5, 1e-6);
    CHECK_CLOSE(result[1][0][0], 0.5, 1e-6);
}

template <size_t dim, size_t R, size_t C>
void integral_term_test(const IntegrationMethodI<dim,R,C>& mthd,
    const Kernel<dim,R,C>& K, double distance, double exact) 
{
    QuadStrategy<3> quad_strategy(2);
    ObsPt<3> obs{0.01, {0.5, 0.1, distance}, {1.0, 0.0, 0.0}, {0.0, 0.0, 1.0}};
    auto facet_info = FacetInfo<3>::build({{{0,0,0},{2,0,0},{0,1,0}}});
    IntegralTerm<dim,R,C> term{K, obs, facet_info};
    auto nearest_pt = FarNearLogic<3>{3.0, 1.0}.decide(obs.loc, facet_info);
    auto result = mthd.compute_term(term, nearest_pt);
    auto est = sum(result); 
    CHECK_CLOSE(est[0][0], exact, 1e-3);
}

void integral_laplace_single(const IntegrationMethodI<3,1,1>& mthd) {
    LaplaceSingle<3> single_kernel;
    integral_term_test(mthd, single_kernel, 20.0, 0.00398);
    integral_term_test(mthd, single_kernel, 2.0, 0.0381);
    integral_term_test(mthd, single_kernel, 1e-1, 0.194);
    integral_term_test(mthd, single_kernel, 1e-6, 0.235);
}

TEST(IntegralLaplaceSingle) {
    QuadStrategy<3> qs(2);
    AdaptiveIntegrationMethod<3,1,1> mthd_adapt(qs);
    integral_laplace_single(mthd_adapt);
    SinhIntegrationMethod<3,1,1> mthd_sinh(qs);
    integral_laplace_single(mthd_sinh);
}

void integral_laplace_double(const IntegrationMethodI<3,1,1>& mthd) {
    LaplaceDouble<3> double_kernel;
    integral_term_test(mthd, double_kernel, 20.0, -0.00020);
    integral_term_test(mthd, double_kernel, 1.0, -0.0549);
    integral_term_test(mthd, double_kernel, 1e-1, -0.336);
    integral_term_test(mthd, double_kernel, 1e-6, -0.500);
}

TEST(IntegralLaplaceDouble) {
    QuadStrategy<3> qs(2);
    AdaptiveIntegrationMethod<3,1,1> mthd_adapt(qs);
    integral_laplace_double(mthd_adapt);
    SinhIntegrationMethod<3,1,1> mthd_sinh(qs);
    integral_laplace_double(mthd_sinh);
}

TEST(IntegralElasticDisplacement) {
    ElasticDisplacement<3> k(1.0, 0.25);
    QuadStrategy<3> qs(2);
    SinhIntegrationMethod<3,3,3> mthd(qs);
    integral_term_test(mthd, k, 20.0, 0.00265);
    integral_term_test(mthd, k, 1.0, 0.0495);
    integral_term_test(mthd, k, 1e-1, 0.1607);
    integral_term_test(mthd, k, 1e-6, 0.2014);
}

template <size_t dim, size_t R, size_t C>
void sinh_sufficient_accuracy(const Kernel<dim,R,C>& K) {
    QuadStrategy<dim> qs(2, 2, 8, 3.0, 1e-4);
    AdaptiveIntegrationMethod<dim,R,C> mthd_adapt(qs);
    SinhIntegrationMethod<dim,R,C> mthd_sinh(qs);

    double max_x = 1.0;
    double max_y = 1.0;
    auto facet_info = FacetInfo<dim>::build({{{0,0,0},{max_x,0,0},{0,max_y,0}}});
#pragma omp parallel for
    for (size_t ix = 1; ix < 10; ix++) {
        double x = (max_x / 10.0) * ix;
        if (x == 0.5) continue;
        for (double y_hat = 0.1; y_hat < 1.0; y_hat += 0.1) {
            double y = y_hat * max_y * ((max_x - x) / max_x);
            for (double log_z = 2; log_z > -6; log_z -= 1) {
                std::cout << Vec<double,3>{x, y, log_z} << std::endl;
                double z = std::pow(10.0, log_z);
                ObsPt<3> obs{0.01, {x, y, z}, {0.0, 0.0, 1.0}, {0.0, 0.0, 1.0}};
                IntegralTerm<dim,R,C> term{K, obs, facet_info};
                auto nearest_pt = FarNearLogic<dim>{3.0, 1.0}.decide(obs.loc, facet_info);
                auto sinh_eval = mthd_sinh.compute_term(term, nearest_pt);
                // auto adapt_eval = mthd_adapt.compute_term(term, nearest_pt);
                // auto error = fabs(sinh_eval - adapt_eval) / fabs(adapt_eval);
                // CHECK_CLOSE(sinh_eval, adapt_eval, 1e-3);
            }
        }
    }
}

TEST(SinhSufficientAccuracy) {
    // sinh_sufficient_accuracy(LaplaceSingle<3>());
    // sinh_sufficient_accuracy(ElasticTraction<3>(1.0, 0.25));
    sinh_sufficient_accuracy(ElasticHypersingular<3>(1.0, 0.25));
}

TEST(RichardsonExtrapolate) {
    // For n values in a sequence, Richardson extrapolation should be 
    // exact for polynomial sequences with degree n - 1.
    size_t n = 4;
    double allowed_error = 1e-6;
    std::vector<double> x = {
        1.0,
        1.0 / 2,
        1.0 / 4,
        1.0 / 8
    };
    for (size_t i = 1; i < n; i++) {
        std::vector<double> input;
        for (size_t j = 0; j < n; j++) {
            input.push_back(std::pow(x[j], i) - 1.0);
        }
        double result = richardson_limit(input);
        CHECK_CLOSE(result, -1.0, allowed_error);
    }
}

TEST(TensorKernel) {
    ElasticDisplacement<2> k(1.0, 0.25);
    QuadStrategy<2> qs(2);
    auto facet_info = FacetInfo<2>::build({{{-1.0, 0.0}, {1.0, 0.0}}});
    ObsPt<2> obs{0.1, {0.0, 1.0}, {0.0, 0.0}, {0.0, 0.0}};
    IntegralTerm<2,2,2> term{k, obs, facet_info};
    auto result = term.eval_point_influence({0.0});
    CHECK_CLOSE(result[0][1][1], 0.0265258, 1e-6);
    CHECK_CLOSE(result[1][1][1], 0.0265258, 1e-6);
    CHECK_EQUAL(result[0][0][0], 0.0); CHECK_EQUAL(result[0][0][1], 0.0);
    CHECK_EQUAL(result[0][1][0], 0.0); CHECK_EQUAL(result[1][0][0], 0.0);
    CHECK_EQUAL(result[1][0][1], 0.0); CHECK_EQUAL(result[1][1][0], 0.0);
}

int main(int, char const *[])
{
    return UnitTest::RunAllTests();
}

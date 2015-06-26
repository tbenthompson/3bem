#include "catch.hpp"
#include "mesh_gen.h"
#include "integral_term.h"
#include "laplace_kernels.h"
#include "elastic_kernels.h"
#include "identity_kernels.h"
#include "util.h"
#include "richardson.h"

using namespace tbem;

struct NearfieldIntegratorFake: public NearfieldIntegratorI<2,1,1> {
    size_t count = 0;
    std::vector<Vec<double,2>> obs_locs;

    virtual Vec<Vec<Vec<double,1>,1>,2> 
    compute_nearfield(const Kernel<2,1,1>&, const IntegralTerm<2,1,1>& term, 
        const NearestPoint<2>& np) 
    {
        (void)np;
        count++; 
        obs_locs.push_back(term.obs.loc);
        return zeros<Vec<Vec<Vec<double,1>,1>,2>>::make();
    }
};

size_t integral_count(const IntegrationStrategy<2,1,1>& integrator)
{
    auto fake = dynamic_cast<NearfieldIntegratorFake*>(
        integrator.nearfield_integrator.get()
    );
    return fake->count;
}

std::vector<Vec<double,2>> integral_locs(const IntegrationStrategy<2,1,1>& integrator)
{
    auto fake = dynamic_cast<NearfieldIntegratorFake*>(
        integrator.nearfield_integrator.get()
    );
    return fake->obs_locs;
}

IntegrationStrategy<2,1,1> fake_integral(const ObsPt<2>& pt)
{
    auto in = std::unique_ptr<NearfieldIntegratorI<2,1,1>>(new NearfieldIntegratorFake);
    auto integrator = make_integrator(in, 2, 8, 3.0, LaplaceSingle<2>());

    Facet<2> f{{{0, 0}, {1, 0}}};
    auto facet_info = FacetInfo<2>::build(f);
    IntegralTerm<2,1,1> term{pt, facet_info};
    integrator.compute_term(term);
    return integrator;
}

TEST_CASE("nearfield integral count", "[integral_term]") 
{
    SECTION("farfield") {
        REQUIRE(integral_count(fake_integral({{0, 10}, {0, 1}, {0, 0}})) == 0);
    }

    SECTION("nearfield") {
        REQUIRE(integral_count(fake_integral({{0, 2}, {0, 1}, {0, 0}})) == 1);
    }

    SECTION("very close but no limit") {
        REQUIRE(integral_count(fake_integral({{0, 0.1}, {0, 1}, {0, 0}})) == 1);
    }

    SECTION("very close with numerical noise in limit direction") {
        REQUIRE(integral_count(fake_integral({{0, 0.1}, {0, 1}, {0, 1e-12}})) == 1);
    }

    SECTION("singular") {
        REQUIRE(integral_count(fake_integral({{0, 0.0}, {0, 1}, {0, 1}})) == 8);
    }

    SECTION("nonsingular with limit") {
        REQUIRE(integral_count(fake_integral({{0, 0.1}, {0, 1}, {0, 1}})) == 8);
    }
}

TEST_CASE("nearfield observation locations", "[integral_term]") 
{
    SECTION("nearfield") {
        ObsPt<2> pt{{0, 0.1}, {0, 1}, {0, 0}};
        auto locs = integral_locs(fake_integral(pt));
        REQUIRE(locs.size() == 1);
        REQUIRE(locs[0] == pt.loc);
    }

    SECTION("singular") {
        ObsPt<2> pt{{0, 0.0}, {0, 1}, {0, 1}};
        auto locs = integral_locs(fake_integral(pt));
        REQUIRE(locs.size() == 8);
        for (size_t i = 0; i < 8; i++) {
            REQUIRE_CLOSE(locs[i][0], 0.0, 1e-15);
            REQUIRE_CLOSE(locs[i][1], 1.0 / std::pow(2, i), 1e-12);
        }
    }
}

template <size_t dim, size_t R, size_t C>
void integral_term_test(const IntegrationStrategy<dim,R,C>& mthd,
    double distance, double exact) 
{
    ObsPt<3> obs{{0.5, 0.1, distance}, {1.0, 0.0, 0.0}, {0.0, 0.0, 0.01}};
    auto facet_info = FacetInfo<3>::build({{{0,0,0},{2,0,0},{0,1,0}}});
    IntegralTerm<dim,R,C> term{obs, facet_info};
    auto result = mthd.compute_term(term);
    auto est = sum(result); 
    REQUIRE_CLOSE(est[0][0], exact, 1e-3);
}

void integral_laplace_single(const IntegrationStrategy<3,1,1>& mthd) 
{
    integral_term_test(mthd, 20.0, 0.00398);
    integral_term_test(mthd, 2.0, 0.0381);
    integral_term_test(mthd, 1e-1, 0.194);
    integral_term_test(mthd, 1e-6, 0.235);
}

TEST_CASE("IntegralLaplaceSingle", "[integral_term]") 
{
    LaplaceSingle<3> single_kernel;
    auto mthd_adapt = make_adaptive_integrator(1e-4, 2, 8, 3.0, single_kernel);
    integral_laplace_single(mthd_adapt);
    auto mthd_sinh = make_sinh_integrator(3, 2, 8, 3.0, single_kernel);
    integral_laplace_single(mthd_sinh);
}

void integral_laplace_double(const IntegrationStrategy<3,1,1>& mthd) 
{
    integral_term_test(mthd, 20.0, -0.00020);
    integral_term_test(mthd, 1.0, -0.0549);
    integral_term_test(mthd, 1e-1, -0.336);
    integral_term_test(mthd, 1e-6, -0.500);
}

TEST_CASE("IntegralLaplaceDouble", "[integral_term]") 
{
    LaplaceDouble<3> double_kernel;
    auto mthd_adapt = make_adaptive_integrator(1e-4, 2, 8, 3.0, double_kernel);
    integral_laplace_double(mthd_adapt);
    auto mthd_sinh = make_sinh_integrator(4, 2, 8, 3.0, double_kernel);
    integral_laplace_double(mthd_sinh);
}

TEST_CASE("IntegralElasticDisplacement", "[integral_term]") 
{
    ElasticDisplacement<3> k(1.0, 0.25);
    auto mthd = make_sinh_integrator(3, 2, 8, 3.0, k);
    integral_term_test(mthd, 20.0, 0.00265);
    integral_term_test(mthd, 1.0, 0.0495);
    integral_term_test(mthd, 1e-1, 0.1607);
    integral_term_test(mthd, 1e-6, 0.2014);
}

template <size_t dim, size_t R, size_t C>
void sinh_sufficient_accuracy(const Kernel<dim,R,C>& K) 
{
    double accuracy = 1e-2;
    auto mthd_adapt = make_adaptive_integrator(accuracy / 1000.0, 2, 8, 3.0, K);
    auto mthd_sinh = make_sinh_integrator(8, 2, 8, 3.0, K);

    auto facet = zeros<Vec<Vec<double,dim>,dim>>::make();
    for (size_t d = 0; d < dim - 1; d++) {
        facet[d + 1][d] = 1.0;
    }
    auto facet_info = FacetInfo<dim>::build(facet);
    for (size_t i = 0; i < 10; i++) {
        //generate random point on facet
        auto on_facet_pt = random_pt_facet<dim>(facet);
        //loop over log-scaled distances from the point
        for (double log_z = 2; log_z > -6; log_z -= 2) {
            double z = std::pow(10.0, log_z);
            
            auto loc = on_facet_pt + z * facet_info.normal;
            ObsPt<dim> obs{loc, facet_info.normal, facet_info.normal * 0.3};
            IntegralTerm<dim,R,C> term{obs, facet_info};

            auto sinh_eval = mthd_sinh.compute_term(term);
            auto adapt_eval = mthd_adapt.compute_term(term);
            auto error = fabs(sinh_eval - adapt_eval) / fabs(adapt_eval);
            size_t n_vals = dim * R * C;
            double l2_adapt = std::sqrt(
                sum(sum(sum(adapt_eval * adapt_eval))) / n_vals
            );
            double l2_error = std::sqrt(
                sum(sum(sum(error * error))) / n_vals
            );
            double l2_error_ratio = l2_error / l2_adapt;
            REQUIRE_CLOSE(l2_error_ratio, 0.0, 1e-2);
        }
    }
}

TEST_CASE("sinh sufficient accuracy 2D", "[integral_term]") 
{
    sinh_sufficient_accuracy(ElasticDisplacement<2>(1.0, 0.25));
    sinh_sufficient_accuracy(ElasticTraction<2>(1.0, 0.25));
    sinh_sufficient_accuracy(ElasticHypersingular<2>(1.0, 0.25));
    // sinh_sufficient_accuracy(ElasticDisplacement<3>(1.0, 0.25));
    // sinh_sufficient_accuracy(ElasticTraction<3>(1.0, 0.25));
    // sinh_sufficient_accuracy(ElasticAdjointTraction<3>(1.0, 0.25));
    // sinh_sufficient_accuracy(ElasticHypersingular<3>(1.0, 0.25));
}

TEST_CASE("tensor kernel", "[integral_term]") 
{
    ElasticDisplacement<2> k(1.0, 0.25);
    auto facet_info = FacetInfo<2>::build({{{-1.0, 0.0}, {1.0, 0.0}}});
    ObsPt<2> obs{{0.0, 1.0}, {0.0, 0.0}, {0.0, 0.0}};
    IntegralTerm<2,2,2> term{obs, facet_info};
    auto result = term.eval_point_influence(k, {0.0});
    REQUIRE_CLOSE(result[0][1][1], 0.0265258, 1e-6);
    REQUIRE_CLOSE(result[1][1][1], 0.0265258, 1e-6);
    REQUIRE(result[0][0][0] == 0.0);
    REQUIRE(result[0][0][1] == 0.0);
    REQUIRE(result[0][1][0] == 0.0);
    REQUIRE(result[1][0][0] == 0.0);
    REQUIRE(result[1][0][1] == 0.0);
    REQUIRE(result[1][1][0] == 0.0);
}

TEST_CASE("sinh integration -- scale shouldn't matter 2D", "[integral_term]")
{
    for (size_t steps = 2; steps < 10; steps++) {
        double L = std::pow(10, steps - 2);
        auto mthd = make_sinh_integrator(
            12, 2, 8, 3.0, LaplaceDouble<2>()
        );
        auto facet_info = FacetInfo<2>::build({{{0, 0}, {L, 0}}});
        ObsPt<2> obs{{L / 2, 0}, {0.0, 1.0}, {0.0, L / 3.0}};
        IntegralTerm<2,1,1> term{obs, facet_info};
        auto result = mthd.compute_term(term);
        REQUIRE_CLOSE(result[0][0][0], -0.25, 1e-7);
    }
}

TEST_CASE("sinh integration -- scale shouldn't matter 3D", "[integral_term]")
{
    for (size_t steps = 2; steps < 10; steps += 3) {
        double L = std::pow(30, steps - 2);
        auto mthd = make_sinh_integrator(
            5, 2, 8, 3.0, LaplaceDouble<3>()
        );
        auto facet_info = FacetInfo<3>::build({{{0, 0, 0}, {L, 0, 0}, {0, L, 0}}});
        ObsPt<3> obs{{L / 4, L / 4, 0}, {0.0, 0.0, 1.0}, {0.0, 0.0, L / 10.0}};
        IntegralTerm<3,1,1> term{obs, facet_info};
        auto result = mthd.compute_term(term);
        REQUIRE_CLOSE(result[0][0][0], -0.25, 1e-7);
    }
}

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
    auto integrator = make_integrator(in, 2, 2, 2, 8, 3.0, LaplaceSingle<2>());

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
        REQUIRE(integral_count(fake_integral({{0, 3.01}, {0, 1}, {0, 0}})) == 0);
    }

    SECTION("nearfield") {
        REQUIRE(integral_count(fake_integral({{0, 2.99}, {0, 1}, {0, 0}})) == 1);
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
            12, 2, 2, 2, 8, 3.0, LaplaceDouble<2>()
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
            5, 2, 2, 2, 8, 3.0, LaplaceDouble<3>()
        );
        auto facet_info = FacetInfo<3>::build({{{0, 0, 0}, {L, 0, 0}, {0, L, 0}}});
        ObsPt<3> obs{{L / 4, L / 4, 0}, {0.0, 0.0, 1.0}, {0.0, 0.0, L / 10.0}};
        IntegralTerm<3,1,1> term{obs, facet_info};
        auto result = mthd.compute_term(term);
        REQUIRE_CLOSE(result[0][0][0], -0.25, 1e-7);
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

TEST_CASE("IntegralLaplaceSingle", "[convergence]") 
{
    LaplaceSingle<3> single_kernel;
    auto mthd_adapt = make_adaptive_integrator(1e-4, 2, 2, 2, 8, 3.0, single_kernel);
    integral_laplace_single(mthd_adapt);
    auto mthd_sinh = make_sinh_integrator(3, 2, 2, 2, 8, 3.0, single_kernel);
    integral_laplace_single(mthd_sinh);
}

void integral_laplace_double(const IntegrationStrategy<3,1,1>& mthd) 
{
    integral_term_test(mthd, 20.0, -0.00020);
    integral_term_test(mthd, 1.0, -0.0549);
    integral_term_test(mthd, 1e-1, -0.336);
    integral_term_test(mthd, 1e-6, -0.500);
}

TEST_CASE("IntegralLaplaceDouble", "[convergence]") 
{
    LaplaceDouble<3> double_kernel;
    auto mthd_adapt = make_adaptive_integrator(1e-4, 2, 2, 2, 8, 3.0, double_kernel);
    integral_laplace_double(mthd_adapt);
    auto mthd_sinh = make_sinh_integrator(4, 2, 2, 2, 8, 3.0, double_kernel);
    integral_laplace_double(mthd_sinh);
}

TEST_CASE("IntegralElasticDisplacement", "[convergence]") 
{
    ElasticDisplacement<3> k(1.0, 0.25);
    auto mthd = make_sinh_integrator(3, 2, 2, 2, 8, 3.0, k);
    integral_term_test(mthd, 20.0, 0.00265);
    integral_term_test(mthd, 1.0, 0.0495);
    integral_term_test(mthd, 1e-1, 0.1607);
    integral_term_test(mthd, 1e-6, 0.2014);
}


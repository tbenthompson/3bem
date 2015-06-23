#include "catch.hpp"
#include "mesh_gen.h"
#include "integral_term.h"
#include "laplace_kernels.h"
#include "elastic_kernels.h"
#include "identity_kernels.h"
#include "util.h"
#include "richardson.h"

using namespace tbem;

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

TEST_CASE("IntegralOne", "[integral_term]") 
{
    IdentityScalar<2> identity;
    auto mthd = make_adaptive_integrator(1e-4, 2, 8, 3.0, identity);
    ObsPt<2> obs{{0.0, 0.0}, {0.0, 0.0}, {0.0, 0.01}};
    auto facet_info = FacetInfo<2>::build({{{0,0},{1,0}}});
    IntegralTerm<2,1,1> term{obs, facet_info};
    NearestPoint<2> nearest_pt{{0.0}, {0.0, 0.0}, 0.0, FarNearType::Farfield};
    auto result = mthd.compute_term(term, nearest_pt);
    REQUIRE_CLOSE(result[0][0][0], 0.5, 1e-6);
    REQUIRE_CLOSE(result[1][0][0], 0.5, 1e-6);
}

template <size_t dim, size_t R, size_t C>
void integral_term_test(const IntegrationStrategy<dim,R,C>& mthd,
    double distance, double exact) 
{
    ObsPt<3> obs{{0.5, 0.1, distance}, {1.0, 0.0, 0.0}, {0.0, 0.0, 0.01}};
    auto facet_info = FacetInfo<3>::build({{{0,0,0},{2,0,0},{0,1,0}}});
    IntegralTerm<dim,R,C> term{obs, facet_info};
    auto nearest_pt = FarNearLogic<3>{3.0, 1.0}.decide(obs.loc, facet_info);
    auto result = mthd.compute_term(term, nearest_pt);
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
    auto mthd_sinh = make_sinh_integrator(10, 2, 8, 3.0, single_kernel);
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
    auto mthd_sinh = make_sinh_integrator(10, 2, 8, 3.0, double_kernel);
    integral_laplace_double(mthd_sinh);
}

TEST_CASE("IntegralElasticDisplacement", "[integral_term]") 
{
    ElasticDisplacement<3> k(1.0, 0.25);
    auto mthd = make_sinh_integrator(10, 2, 8, 3.0, k);
    integral_term_test(mthd, 20.0, 0.00265);
    integral_term_test(mthd, 1.0, 0.0495);
    integral_term_test(mthd, 1e-1, 0.1607);
    integral_term_test(mthd, 1e-6, 0.2014);
}

template <size_t dim, size_t R, size_t C>
void sinh_sufficient_accuracy(const Kernel<dim,R,C>& K) 
{
    auto mthd_adapt = make_adaptive_integrator(1e-5, 2, 8, 3.0, K);
    auto mthd_sinh = make_sinh_integrator(10, 2, 8, 3.0, K);

    double max_x = 1.0;
    double max_y = 1.0;
    auto facet_info = FacetInfo<dim>::build({{{0,0,0},{max_x,0,0},{0,max_y,0}}});
    for (size_t ix = 1; ix < 10; ix++) {
        double x = (max_x / 10.0) * ix;
        for (double y_hat = 0.1; y_hat < 1.0; y_hat += 0.1) {
            double y = y_hat * max_y * ((max_x - x) / max_x);
            for (double log_z = 2; log_z > -6; log_z -= 1) {
                double z = std::pow(10.0, log_z);
                ObsPt<3> obs{{x, y, z}, {0.0, 0.0, 1.0}, {0.0, 0.0, 1.0}};
                IntegralTerm<dim,R,C> term{obs, facet_info};
                auto nearest_pt = FarNearLogic<dim>{3.0, 1.0}.decide(obs.loc, facet_info);
                auto sinh_eval = mthd_sinh.compute_term(term, nearest_pt);
                auto adapt_eval = mthd_adapt.compute_term(term, nearest_pt);
                auto error = fabs(sinh_eval - adapt_eval) / fabs(adapt_eval);
                REQUIRE_ARRAY_CLOSE(
                    &sinh_eval[0][0][0], &adapt_eval[0][0][0], dim * R * C, 1e-3
                );
            }
        }
    }
}

TEST_CASE("sinh sufficient accuracy", "[integral_term]") 
{
    // sinh_sufficient_accuracy(LaplaceSingle<3>());
    // sinh_sufficient_accuracy(ElasticTraction<3>(1.0, 0.25));
    /* sinh_sufficient_accuracy(ElasticHypersingular<3>(1.0, 0.25)); */
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

TEST_CASE("kernel inline in mthd creation call", "[integral_term]")
{
    // When kernels were not copied, this test would segfault because the
    // LaplaceDouble<3> object would go out of scope after the make_..._mthd
    // call and be destroyed
    auto mthd_adapt = make_adaptive_integrator(
        1e-4, 2, 8, 3.0, LaplaceDouble<3>()
    );
    integral_laplace_double(mthd_adapt);
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
        auto nearest_pt = FarNearLogic<2>{3.0, 1.0}.decide(obs.loc, facet_info);
        auto result = mthd.compute_term(term, nearest_pt);
        REQUIRE_CLOSE(result[0][0][0], -0.25, 1e-7);
    }
}

TEST_CASE("sinh integration -- scale shouldn't matter 3D", "[integral_term]")
{
    for (size_t steps = 2; steps < 10; steps++) {
        double L = std::pow(30, steps - 2);
        auto mthd = make_sinh_integrator(
            5, 2, 8, 3.0, LaplaceDouble<3>()
        );
        auto facet_info = FacetInfo<3>::build({{{0, 0, 0}, {L, 0, 0}, {0, L, 0}}});
        ObsPt<3> obs{{L / 4, L / 4, 0}, {0.0, 0.0, 1.0}, {0.0, 0.0, L / 10.0}};
        IntegralTerm<3,1,1> term{obs, facet_info};
        auto nearest_pt = FarNearLogic<3>{3.0, 1.0}.decide(obs.loc, facet_info);
        auto result = mthd.compute_term(term, nearest_pt);
        REQUIRE_CLOSE(result[0][0][0], -0.25, 1e-7);
    }
}

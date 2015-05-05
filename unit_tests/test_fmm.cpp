#include "catch.hpp"
#include "fmm.h"
#include "laplace_kernels.h"
#include "elastic_kernels.h"
#include "identity_kernels.h"
#include "nbody_operator.h"

using namespace tbem;

template <size_t dim>
NBodyData<dim> ones_data(size_t n) 
{
    auto src_pts = random_pts<dim>(n);
    auto obs_pts = random_pts<dim>(n);
    auto normals = random_pts<dim>(n);
    std::vector<double> weights(n, 1.0);
    return NBodyData<dim>{obs_pts, normals, src_pts, normals, weights};
}

TEST_CASE("MakeSurroundingSurface", "[fmm]") 
{
    auto surface = make_surrounding_surface<2>(4);
    REQUIRE_CLOSE(surface.pts[0], (Vec<double,2>{1.0, 0.0}), 1e-12);
    REQUIRE_CLOSE(surface.pts[1], (Vec<double,2>{0.0, 1.0}), 1e-12);
    REQUIRE_CLOSE(surface.pts[3], (Vec<double,2>{0.0, -1.0}), 1e-12);
    REQUIRE_CLOSE(surface.normals[0], (Vec<double,2>{1.0, 0.0}), 1e-12);
    REQUIRE_CLOSE(surface.normals[1], (Vec<double,2>{0.0, 1.0}), 1e-12);
} 

TEST_CASE("IdentityP2M", "[fmm]") 
{
    size_t n = 500;
    std::vector<double> x(n, 1.0);
    std::vector<double> weights(n, 1.0);
    NBodyData<2> data{
        random_pts<2>(n, 1, 2), random_pts<2>(n, -1, 1),
        random_pts<2>(n, -1, 1), random_pts<2>(n, -1, 1),
        weights
    };
    IdentityScalar<2> K;
    FMMOperator<2,1,1> tree(K, data, {0.3, 1, 1, 0.05});
    CheckToEquiv up_check_to_equiv;
    CheckToEquiv down_check_to_equiv;
    tree.build_check_to_equiv(tree.src_oct, up_check_to_equiv, down_check_to_equiv);

    double M_coeff;
    tree.P2M(tree.src_oct, up_check_to_equiv[0], x, &M_coeff);
    REQUIRE(M_coeff == static_cast<double>(n));

    double L_coeff = 0.0;
    tree.M2L(tree.obs_oct, tree.src_oct, down_check_to_equiv[0], &M_coeff, &L_coeff);
    REQUIRE(L_coeff == static_cast<double>(n));

    std::vector<double> L_child(4);
    std::vector<double*> child_ptr(4);
    for (size_t i = 0; i < 4; i++) {
        child_ptr[i] = &L_child[i];
    }
    tree.L2L(tree.src_oct, down_check_to_equiv[1], &L_coeff, child_ptr);
    for (size_t i = 0; i < 4; i++) {
        REQUIRE(L_child[i] == static_cast<double>(n));
    }

    std::vector<double> out(n, 0.0);
    tree.L2P(tree.src_oct, &L_coeff, out);
    for (size_t i = 0; i < n; i++) {
        REQUIRE(out[i] == static_cast<double>(n));
    }
}

template <size_t dim, size_t R, size_t C>
void test_kernel(const NBodyData<dim>& data, const Kernel<dim,R,C>& K,
    size_t order, double allowed_error) 
{
    size_t n_per_cell = 1;
    std::vector<double> x(C * data.src_locs.size(), 1.0);
    FMMOperator<dim,R,C> tree(K, data, {0.35, order, n_per_cell, 0.05});
    TIC
    auto out = tree.apply(x);
    TOC("FMM");

    BlockDirectNBodyOperator<dim,R,C> exact_op{data, K};
    auto exact = exact_op.apply(x);
    for (size_t d = 0; d < R; d++) {
        double average_magnitude = 0.0;
        for (size_t i = 0; i < data.obs_locs.size(); i++) {
            average_magnitude += exact[d * data.obs_locs.size() + i];
        }
        average_magnitude /= data.obs_locs.size();
        for (size_t i = 0; i < data.obs_locs.size(); i++) {
            auto out_val = out[d * data.obs_locs.size() + i];
            auto exact_val = exact[d * data.obs_locs.size() + i];
            auto error1 = std::fabs((out_val - exact_val) / exact_val);
            auto error2 = std::fabs((out_val - exact_val) / average_magnitude);
            auto error = std::min(error1, error2);
            // CHECK_CLOSE(error, 0, allowed_error);
        }
    }
}


template <size_t dim, size_t R, size_t C>
void test_kernel(const Kernel<dim,R,C>& K, size_t order, double allowed_error) 
{
    size_t n = 500;
    auto data = ones_data<dim>(n);
    return test_kernel(data, K, order, allowed_error);    
}

TEST_CASE("IdentityFMM", "[fmm]") 
{
    test_kernel(IdentityScalar<2>(), 1, 1e-4);
}

TEST_CASE("SingleLayer2DFMM", "[fmm]") 
{
    test_kernel(LaplaceSingle<2>(), 5, 1e-2);
}

TEST_CASE("DoubleLayer2DFMM", "[fmm]") 
{
    test_kernel(LaplaceDouble<2>(), 20, 1e-4);
}

TEST_CASE("HypersingularLayer2DFMM", "[fmm]") 
{
    test_kernel(LaplaceHypersingular<2>(), 20, 1e-4);
}

TEST_CASE("SingleLayer3DFMM", "[fmm]") 
{
    test_kernel(LaplaceSingle<3>(), 20, 1e-4);
}

TEST_CASE("DoubleLayer3DFMM", "[fmm]") 
{
    test_kernel(LaplaceDouble<3>(), 80, 1e-4);
}

TEST_CASE("ElasticDisplacement2DFMM", "[fmm]")
{
    test_kernel(ElasticDisplacement<2>(30e9, 0.25), 10, 1e-4);
}
TEST_CASE("ElasticTraction2DFMM", "[fmm]")
{
    test_kernel(ElasticTraction<2>(30e9, 0.25), 30, 1e-4);
}
TEST_CASE("ElasticHypersingular2DFMM", "[fmm]")
{
    test_kernel(ElasticHypersingular<2>(30e9, 0.25), 25, 1e-4);
}

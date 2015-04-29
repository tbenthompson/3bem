#include "UnitTest++.h"
#include "fmm.h"
#include "laplace_kernels.h"
#include "elastic_kernels.h"
#include "identity_kernels.h"
#include "nbody_operator.h"
#include "test_shared.h"

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

TEST(MakeSurroundingSurface) 
{
    auto surface = make_surrounding_surface<2>(4);
    CHECK_CLOSE(surface.pts[0], (Vec<double,2>{1.0, 0.0}), 1e-12);
    CHECK_CLOSE(surface.pts[1], (Vec<double,2>{0.0, 1.0}), 1e-12);
    CHECK_CLOSE(surface.pts[3], (Vec<double,2>{0.0, -1.0}), 1e-12);
    CHECK_CLOSE(surface.normals[0], (Vec<double,2>{1.0, 0.0}), 1e-12);
    CHECK_CLOSE(surface.normals[1], (Vec<double,2>{0.0, 1.0}), 1e-12);
} 

TEST(IdentityP2M) 
{
    size_t n = 500;
    BlockVectorX x(1, VectorX(n, 1.0));
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
    CHECK_EQUAL(M_coeff, static_cast<double>(n));

    double L_coeff = 0.0;
    tree.M2L(tree.obs_oct, tree.src_oct, down_check_to_equiv[0], &M_coeff, &L_coeff);
    CHECK_EQUAL(L_coeff, static_cast<double>(n));

    std::vector<double> L_child(4);
    std::vector<double*> child_ptr(4);
    for (size_t i = 0; i < 4; i++) {
        child_ptr[i] = &L_child[i];
    }
    tree.L2L(tree.src_oct, down_check_to_equiv[1], &L_coeff, child_ptr);
    for (size_t i = 0; i < 4; i++) {
        CHECK_EQUAL(L_child[i], static_cast<double>(n));
    }

    BlockVectorX out(1, VectorX(n, 0.0));
    tree.L2P(tree.src_oct, &L_coeff, out);
    for (size_t i = 0; i < n; i++) {
        CHECK_EQUAL(out[0][i], static_cast<double>(n));
    }
}

template <size_t dim, size_t R, size_t C>
void test_kernel(const NBodyData<dim>& data, const Kernel<dim,R,C>& K,
    size_t order, double allowed_error) 
{
    size_t n_per_cell = std::max<size_t>(20, order);
    BlockVectorX x(C, VectorX(data.src_weights));
    FMMOperator<dim,R,C> tree(K, data, {0.35, order, n_per_cell, 0.05});
    TIC
    auto out = tree.apply(x);
    TOC("FMM");

    BlockDirectNBodyOperator<dim,R,C> exact_op{data, K};
    auto exact = exact_op.apply(x);
    for (size_t d = 0; d < R; d++) {
        double average_magnitude = 0.0;
        for (size_t i = 0; i < data.obs_locs.size(); i++) {
            average_magnitude += exact[d][i];
        }
        average_magnitude /= data.obs_locs.size();
        for (size_t i = 0; i < data.obs_locs.size(); i++) {
            auto error1 = std::fabs((out[d][i] - exact[d][i]) / exact[d][i]);
            auto error2 = std::fabs((out[d][i] - exact[d][i]) / average_magnitude);
            auto error = std::min(error1, error2);
            CHECK_CLOSE(error, 0, allowed_error);
        }
    }
}


template <size_t dim, size_t R, size_t C>
void test_kernel(const Kernel<dim,R,C>& K, size_t order, double allowed_error) 
{
    size_t n = 10000;
    auto data = ones_data<dim>(n);
    return test_kernel(data, K, order, allowed_error);    
}

TEST(IdentityFMM) 
{
    test_kernel(IdentityScalar<2>(), 1, 1e-4);
}

TEST(SingleLayer2DFMM) 
{
    test_kernel(LaplaceSingle<2>(), 15, 1e-2);
}

TEST(DoubleLayer2DFMM) 
{
    test_kernel(LaplaceDouble<2>(), 20, 1e-4);
}

TEST(HypersingularLayer2DFMM) 
{
    test_kernel(LaplaceHypersingular<2>(), 20, 1e-4);
}

TEST(SingleLayer3DFMM) 
{
    test_kernel(LaplaceSingle<3>(), 20, 1e-4);
}

TEST(DoubleLayer3DFMM) 
{
    test_kernel(LaplaceDouble<3>(), 80, 1e-4);
}

TEST(ElasticDisplacement2DFMM)
{
    test_kernel(ElasticDisplacement<2>(30e9, 0.25), 10, 1e-4);
}
TEST(ElasticTraction2DFMM)
{
    test_kernel(ElasticTraction<2>(30e9, 0.25), 30, 1e-4);
}
TEST(ElasticHypersingular2DFMM)
{
    test_kernel(ElasticHypersingular<2>(30e9, 0.25), 25, 1e-4);
}

int main(int, char const *[])
{
    // return UnitTest::RunAllTests();
    return RunOneTest("SingleLayer2DFMM");
}

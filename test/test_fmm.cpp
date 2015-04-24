#include "UnitTest++.h"
#include "fmm.h"
#include "laplace_kernels.h"
#include "elastic_kernels.h"
#include "identity_kernels.h"
#include "nbody_operator.h"
#include "test_shared.h"

using namespace tbem;

TEST(LU) {
    std::vector<double> matrix{
        2, 1, -1, 0.5
    };
    auto lu = LU_decompose(matrix);
    auto soln = LU_solve(lu, {1,1});
    std::vector<double> correct{
        -0.25, 1.5
    };
    CHECK_ARRAY_CLOSE(soln, correct, 2, 1e-14);
}

template <size_t dim>
NBodyData<dim> ones_data(size_t n) 
{
    auto src_pts = random_pts<dim>(n);
    auto obs_pts = random_pts<dim>(n);
    auto normals = random_pts<dim>(n);
    std::vector<double> weights(n, 1.0);
    return NBodyData<dim>{src_pts, normals, obs_pts, normals, weights};
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

template <size_t dim, size_t R, size_t C>
void test_kernel(const Kernel<dim,R,C>& K, size_t order, double allowed_error) 
{
    size_t n = 10000;
    size_t n_per_cell = std::max<size_t>(50, order);
    auto data = ones_data<dim>(n);
    BlockVectorX x(C, VectorX(data.src_weights));
    TreeNBodyOperator<dim,R,C> tree(K, data, n_per_cell, order, 3.5);
    TIC
    auto out = tree.apply(x);
    TOC("ABAS");

    BlockDirectNBodyOperator<dim,R,C> exact_op{data, K};
    auto exact = exact_op.apply(x);
    for (size_t d = 0; d < R; d++) {
        for (size_t i = 0; i < n; i++) {
            auto error = std::fabs((out[d][i] - exact[d][i]) / exact[d][i]);
            CHECK_CLOSE(error, 0, allowed_error);
        }
    }
}

TEST(SingleLayer2DFMM) 
{
    test_kernel(LaplaceSingle<2>(), 15, 1e-4);
}

TEST(DoubleLayer2DFMM) 
{
    test_kernel(LaplaceDouble<2>(), 30, 1e-4);
}

TEST(HypersingularLayer2DFMM) 
{
    test_kernel(LaplaceHypersingular<2>(), 30, 1e-4);
}

TEST(SingleLayer3DFMM) 
{
    test_kernel(LaplaceSingle<3>(), 35, 1e-4);
}

TEST(DoubleLayer3DFMM) 
{
    test_kernel(LaplaceDouble<3>(), 200, 1e-4);
}

TEST(ElasticDisplacement2DFMM)
{
    test_kernel(ElasticDisplacement<2>(30e9, 0.25), 15, 1e-4);
}
TEST(ElasticTraction2DFMM)
{
    test_kernel(ElasticTraction<2>(30e9, 0.25), 40, 1e-4);
}
TEST(ElasticHypersingular2DFMM)
{
    test_kernel(ElasticHypersingular<2>(30e9, 0.25), 40, 1e-4);
}

int main(int, char const *[])
{
    return UnitTest::RunAllTests();
}

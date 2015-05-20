#include "catch.hpp"
#include "nearfield_operator.h"
#include "continuity_builder.h"
#include "quadrature.h"
#include "laplace_kernels.h"
#include "mesh_gen.h"
#include "integral_operator.h"

using namespace tbem;

TEST_CASE("nearfield obs pts", "[nearfield_operator]")
{
    auto m = circle_mesh({0, 0}, 1.0, 1);
    REQUIRE(m.facets.size() == 8);
    auto pts = nearfield_obs_pts(m, gauss(2), m);
    REQUIRE(pts.size() == 16);
}

TEST_CASE("Constrained nearfield matrix", "[nearfield_operator]")
{
    auto m = circle_mesh({0, 0}, 1.0, 2);
    // Use a very large nearfield range (300000) to force everything to
    // be nearfield.
    QuadStrategy<2> qs(3, 8, 300000, 1e-4);
    LaplaceDouble<2> k;
    auto mthd = make_adaptive_integration_mthd(qs, k);
    auto galerkin = make_galerkin_operator(1, m, mthd.get_obs_quad());
    auto nearfield = make_nearfield_operator(m, m, mthd, {m});
    auto dense_matrix = galerkin.right_multiply_with_dense(nearfield.to_dense());
    std::vector<MatrixEntry> sparsified;
    for (size_t i = 0; i < dense_matrix.n_rows(); i++) {
        for (size_t j = 0; j < dense_matrix.n_cols(); j++) {
            auto val = dense_matrix[i * dense_matrix.n_cols() + j];
            if (val == 0.0) {
                continue;
            }
            sparsified.push_back({i, j, val});
        }
    }
    auto matrix = SparseOperator::csr_from_coo(
        dense_matrix.n_rows(), dense_matrix.n_cols(), sparsified
    );

    SECTION("Dense operator equals galerkin nearfield") {
        auto dense_matrix2 = *dense_integral_operator(m, m, mthd, {m}).storage;
        REQUIRE_ARRAY_CLOSE(dense_matrix2, dense_matrix, dense_matrix2.size(), 1e-12);
    }

    SECTION("Dense condensation and sparse condensation are the same") {
        auto cm = from_constraints(
            convert_to_constraints(mesh_continuity(m.begin()))
        );
        auto test = condense_matrix(cm, cm, matrix).to_dense().data();
        auto correct = condense_matrix(cm, cm, dense_matrix).data();
        REQUIRE(test.size() == correct.size());
        REQUIRE_ARRAY_CLOSE(test, correct, correct.size(), 1e-12);
    }
}

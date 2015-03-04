#include "UnitTest++.h"
#include "matrix_entry.h"
#include "petsc_facade.h"
#include "mesh_gen.h"
#include "util.h"
#include "identity_kernels.h"
#include "dense_builder.h"
#include "matrix_free_farfield_builder.h"
#include "laplace_kernels.h"

using namespace tbem;

TEST(MatrixFreeMassOperator) {
    auto m = circle_mesh({0.0, 0.0}, 1.0, 4);
    BlockVectorX in{random_list(m.n_dofs()), random_list(m.n_dofs())};
    IdentityTensor<2,2,2> id;
    QuadStrategy<2> qs(2);
    auto problem = make_problem(m, m, id);
    auto dense_op = mass_operator(problem, qs);
    auto correct = dense_op.apply({in});
    auto mf_op = matrix_free_mass_operator(problem, qs);
    auto to_test = mf_op.apply({in});
    for (size_t d = 0; d < to_test.size(); d++) {
        CHECK_ARRAY_CLOSE(correct[d], to_test[d], to_test[d].size(), 1e-12);
    }
}

TEST(MatrixFreeIntegralOperator) {
    auto m = circle_mesh({0.0, 0.0}, 1.0, 4);
    LaplaceDouble<2> k;
    QuadStrategy<2> qs(2);
    auto problem = make_problem(m, m, k);
    auto op = make_matrix_free(problem, qs);
    auto dense_op = mesh_to_mesh_operator(problem, qs);

    BlockVectorX in{random_list(m.n_dofs())};
    auto to_test = op.apply(in);
    auto correct = dense_op.apply(in);
    for (size_t d = 0; d < to_test.size(); d++) {
        CHECK_ARRAY_CLOSE(correct[d], to_test[d], to_test[d].size(), 1e-12);
    }
}

int main() {
    return UnitTest::RunAllTests();
}

#include "UnitTest++.h"
#include "mesh_gen.h"
#include "util.h"
#include "identity_kernels.h"
#include "dense_builder.h"
#include "matrix_free_farfield_builder.h"
#include "mass_operator.h"
#include "laplace_kernels.h"

using namespace tbem;

TEST(MassTerm) {
    auto sphere = sphere_mesh({0,0,0}, 1.0, 0);
    std::vector<double> str(sphere.n_dofs(), 1.0);
    for (std::size_t i = 0; i < sphere.n_facets(); i++) {
        for (int d = 0; d < 3; d++) {
            if (sphere.facets[i][d][0] > 0.5) {
                str[3 * i + d] = 0.0;
            }
        }
    }
    IdentityScalar<3> identity;
    auto p = make_boundary_integral<3>(sphere, sphere, identity);
    QuadStrategy<3> qs(2);
    auto mass_op = mass_operator(p, qs);
    auto res = mass_op.apply({str})[0];
    CHECK_EQUAL(res.size(), sphere.n_dofs());
    double true_area = 0.0;
    for (auto f: sphere.facets) {
        true_area += tri_area(f);
    }
    double mass_area = 0.0;
    for (auto r: res) {
        mass_area += r;
    }   
    CHECK_CLOSE(mass_area, (5.0 / 6.0) * true_area, 1e-12);
}

TEST(TensorMassTerm) {
    auto sphere = sphere_mesh({0,0,0}, 1.0, 1);
    BlockVectorX str(3, VectorX(sphere.n_dofs(), 1.0));
    IdentityTensor<3,3,3> identity;
    auto p = make_boundary_integral<3>(sphere, sphere, identity);
    QuadStrategy<3> qs(2);
    auto mass_op = mass_operator(p, qs);
    CHECK_EQUAL(mass_op.ops.size(), 9);
    auto res = mass_op.apply({str});
    CHECK_EQUAL(res.size(), 3);
}

TEST(MatrixFreeIntegralOperator) {
    auto m = circle_mesh({0.0, 0.0}, 1.0, 4);
    LaplaceDouble<2> k;
    QuadStrategy<2> qs(2);
    auto problem = make_boundary_integral(m, m, k);
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

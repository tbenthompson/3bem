#include "UnitTest++.h"
#include "constraint.h"
#include "mesh.h"
#include "mesh_gen.h"
#include <iostream>

using namespace tbem;

TEST(ConstraintsAreCreated) {
    auto c = continuity_constraint(1, 2);
    Constraint correct = {
        {DOFWeight{1, 1.0f}, DOFWeight{2, -1.0f}},
        0.0
    };
    CHECK(c.dof_constraints == correct.dof_constraints);
    CHECK(c.rhs_value == correct.rhs_value);

    auto c2 = offset_constraint(3, 4, 5.0f);
    Constraint correct2 = {
        {DOFWeight{3, 1.0f}, DOFWeight{4, -1.0f}},
        -5.0
    };
    CHECK(c2.dof_constraints == correct2.dof_constraints);
    CHECK(c2.rhs_value == correct2.rhs_value);
}

TEST(ConstraintMatrixIsAppendedTo) {
    auto c = continuity_constraint(1, 2);
    auto cm = ConstraintMatrix::from_constraints({c});
    auto res = cm.c_map.at(2);
    CHECK(res.dof_constraints[0].first == 1);
    CHECK(res.dof_constraints[0].second == 1);
    CHECK(res.dof_constraints[1].first == 2);
    CHECK(res.dof_constraints[1].second == -1);
    CHECK(res.rhs_value == 0.0);
}

TEST(ConstraintMatrixGetAll) {
    auto c0 = boundary_condition(1, 4.0);
    auto c1 = continuity_constraint(1, 2);
    auto c2 = continuity_constraint(2, 3);
    auto cm = ConstraintMatrix::from_constraints({c0, c1, c2});
    auto in = cm.get_reduced({2.0, 4.0, 4.0, 4.0});
    auto res = cm.get_all(in, 4);
    double res_exact[4] = {2.0, 4.0, 4.0, 4.0};
    CHECK_ARRAY_CLOSE(res, res_exact, 4, 1e-13);
}

TEST(ConstraintMesh) {
    auto sphere = sphere_mesh({0, 0, 0}, 1).refine_repeatedly(2);
    auto constraints = mesh_continuity<3>(sphere);
    auto cons_mat = ConstraintMatrix::from_constraints(constraints);
    CHECK_EQUAL(3 * sphere.facets.size(), 384);
    CHECK_EQUAL(cons_mat.c_map.size(), 318);
    auto my_c = cons_mat.c_map.begin()->second.dof_constraints;
    CHECK_EQUAL(my_c[0].second, 1);
    CHECK_EQUAL(my_c[1].second, -1);
}

TEST(Condense) {
    auto sphere = sphere_mesh({0, 0, 0}, 1).refine_repeatedly(0);
    auto constraints = mesh_continuity<3>(sphere);
    auto cons_mat = ConstraintMatrix::from_constraints(constraints);
    std::vector<double> all_dofs(3 * sphere.facets.size(), 0.0);
    for (std::size_t i = 0; i < all_dofs.size(); i++) {
        cons_mat.add_vec_with_constraints({i, 1.0}, all_dofs);
    }
    int n_vertices = 0;
    for (auto c: all_dofs) {
        // Either a constraint is condensed out and so c == 1
        // or a constraint represents a vertex with 4 faces and c == 4
        CHECK(c == 0 || c == 4);
        if (c == 4) {
            n_vertices++;
        }
    }
    CHECK_EQUAL(n_vertices, 6);
}

TEST(MeshContinuity2D) {
    auto circle = circle_mesh({0,0},1).refine_repeatedly(4);
    auto constraints =
        ConstraintMatrix::from_constraints(mesh_continuity<2>(circle));
    std::vector<double> reduced(circle.facets.size(), 1.0);
    auto all_vec = constraints.get_all(reduced, 2 * circle.facets.size());
    for (auto a: all_vec) {
        CHECK_CLOSE(a, 1.0, 1e-12);
    }
    for (auto it = constraints.c_map.begin(); it != constraints.c_map.end(); it++) {
        // Of a pair of DOFs, the higher index should be even,
        // Because the last DOF is constrained, 
        // the constrained dofs should be even.

        // But, we don't check the last DOF which will be constrained by the first
        // DOF (connecting the circle) since then the evenness will be broken.
        if (it->first != (int)(2 * circle.facets.size() - 1)) {
            CHECK_EQUAL(it->first % 2, 0);
        }
    }
}

int main(int, char const *[])
{
    return UnitTest::RunAllTests();
}

#include "UnitTest++.h"
#include "constraint.h"
#include <iostream>

TEST(ConstraintsAreCreated) {
    auto c = continuity_constraint(1, 2);
    Constraint correct = {
        {DOFWeight{1, 1.0f}, DOFWeight{2, -1.0f}},
        0.0
    };
    CHECK(c.dof_constraints == correct.dof_constraints);
    CHECK(c.rhs_value == correct.rhs_value);

    auto c2 = offset_constraint(3, 4, 5.0f);
    correct = {
        {DOFWeight{3, 1.0f}, DOFWeight{4, -1.0f}},
        -5.0
    };
    CHECK(c2.dof_constraints == correct.dof_constraints);
    CHECK(c2.rhs_value == correct.rhs_value);
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
    auto in = cm.get_unconstrained({2.0, 4.0, 4.0, 4.0});
    auto res = cm.get_all(in, 4);
    double res_exact[4] = {2.0, 4.0, 4.0, 4.0};
    CHECK_ARRAY_CLOSE(res, res_exact, 4, 1e-13);
}

int main(int, char const *[])
{
    return UnitTest::RunAllTests();
}

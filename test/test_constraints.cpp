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

TEST(ConstraintMatrixApply) {
    auto c0 = boundary_condition(1, 4.0);
    auto c1 = continuity_constraint(1, 2);
    auto c2 = continuity_constraint(2, 3);
    auto cm = ConstraintMatrix::from_constraints({c0, c1, c2});
    auto res = cm.apply({2.0, 0.0, 0.0, 0.0});
    double res2[4] = {2.0, 4.0, 4.0, 4.0};
    CHECK_ARRAY_CLOSE(res, res2, 4, 1e-13);
}

// TEST(AddWithConstraintsSimple) {
//     arma::Mat<double> mymat = arma::zeros<mat>(4,4);
//     arma::Col<double> myrhs = arma::zeros<vec>(4);
//     ConstraintMatrix cm;
//     add_constraint(cm, continuity_constraint(2, 3));
//     add_mat_with_constraints(MatrixEntry({1, 1, 2}), mymat, myrhs, cm);
//     add_mat_with_constraints(MatrixEntry({2, 1, 2}), mymat, myrhs, cm);
//     CHECK_CLOSE(mymat(1,1), 2.0, 1e-6);
//     CHECK_CLOSE(mymat(3,1), 2.0, 1e-6);
// 
//     add_mat_with_constraints(MatrixEntry({2, 2, 3.7}), mymat, myrhs, cm);
//     CHECK_CLOSE(mymat(3,3), 3.7, 1e-6);
// 
//     add_rhs_with_constraints(DOFWeight({3, 2.0}), myrhs, cm);
//     CHECK_CLOSE(myrhs(3), 2.0, 1e-6);
// 
//     add_rhs_with_constraints(DOFWeight({2, 2.0}), myrhs, cm);
//     CHECK_CLOSE(myrhs(3), 4.0, 1e-6);
// }
// 
// TEST(AddWithConstraintsComplex) {
//     arma::Mat<double> mymat = arma::zeros<mat>(4,4);
//     arma::Col<double> myrhs = arma::zeros<vec>(4);
// 
//     ConstraintMatrix cm;
//     add_constraint(cm, continuity_constraint(2, 3));
// }
// 
// TEST(AddWithInhomogenousConstraints) {
//     arma::Mat<double> mymat = arma::zeros<mat>(4,4);
//     arma::Col<double> myrhs = arma::zeros<vec>(4);
// 
//     ConstraintMatrix cm;
//     add_constraint(cm, offset_constraint(2, 3, 0.5));
// 
//     add_mat_with_constraints(MatrixEntry({1, 2, 2.0}), mymat, myrhs, cm);
// 
//     CHECK_CLOSE(mymat(1,3), 2.0, 1e-6);
//     CHECK_CLOSE(myrhs(1), -1.0, 1e-6);
// }
// 
int main(int, char const *[])
{
    return UnitTest::RunAllTests();
}

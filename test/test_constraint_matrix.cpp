#include "UnitTest++.h"
#include "constraint_matrix.h"

using namespace tbem;

ConstraintMatrix two_bcs_constraint_map() {
    ConstraintEQ eqtn0{{LinearTerm{3,1}}, 4.0};
    ConstraintEQ eqtn1{{LinearTerm{1,1}}, 2.0};
    ConstraintMatrix constraint_set;
    constraint_set.insert(std::make_pair(1, isolate_term_on_lhs(eqtn1, 0)));
    constraint_set.insert(std::make_pair(3, isolate_term_on_lhs(eqtn0, 0)));
    return constraint_set;
}

TEST(IsConstrained) {
    auto constraint_set = two_bcs_constraint_map();
    CHECK(is_constrained(constraint_set, 0) == false);
    CHECK(is_constrained(constraint_set, 1) == true);
    CHECK(is_constrained(constraint_set, 2) == false);
    CHECK(is_constrained(constraint_set, 3) == true);
}

TEST(MakeLowerTriangular) {
    auto constraint_set = two_bcs_constraint_map();    
    //4x_2 - x_3 = 0 combined with x_3 = 4.0 -->
    //4x_2 - 4.0 = 0 -->
    //4x_2 = 4.0
    ConstraintEQ in{{LinearTerm{2,4}, LinearTerm{3,-1}}, 0.0};
    auto c_lower_tri = make_lower_triangular(in, constraint_set);
    CHECK_EQUAL(c_lower_tri.constrained_dof, 2);
    CHECK_EQUAL(c_lower_tri.terms.size(), 0);
    CHECK_EQUAL(c_lower_tri.rhs, 1);
}

void check_distribute_vector(const ConstraintMatrix& cm, 
                   const std::vector<double>& condensed,
                   const std::vector<double>& correct) {
    size_t n = correct.size();
    auto all_vals = distribute_vector(cm, condensed, n);
    CHECK_ARRAY_EQUAL(&all_vals[0], &correct[0], n);
}

void check_expands_to_all_ones(const ConstraintMatrix& cm, int n) {
    check_distribute_vector(cm, {1.0}, std::vector<double>(n, 1.0));
}

TEST(AEqualsB) {
    auto c0 = continuity_constraint(0, 1);
    auto cm = from_constraints({c0});
    check_expands_to_all_ones(cm, 2);
}

TEST(AEqualsBEqualsC) {
    auto cm = from_constraints({
        continuity_constraint(0, 1),
        continuity_constraint(1, 2)
    });
    CHECK_EQUAL(cm.size(), 2);
    check_expands_to_all_ones(cm, 3);
}

TEST(AEqualsBEqualsCEqualsD) {
    auto cm = from_constraints({
        continuity_constraint(0, 1),
        continuity_constraint(2, 1),
        continuity_constraint(2, 3)
    });
    check_expands_to_all_ones(cm, 4);
}

TEST(AEqualsBPlusC) {
    auto cm = from_constraints({
        ConstraintEQ{{LinearTerm{0, 1}, LinearTerm{1, -1}, LinearTerm{2, -1}}, 0}
    });
    check_distribute_vector(cm, {1.0, 1.0}, {1.0, 1.0, 0.0});
}

TEST(AEqualsBPlusCAndCEqualsD) {
    auto cm = from_constraints({
        ConstraintEQ{{LinearTerm{0, 1}, LinearTerm{1, -1}, LinearTerm{2, -1}}, 0},
        continuity_constraint(2, 3)
    });
    check_distribute_vector(cm, {1.0, 0.5}, {1.0, 0.5, 0.5, 0.5});
}

TEST(ZeroWeightConstraint) {
    auto cm = from_constraints({
        ConstraintEQ{{{0, 1.0}, {1, 0.0}}, 0.0},
    });
    CHECK_EQUAL(cm.at(0).terms.size(), 0);
}

TEST(EmptyConstraint) {
    auto cm = from_constraints({
        ConstraintEQ{{{1, 0.0}}, 0.0},
    });
    CHECK_EQUAL(cm.size(), 0);
}

TEST(CircularConstraints) {
    auto cm = from_constraints({
        continuity_constraint(0, 1),
        continuity_constraint(1, 0)
    });
    CHECK_EQUAL(cm.size(), 1);
}

TEST(CondenseUnconstrained) {
    auto cm = from_constraints({
        continuity_constraint(0, 2)
    });
    auto result = condense_vector(cm, {0, -1, 0});
    std::vector<double> correct{0, -1};
    CHECK_EQUAL(result.size(), 2);
    CHECK_ARRAY_EQUAL(result, correct, 2);
}

TEST(CondenseEmpty) {
    auto cm = from_constraints({
        boundary_condition(0, 4.0)
    });
    auto result = condense_vector(cm, {2.0});
    CHECK_EQUAL(result.size(), 0);
}

ConstraintMatrix bcs1_and_continuity0234(int bc_dof, double bc_val) {
    return from_constraints({
        boundary_condition(bc_dof, bc_val),
        continuity_constraint(0, 2),
        continuity_constraint(2, 3),
        continuity_constraint(3, 4)
    });
}

TEST(CondenseRecurse) {
    auto cm = bcs1_and_continuity0234(1, 4.0);
    auto result = condense_vector(cm, {0.0, 0.0, 0.0, 0.0, 4.0});
    CHECK_EQUAL(result.size(), 1);
    CHECK_EQUAL(result[0], 4.0);
}

TEST(CondenseThenDistribute) {
    auto cm = bcs1_and_continuity0234(1, 4.0);
    auto in = condense_vector(cm, std::vector<double>{2.0, 4.0, 4.0, 4.0, 4.0});
    CHECK_EQUAL(in[0], 14.0);
    auto res = distribute_vector(cm, in, 5);
    double res_exact[5] = {in[0], 4.0, in[0], in[0], in[0]};
    CHECK_ARRAY_CLOSE(res, res_exact, 5, 1e-13);
}

TEST(CondenseWithFullyDeterminedSubset) {
    auto cm = bcs1_and_continuity0234(2, 2.0);
    auto in = condense_vector(cm, std::vector<double>{2.0, 4.0, 4.0, 4.0, 4.0});
    CHECK_EQUAL(in[0], 4.0);
}

TEST(CondenseMatrixContinuity) {
    auto cm = from_constraints({
        continuity_constraint(0, 2),
        continuity_constraint(1, 2)
    });
    std::vector<double> matrix{
        {1,0,0  ,  0,1,0  ,  0,0,1}
    };
    auto result = condense_matrix(cm, cm, make_operator(3, 3, matrix));
    CHECK_EQUAL(result.n_rows * result.n_cols, 1);
    CHECK_EQUAL(result[0], 3);
}

TEST(CondenseMatrixContinuityPartial) {
    auto cm = from_constraints({
        continuity_constraint(1, 2)
    });
    std::vector<double> matrix{
        {1,0,0  ,  0,1,0  ,  0,0,1}
    };
    auto result = condense_matrix(cm, cm, make_operator(3, 3, matrix));
    CHECK_EQUAL(result.n_rows * result.n_cols, 4);
    std::vector<double> exact{1, 0, 0, 2};
    CHECK_ARRAY_EQUAL(result.data(), exact, 4);
}

TEST(CondenseMatrixBoundaryCondition) {
    auto cm = from_constraints({
        boundary_condition(1, 4.0)
    });
    std::vector<double> matrix{
        {1,0,  0,1}
    };
    auto result = condense_matrix(cm, cm, make_operator(2, 2, matrix));
    CHECK_EQUAL(result.n_rows * result.n_cols, 1);
    CHECK_EQUAL(result[0], 1.0);
}

TEST(CondenseNonSquareMatrixContinuity) {
    auto row_cm = from_constraints({
        continuity_constraint(0, 1),
    });
    auto col_cm = from_constraints({
        continuity_constraint(0, 2),
        continuity_constraint(1, 2)
    });
    std::vector<double> matrix{
        {1,0,0  ,  0,1,0}
    };
    auto result = condense_matrix(row_cm, col_cm, make_operator(2, 3, matrix));
    CHECK_EQUAL(result.n_rows * result.n_cols, 1);
    CHECK_EQUAL(result[0], 2);
}

BlockOperator condense_block_operator(const std::vector<ConstraintMatrix>& row_cms,
    const std::vector<ConstraintMatrix>& col_cms, const BlockOperator& op);

int main(int, char const *[])
{
    return UnitTest::RunAllTests();
}

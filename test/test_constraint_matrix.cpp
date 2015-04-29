#include "catch.hpp"
#include "constraint_matrix.h"
#include "vectorx.h"

using namespace tbem;

ConstraintMatrix two_bcs_constraint_map() {
    ConstraintEQ eqtn0{{LinearTerm{3,1}}, 4.0};
    ConstraintEQ eqtn1{{LinearTerm{1,1}}, 2.0};
    ConstraintMatrix constraint_set;
    constraint_set.insert(std::make_pair(1, isolate_term_on_lhs(eqtn1, 0)));
    constraint_set.insert(std::make_pair(3, isolate_term_on_lhs(eqtn0, 0)));
    return constraint_set;
}

TEST_CASE("IsConstrained", "[constraint_matrix]") {
    auto constraint_set = two_bcs_constraint_map();
    REQUIRE(is_constrained(constraint_set, 0) == false);
    REQUIRE(is_constrained(constraint_set, 1) == true);
    REQUIRE(is_constrained(constraint_set, 2) == false);
    REQUIRE(is_constrained(constraint_set, 3) == true);
}

TEST_CASE("MakeLowerTriangular", "[constraint_matrix]") {
    auto constraint_set = two_bcs_constraint_map();    
    //4x_2 - x_3 = 0 combined with x_3 = 4.0 -->
    //4x_2 - 4.0 = 0 -->
    //4x_2 = 4.0
    ConstraintEQ in{{LinearTerm{2,4}, LinearTerm{3,-1}}, 0.0};
    auto c_lower_tri = make_lower_triangular(in, constraint_set);
    REQUIRE(c_lower_tri.constrained_dof == 2);
    REQUIRE(c_lower_tri.terms.size() == 0);
    REQUIRE(c_lower_tri.rhs == 1);
}

void check_distribute_vector(const ConstraintMatrix& cm, 
                   const std::vector<double>& condensed,
                   const std::vector<double>& correct) {
    size_t n = correct.size();
    auto all_vals = distribute_vector(cm, condensed, n);
    REQUIRE_ARRAY_EQUAL(&all_vals[0], &correct[0], n);
}

void check_expands_to_all_ones(const ConstraintMatrix& cm, int n) {
    check_distribute_vector(cm, {1.0}, std::vector<double>(n, 1.0));
}

TEST_CASE("AEqualsB", "[constraint_matrix]") {
    auto c0 = continuity_constraint(0, 1);
    auto cm = from_constraints({c0});
    check_expands_to_all_ones(cm, 2);
}

TEST_CASE("AEqualsBEqualsC", "[constraint_matrix]") {
    auto cm = from_constraints({
        continuity_constraint(0, 1),
        continuity_constraint(1, 2)
    });
    REQUIRE(cm.size() == 2);
    check_expands_to_all_ones(cm, 3);
}

TEST_CASE("AEqualsBEqualsCEqualsD", "[constraint_matrix]") {
    auto cm = from_constraints({
        continuity_constraint(0, 1),
        continuity_constraint(2, 1),
        continuity_constraint(2, 3)
    });
    check_expands_to_all_ones(cm, 4);
}

TEST_CASE("AEqualsBPlusC", "[constraint_matrix]") {
    auto cm = from_constraints({
        ConstraintEQ{{LinearTerm{0, 1}, LinearTerm{1, -1}, LinearTerm{2, -1}}, 0}
    });
    check_distribute_vector(cm, {1.0, 1.0}, {1.0, 1.0, 0.0});
}

TEST_CASE("AEqualsBPlusCAndCEqualsD", "[constraint_matrix]") {
    auto cm = from_constraints({
        ConstraintEQ{{LinearTerm{0, 1}, LinearTerm{1, -1}, LinearTerm{2, -1}}, 0},
        continuity_constraint(2, 3)
    });
    check_distribute_vector(cm, {1.0, 0.5}, {1.0, 0.5, 0.5, 0.5});
}

TEST_CASE("ZeroWeightConstraint", "[constraint_matrix]") {
    auto cm = from_constraints({
        ConstraintEQ{{{0, 1.0}, {1, 0.0}}, 0.0},
    });
    REQUIRE(cm.at(0).terms.size() == 0);
}

TEST_CASE("EmptyConstraint", "[constraint_matrix]") {
    auto cm = from_constraints({
        ConstraintEQ{{{1, 0.0}}, 0.0},
    });
    REQUIRE(cm.size() == 0);
}

TEST_CASE("CircularConstraints", "[constraint_matrix]") {
    auto cm = from_constraints({
        continuity_constraint(0, 1),
        continuity_constraint(1, 0)
    });
    REQUIRE(cm.size() == 1);
}

TEST_CASE("CondenseUnconstrained", "[constraint_matrix]") {
    auto cm = from_constraints({
        continuity_constraint(0, 2)
    });
    auto result = condense_vector(cm, {0, -1, 0});
    std::vector<double> correct{0, -1};
    REQUIRE(result.size() == 2);
    REQUIRE_ARRAY_EQUAL(result, correct, 2);
}

TEST_CASE("CondenseEmpty", "[constraint_matrix]") {
    auto cm = from_constraints({
        boundary_condition(0, 4.0)
    });
    auto result = condense_vector(cm, {2.0});
    REQUIRE(result.size() == 0);
}

ConstraintMatrix bcs1_and_continuity0234(int bc_dof, double bc_val) {
    return from_constraints({
        boundary_condition(bc_dof, bc_val),
        continuity_constraint(0, 2),
        continuity_constraint(2, 3),
        continuity_constraint(3, 4)
    });
}

TEST_CASE("CondenseRecurse", "[constraint_matrix]") {
    auto cm = bcs1_and_continuity0234(1, 4.0);
    auto result = condense_vector(cm, {0.0, 0.0, 0.0, 0.0, 4.0});
    REQUIRE(result.size() == 1);
    REQUIRE(result[0] == 4.0);
}

TEST_CASE("CondenseThenDistribute", "[constraint_matrix]") {
    auto cm = bcs1_and_continuity0234(1, 4.0);
    auto in = condense_vector(cm, std::vector<double>{2.0, 4.0, 4.0, 4.0, 4.0});
    REQUIRE(in[0] == 14.0);
    auto res = distribute_vector(cm, in, 5);
    double res_exact[5] = {in[0], 4.0, in[0], in[0], in[0]};
    REQUIRE_ARRAY_CLOSE(res, res_exact, 5, 1e-13);
}

TEST_CASE("CondenseWithFullyDeterminedSubset", "[constraint_matrix]") {
    auto cm = bcs1_and_continuity0234(2, 2.0);
    auto in = condense_vector(cm, std::vector<double>{2.0, 4.0, 4.0, 4.0, 4.0});
    REQUIRE(in[0] == 4.0);
}

TEST_CASE("CondenseMatrixContinuity", "[constraint_matrix]") {
    auto cm = from_constraints({
        continuity_constraint(0, 2),
        continuity_constraint(1, 2)
    });
    std::vector<double> matrix{
        {1,0,0  ,  0,1,0  ,  0,0,1}
    };
    auto result = condense_matrix(cm, cm, DenseOperator(3, 3, matrix));
    REQUIRE(result.n_elements() == 1);
    REQUIRE(result[0] == 3);
}

TEST_CASE("CondenseMatrixContinuityPartial", "[constraint_matrix]") {
    auto cm = from_constraints({
        continuity_constraint(1, 2)
    });
    std::vector<double> matrix{
        {1,0,0  ,  0,1,0  ,  0,0,1}
    };
    auto result = condense_matrix(cm, cm, DenseOperator(3, 3, matrix));
    REQUIRE(result.n_elements() == 4);
    std::vector<double> exact{1, 0, 0, 2};
    REQUIRE_ARRAY_EQUAL(result.data(), exact, 4);
}

TEST_CASE("CondenseMatrixBoundaryCondition", "[constraint_matrix]") {
    auto cm = from_constraints({
        boundary_condition(1, 4.0)
    });
    std::vector<double> matrix{
        {1,0,  0,1}
    };
    auto result = condense_matrix(cm, cm, DenseOperator(2, 2, matrix));
    REQUIRE(result.n_elements() == 1);
    REQUIRE(result[0] == 1.0);
}

TEST_CASE("CondenseNonSquareMatrixContinuity", "[constraint_matrix]") {
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
    auto result = condense_matrix(row_cm, col_cm, DenseOperator(2, 3, matrix));
    REQUIRE(result.n_elements() == 1);
    REQUIRE(result[0] == 2);
}

BlockDenseOperator condense_block_operator(const std::vector<ConstraintMatrix>& row_cms,
    const std::vector<ConstraintMatrix>& col_cms, const BlockDenseOperator& op);

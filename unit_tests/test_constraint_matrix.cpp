#include "catch.hpp"
#include "constraint_matrix.h"

using namespace tbem;

ConstraintMatrix two_bcs_constraint_map() 
{
    ConstraintEQ eqtn0{{LinearTerm{3,1}}, 4.0};
    ConstraintEQ eqtn1{{LinearTerm{1,1}}, 2.0};
    ConstraintMatrixData constraint_set;
    constraint_set.insert(std::make_pair(1, isolate_term_on_lhs(eqtn1, 0)));
    constraint_set.insert(std::make_pair(3, isolate_term_on_lhs(eqtn0, 0)));
    return {constraint_set};
}

TEST_CASE("IsConstrained", "[constraint_matrix]") 
{
    auto constraint_set = two_bcs_constraint_map();
    REQUIRE(is_constrained(constraint_set.map, 0) == false);
    REQUIRE(is_constrained(constraint_set.map, 1) == true);
    REQUIRE(is_constrained(constraint_set.map, 2) == false);
    REQUIRE(is_constrained(constraint_set.map, 3) == true);
}

TEST_CASE("MakeLowerTriangular", "[constraint_matrix]") 
{
    auto constraint_set = two_bcs_constraint_map();    
    //4x_2 - x_3 = 0 combined with x_3 = 4.0 -->
    //4x_2 - 4.0 = 0 -->
    //4x_2 = 4.0
    ConstraintEQ in{{LinearTerm{2,4}, LinearTerm{3,-1}}, 0.0};
    auto c_lower_tri = make_lower_triangular(in, constraint_set.map);
    REQUIRE(c_lower_tri.terms.size() == 1);
    REQUIRE(c_lower_tri.terms[0].dof == 2);
    REQUIRE(c_lower_tri.terms[0].weight == 4.0);
    REQUIRE(c_lower_tri.rhs == 4.0);
}

void check_distribute_vector(const ConstraintMatrix& cm, 
                   const std::vector<double>& condensed,
                   const std::vector<double>& correct) 
{
    size_t n = correct.size();
    auto all_vals = distribute_vector(cm, condensed, n);
    REQUIRE_ARRAY_EQUAL(&all_vals[0], &correct[0], n);
}

void check_expands_to_all_ones(const ConstraintMatrix& cm, int n) 
{
    check_distribute_vector(cm, {1.0}, std::vector<double>(n, 1.0));
}

TEST_CASE("AEqualsB", "[constraint_matrix]") 
{
    auto c0 = continuity_constraint(0, 1);
    auto cm = from_constraints({c0});
    check_expands_to_all_ones(cm, 2);
}

TEST_CASE("AEqualsBEqualsC", "[constraint_matrix]") 
{
    auto cm = from_constraints({
        continuity_constraint(0, 1),
        continuity_constraint(1, 2)
    });
    REQUIRE(cm.size() == 2);
    check_expands_to_all_ones(cm, 3);
}

TEST_CASE("AEqualsBEqualsCEqualsD", "[constraint_matrix]") 
{
    auto cm = from_constraints({
        continuity_constraint(0, 1),
        continuity_constraint(2, 1),
        continuity_constraint(2, 3)
    });
    check_expands_to_all_ones(cm, 4);
}

TEST_CASE("AEqualsBPlusC", "[constraint_matrix]") 
{
    auto cm = from_constraints({
        ConstraintEQ{{LinearTerm{0, 1}, LinearTerm{1, -1}, LinearTerm{2, -1}}, 0}
    });
    check_distribute_vector(cm, {1.0, 1.0}, {1.0, 1.0, 0.0});
}

TEST_CASE("AEqualsBPlusCAndCEqualsD", "[constraint_matrix]") 
{
    auto cm = from_constraints({
        ConstraintEQ{{LinearTerm{0, 1}, LinearTerm{1, -1}, LinearTerm{2, -1}}, 0},
        continuity_constraint(2, 3)
    });
    check_distribute_vector(cm, {1.0, 0.5}, {1.0, 0.5, 0.5, 0.5});
}

TEST_CASE("ZeroWeightConstraint", "[constraint_matrix]") 
{
    auto cm = from_constraints({
        ConstraintEQ{{{0, 1.0}, {1, 0.0}}, 0.0},
    });
    REQUIRE(cm.map.at(0).terms.size() == 0);
}

TEST_CASE("EmptyConstraint", "[constraint_matrix]") 
{
    auto cm = from_constraints({
        ConstraintEQ{{{1, 0.0}}, 0.0},
    });
    REQUIRE(cm.size() == 0);
}

TEST_CASE("CircularConstraints", "[constraint_matrix]") 
{
    auto cm = from_constraints({
        continuity_constraint(0, 1),
        continuity_constraint(1, 0)
    });
    REQUIRE(cm.size() == 1);
}

TEST_CASE("CondenseUnconstrained", "[constraint_matrix]") 
{
    auto cm = from_constraints({
        continuity_constraint(0, 2)
    });
    auto result = condense_vector(cm, {0, -1, 0});
    std::vector<double> correct{0, -1};
    REQUIRE(result.size() == 2);
    REQUIRE_ARRAY_EQUAL(result, correct, 2);
}

TEST_CASE("CondenseEmpty", "[constraint_matrix]") 
{
    auto cm = from_constraints({
        boundary_condition(0, 4.0)
    });
    auto result = condense_vector(cm, {2.0});
    REQUIRE(result.size() == 0);
}

ConstraintMatrix bcs1_and_continuity0234(int bc_dof, double bc_val) 
{
    return from_constraints({
        boundary_condition(bc_dof, bc_val),
        continuity_constraint(0, 2),
        continuity_constraint(2, 3),
        continuity_constraint(3, 4)
    });
}

TEST_CASE("CondenseRecurse", "[constraint_matrix]") 
{
    auto cm = bcs1_and_continuity0234(1, 4.0);
    auto result = condense_vector(cm, {0.0, 0.0, 0.0, 0.0, 4.0});
    REQUIRE(result.size() == 1);
    REQUIRE(result[0] == 4.0);
}

TEST_CASE("CondenseThenDistribute", "[constraint_matrix]") 
{
    auto cm = bcs1_and_continuity0234(1, 4.0);
    auto in = condense_vector(cm, std::vector<double>{2.0, 4.0, 4.0, 4.0, 4.0});
    REQUIRE(in[0] == 14.0);
    auto res = distribute_vector(cm, in, 5);
    double res_exact[5] = {in[0], 4.0, in[0], in[0], in[0]};
    REQUIRE_ARRAY_CLOSE(res, res_exact, 5, 1e-13);
}

TEST_CASE("CondenseWithFullyDeterminedSubset", "[constraint_matrix]") 
{
    auto cm = bcs1_and_continuity0234(2, 2.0);
    auto in = condense_vector(cm, std::vector<double>{2.0, 4.0, 4.0, 4.0, 4.0});
    REQUIRE(in[0] == 4.0);
}

TEST_CASE("CondenseMatrixContinuity", "[constraint_matrix]") 
{
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

TEST_CASE("CondenseMatrixContinuityPartial", "[constraint_matrix]") 
{
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

TEST_CASE("CondenseMatrixBoundaryCondition", "[constraint_matrix]") 
{
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

TEST_CASE("CondenseNonSquareMatrixContinuity", "[constraint_matrix]") 
{
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

TEST_CASE("Sparse matrix continuity", "[constraint_matrix]") 
{
    auto cm = from_constraints({
        continuity_constraint(0, 2),
        continuity_constraint(1, 2)
    });
    auto matrix = SparseOperator::csr_from_coo(
        3, 3, {{0, 0, 1}, {1, 1, 1}, {2, 2, 1}}
    );
    auto result = condense_matrix(cm, cm, matrix);
    REQUIRE(result.nnz() == 1);
    REQUIRE(result.row_ptrs.size() == 2);
    REQUIRE(result.values[0] == 3);
    REQUIRE(result.column_indices[0] == 0);
}

TEST_CASE("Sparse matrix larger matrix", "[constraint_matrix]") 
{
    auto cm = from_constraints({
        continuity_constraint(1, 2)
    });
    auto matrix = SparseOperator::csr_from_coo(
        5, 5,
        {{0, 0, 1}, {0, 2, 4}, {1, 1, 1}, {2, 2, 1}, {2, 0, 1}, 
        {4, 0, -1}, {4, 3, 1}, {3, 4, 1}}
    );
    auto result = condense_matrix(cm, cm, matrix);
    auto correct = SparseOperator::csr_from_coo(
        4, 4,
        {{0, 0, 1}, {0, 1, 4}, {1, 0, 1}, {1, 1, 2}, {3, 0, -1},
        {3, 2, 1}, {2, 3, 1}}
    );
    auto nnz = correct.nnz();
    REQUIRE_ARRAY_EQUAL(result.values, correct.values, nnz);
    REQUIRE_ARRAY_EQUAL(result.column_indices, correct.column_indices, nnz);
    REQUIRE_ARRAY_EQUAL(result.row_ptrs, correct.row_ptrs, correct.n_rows() + 1);
}

TEST_CASE("homogenize constraint matrix", "[constraint_matrix]")
{
    auto cm = bcs1_and_continuity0234(1, 2.0);
    auto homogenized = homogenize_constraints(cm);
    for (auto it = homogenized.map.begin(); it != homogenized.map.end(); ++it) {
        REQUIRE(it->second.rhs == 0.0);
    }
}

TEST_CASE("average field value constraint", "[constraint_matrix]") 
{
    size_t n = 8;
    std::vector<LinearTerm> terms;
    for (size_t i = 0; i < n; i++) {
        terms.push_back(LinearTerm{i, 1.0});
    }
    ConstraintEQ avg_field_constraints{terms, 0.0};
    auto cm = from_constraints({
        continuity_constraint(0, 1),
        continuity_constraint(2, 3),
        continuity_constraint(4, 5),
        continuity_constraint(6, 7),
        avg_field_constraints
    });
    auto half_minus_one = n / 2 - 1;
    auto result = distribute_vector(cm, std::vector<double>(half_minus_one, 1.0), n);
    double result_sum = 0.0;
    for (size_t i = 0; i < n; i++) {
        result_sum += result[i];
    }
    REQUIRE(result_sum == 0.0);
    REQUIRE(-result[n - 1] == half_minus_one);
    REQUIRE(-result[n - 2] == half_minus_one);
}


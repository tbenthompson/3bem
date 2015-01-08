#include "UnitTest++.h"
#include "constraint.h"
#include <iostream>
#include "shared.h"
#include "mesh_gen.h"

using namespace tbem;

TEST(CreateLinearTerm) {
    LinearTerm dw{1, 2.0};
    CHECK_EQUAL(dw.dof, 1);
    CHECK_EQUAL(dw.weight, 2.0);
}

TEST(CreateConstraintEQ) {
    ConstraintEQ eqtn{{LinearTerm{0,3}, LinearTerm{1,-1}, LinearTerm{2,4}}, 13.7};
}

TEST(ContinuityConstraintsAreCreated) {
    auto c = continuity_constraint(1, 2);
    ConstraintEQ correct = {
        {LinearTerm{1, 1.0f}, LinearTerm{2, -1.0f}},
        0.0
    };
    CHECK(c.terms == correct.terms);
    CHECK(c.rhs == correct.rhs);
}

TEST(RearrangeConstraintEQ) {
    ConstraintEQ eqtn{{LinearTerm{0,3}, LinearTerm{1,-1}, LinearTerm{2,4}}, 13.7};
    auto rearranged = isolate_term_on_lhs(eqtn, 2);
    CHECK_EQUAL(rearranged.constrained_dof, 2);
    CHECK_EQUAL(rearranged.rhs, 13.7 / 4);
    CHECK_EQUAL(rearranged.terms[0], (LinearTerm{0, -3.0 / 4}));
    CHECK_EQUAL(rearranged.terms[1], (LinearTerm{1, 1.0 / 4}));
}
TEST(FindLastDOFIndex) {
    ConstraintEQ eqtn{{LinearTerm{0,3}, LinearTerm{2,4}, LinearTerm{1,-1}}, 13.7};
    int index = find_last_dof_index(eqtn);
    CHECK_EQUAL(index, 1);
}

TEST(FilterZeroTerms) {
    auto term0 = LinearTerm{1,1};
    auto term1 = LinearTerm{0,0};
    auto term2 = LinearTerm{2,2};
    ConstraintEQ c{{term0, term1, term2}, 2.0};
    auto result = filter_zero_terms(c);
    CHECK_EQUAL(result.terms.size(), 2);
    CHECK_EQUAL(result.terms[0], term0);
    CHECK_EQUAL(result.terms[1], term2);
}

void subs_test(const ConstraintEQ& subs_victim,
               const ConstraintEQ& subs_in,
               const ConstraintEQ& correct) {
    auto in_rearranged = isolate_term_on_lhs(subs_in, 0);
    auto result = substitute(subs_victim, 0, in_rearranged);
    CHECK_EQUAL(result.terms.size(), correct.terms.size());
    for (size_t i = 0; i < correct.terms.size(); i++) {
        CHECK_EQUAL(result.terms[i], correct.terms[i]);
    }
    CHECK_EQUAL(result.rhs, correct.rhs);
}

TEST(SubstituteRHS) {
    ConstraintEQ eqtn0{{LinearTerm{1,1}, LinearTerm{3,1}}, 4.0};
    ConstraintEQ eqtn1{{LinearTerm{1,1}}, 2.0};
    ConstraintEQ correct{{LinearTerm{3,1}}, 2.0};
    subs_test(eqtn0, eqtn1, correct);
}

TEST(SubstituteRHSNegation) {
    ConstraintEQ eqtn0{{LinearTerm{1,1}, LinearTerm{3,1}}, 4.0};
    ConstraintEQ eqtn1{{LinearTerm{1,-1}}, 2.0};
    ConstraintEQ correct{{LinearTerm{3,1}}, 6.0};
    subs_test(eqtn0, eqtn1, correct);
}

TEST(SubstituteWithTerms) {
    ConstraintEQ eqtn0{{LinearTerm{1,1}, LinearTerm{3,1}}, 4.0};
    ConstraintEQ eqtn1{{LinearTerm{1,1}, LinearTerm{2,-3}}, 0.0};
    ConstraintEQ correct{{LinearTerm{3,1}, LinearTerm{2,3}}, 4.0};
    subs_test(eqtn0, eqtn1, correct);
}

TEST(SubstituteWithTermsNegation) {
    ConstraintEQ eqtn0{{LinearTerm{1,-1}, LinearTerm{3,1}}, 4.0};
    ConstraintEQ eqtn1{{LinearTerm{1,1}, LinearTerm{2,-3}}, 0.0};
    ConstraintEQ correct{{LinearTerm{3,1}, LinearTerm{2,-3}}, 4.0};
    subs_test(eqtn0, eqtn1, correct);
}

TEST(SubstituteWithTermsAddToPreexistingTerm) {
    ConstraintEQ eqtn0{{LinearTerm{1,1}, LinearTerm{3,1}}, 4.0};
    ConstraintEQ eqtn1{{LinearTerm{1,1}, LinearTerm{3,-3}}, 0.0};
    ConstraintEQ correct{{LinearTerm{3,4}}, 4.0};
    subs_test(eqtn0, eqtn1, correct);
}

ConstraintMapT two_bcs_constraint_map() {
    ConstraintEQ eqtn0{{LinearTerm{3,1}}, 4.0};
    ConstraintEQ eqtn1{{LinearTerm{1,1}}, 2.0};
    ConstraintMapT constraint_set;
    constraint_set[1] = isolate_term_on_lhs(eqtn1, 0);
    constraint_set[3] = isolate_term_on_lhs(eqtn0, 0);
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

TEST(ConstraintEQOutput) {
    ConstraintEQ c{{{0, 1}, {1, 1}}, 0.0};
    std::stringstream output_buf;
    output_buf << c;
    CHECK_EQUAL(output_buf.str(),
                "ConstraintEQ[[(rhs, 0), (0, 1), (1, 1), ]]");
}

TEST(RearrangedConstraintEQOutput) {
    auto rearranged_c = isolate_term_on_lhs(ConstraintEQ{{{0, 1}, {1, 1}}, 0.0}, 0);
    std::stringstream output_buf;
    output_buf << rearranged_c;
    CHECK_EQUAL(output_buf.str(), 
                "RearrangedConstraintEQ[[(constrained_dof=0, 1), (rhs, 0), (1, -1), ]]");
}

TEST(ConstraintMatrixFromConstraints) {
    ConstraintEQ eqtn{{LinearTerm{3,1}}, 4.0};
    auto matrix = ConstraintMatrix::from_constraints({eqtn});
    CHECK_EQUAL(matrix.map.size(), 1);
    const auto rearranged_constraint = matrix.map.find(3);
    CHECK(rearranged_constraint != matrix.map.end());
    CHECK_EQUAL(rearranged_constraint->second.rhs, 4.0);
}

void check_get_all(const ConstraintMatrix& cm, 
                   const std::vector<double>& reduced,
                   const std::vector<double>& correct) {
    size_t n = correct.size();
    auto all_vals = cm.get_all(reduced, n);
    CHECK_ARRAY_EQUAL(&all_vals[0], &correct[0], n);
}

void check_expands_to_all_ones(const ConstraintMatrix& cm, int n) {
    check_get_all(cm, {1.0}, std::vector<double>(n, 1.0));
}

TEST(AEqualsB) {
    auto c0 = continuity_constraint(0, 1);
    auto cm = ConstraintMatrix::from_constraints({c0});
    check_expands_to_all_ones(cm, 2);
}

TEST(AEqualsBEqualsC) {
    auto cm = ConstraintMatrix::from_constraints({
        continuity_constraint(0, 1),
        continuity_constraint(1, 2)
    });
    CHECK_EQUAL(cm.map.size(), 2);
    check_expands_to_all_ones(cm, 3);
}

TEST(AEqualsBEqualsCEqualsD) {
    auto cm = ConstraintMatrix::from_constraints({
        continuity_constraint(0, 1),
        continuity_constraint(2, 1),
        continuity_constraint(2, 3)
    });
    check_expands_to_all_ones(cm, 4);
}

TEST(AEqualsBPlusC) {
    auto cm = ConstraintMatrix::from_constraints({
        ConstraintEQ{{LinearTerm{0, 1}, LinearTerm{1, -1}, LinearTerm{2, -1}}, 0}
    });
    check_get_all(cm, {1.0, 1.0}, {1.0, 1.0, 0.0});
}

TEST(AEqualsBPlusCAndCEqualsD) {
    auto cm = ConstraintMatrix::from_constraints({
        ConstraintEQ{{LinearTerm{0, 1}, LinearTerm{1, -1}, LinearTerm{2, -1}}, 0},
        continuity_constraint(2, 3)
    });
    check_get_all(cm, {1.0, 0.5}, {1.0, 0.5, 0.5, 0.5});
}

TEST(CircularConstraints) {
    auto cm = ConstraintMatrix::from_constraints({
        continuity_constraint(0, 1),
        continuity_constraint(1, 0)
    });
    CHECK_EQUAL(cm.map.size(), 1);
}

TEST(AddTermWithConstraintsUnconstrained) {
    auto cm = ConstraintMatrix::from_constraints({
        continuity_constraint(0, 2)
    });
    auto terms = cm.add_term_with_constraints(LinearTerm{1, -1.0});
    LinearTerm correct{1, -1.0};
    CHECK_EQUAL(terms.size(), 1);
    CHECK_EQUAL(terms[0], correct);
}

TEST(AddTermWithConstraintsEmpty) {
    auto cm = ConstraintMatrix::from_constraints({
        boundary_condition(1, 4.0)
    });
    auto terms = cm.add_term_with_constraints(LinearTerm{1, 2.0});
    CHECK(terms.empty());
}

TEST(AddTermWithConstraintsRecurse) {
    auto cm = ConstraintMatrix::from_constraints({
        boundary_condition(1, 4.0),
        continuity_constraint(0, 2),
        continuity_constraint(2, 3),
        continuity_constraint(3, 4)
    });
    auto terms = cm.add_term_with_constraints(LinearTerm{4, 4.0});
    LinearTerm correct{0, 4.0};
    CHECK_EQUAL(terms.size(), 1);
    CHECK_EQUAL(terms[0], correct);
}

TEST(ConstraintMatrixGetReducedThenGetAll) {
    auto cm = ConstraintMatrix::from_constraints({
        boundary_condition(1, 4.0),
        continuity_constraint(0, 2),
        continuity_constraint(2, 3),
        continuity_constraint(3, 4)
    });
    auto in = cm.get_reduced(std::vector<double>{2.0, 4.0, 4.0, 4.0, 4.0});
    CHECK_EQUAL(in[0], 14.0);
    auto res = cm.get_all(in, 5);
    double res_exact[5] = {in[0], 4.0, in[0], in[0], in[0]};
    CHECK_ARRAY_CLOSE(res, res_exact, 5, 1e-13);
}

TEST(ConstraintMatrixGetAllVec2) {
    auto c0 = boundary_condition(1, 4.0);
    auto c1 = continuity_constraint(1, 2);
    auto c2 = continuity_constraint(2, 3);
    auto cm = ConstraintMatrix::from_constraints({c0, c1, c2});
    auto in = cm.get_reduced(std::vector<Vec2<double>>{
        {2.0,3.0}, {4.0,4.0}, {4.0,1.5}, {4.0,-1.5}
    });
    auto res = cm.get_all(in, 4);
    Vec2<double> res_exact[4] = {{2.0,3.0}, {4.0,4.0}, {4.0,4.0}, {4.0,4.0}};
    CHECK_ARRAY_CLOSE((&res[0][0]), (&res_exact[0][0]), 8, 1e-13);
}

TEST(ConstraintMesh) {
    auto sphere = sphere_mesh({0, 0, 0}, 1).refine_repeatedly(2);
    auto constraints = mesh_continuity<3>(sphere);
    auto matrix = ConstraintMatrix::from_constraints(constraints);
    CHECK_EQUAL(3 * sphere.facets.size(), 384);
    CHECK_EQUAL(matrix.map.size(), 318);
    auto my_c = matrix.map.begin()->second.terms;
    CHECK_EQUAL(my_c[0].weight, 1);
}

TEST(GetReducedToCountTheNumberOfVerticesOnASphereApproximation) {
    auto sphere = sphere_mesh({0, 0, 0}, 1).refine_repeatedly(0);
    auto constraints = mesh_continuity<3>(sphere);
    auto matrix = ConstraintMatrix::from_constraints(constraints);
    std::vector<double> all(3 * sphere.facets.size(), 1.0);
    auto reduced = matrix.get_reduced(all);
    CHECK_EQUAL(reduced.size(), 6);
    CHECK_ARRAY_EQUAL(reduced, (std::vector<double>(6, 4.0)), 6);
}

int main(int, char const *[])
{
    return UnitTest::RunAllTests();
    // return RunOneTest("CircularConstraints");
}

#include "catch.hpp"
#include "constraint.h"

using namespace tbem;

TEST_CASE("ConstraintEquality", "[constraint]") 
{
    ConstraintEQ eqtn0{{LinearTerm{0,3}, LinearTerm{1,-1}, LinearTerm{2,4}}, 13.7};
    ConstraintEQ eqtn1{{LinearTerm{1,-1}}, 1.0};
    CHECK(eqtn0 == eqtn0);
    CHECK(eqtn1 == eqtn1);
}

TEST_CASE("RearrangeConstraintEQ", "[constraint]") 
{
    ConstraintEQ eqtn{{LinearTerm{0,3}, LinearTerm{1,-1}, LinearTerm{2,4}}, 13.7};
    auto rearranged = isolate_term_on_lhs(eqtn, 2);
    REQUIRE(rearranged.constrained_dof == 2);
    REQUIRE(rearranged.rhs == 13.7 / 4);
    REQUIRE(rearranged.terms[0] == (LinearTerm{0, -3.0 / 4}));
    REQUIRE(rearranged.terms[1] == (LinearTerm{1, 1.0 / 4}));
}

TEST_CASE("FindLastDOFIndex", "[constraint]") 
{
    ConstraintEQ eqtn{{LinearTerm{0,3}, LinearTerm{2,4}, LinearTerm{1,-1}}, 13.7};
    int index = find_last_dof_index(eqtn);
    REQUIRE(index == 1);
}

TEST_CASE("FilterZeroTerms", "[constraint]") 
{
    auto term0 = LinearTerm{1,1};
    auto term1 = LinearTerm{0,0};
    auto term2 = LinearTerm{2,2};
    ConstraintEQ c{{term0, term1, term2}, 2.0};
    auto result = filter_zero_terms(c);
    REQUIRE(result.terms.size() == 2);
    REQUIRE(result.terms[0] == term0);
    REQUIRE(result.terms[1] == term2);
}

TEST_CASE("combine terms", "[constraint]") 
{
    ConstraintEQ c{{LinearTerm{1, 1}, LinearTerm{2, 1}, LinearTerm{1, 1}}, 0.0};
    auto result = combine_terms(c);
    ConstraintEQ correct{{LinearTerm{1, 2}, LinearTerm{2, 1}}, 0.0};
    REQUIRE(result == correct);
}

void subs_test(const ConstraintEQ& subs_victim,
               const ConstraintEQ& subs_in,
               const ConstraintEQ& correct) 
{
    auto in_rearranged = isolate_term_on_lhs(subs_in, 0);
    auto result = substitute(subs_victim, 0, in_rearranged);
    REQUIRE(result.terms.size() == correct.terms.size());
    for (size_t i = 0; i < correct.terms.size(); i++) {
        REQUIRE(result.terms[i] == correct.terms[i]);
    }
    REQUIRE(result.rhs == correct.rhs);
}

TEST_CASE("SubstituteRHS", "[constraint]") 
{
    ConstraintEQ eqtn0{{LinearTerm{1,1}, LinearTerm{3,1}}, 4.0};
    ConstraintEQ eqtn1{{LinearTerm{1,1}}, 2.0};
    ConstraintEQ correct{{LinearTerm{3,1}}, 2.0};
    subs_test(eqtn0, eqtn1, correct);
}

TEST_CASE("SubstituteRHSNegation", "[constraint]") 
{
    ConstraintEQ eqtn0{{LinearTerm{1,1}, LinearTerm{3,1}}, 4.0};
    ConstraintEQ eqtn1{{LinearTerm{1,-1}}, 2.0};
    ConstraintEQ correct{{LinearTerm{3,1}}, 6.0};
    subs_test(eqtn0, eqtn1, correct);
}

TEST_CASE("SubstituteWithTerms", "[constraint]") 
{
    ConstraintEQ eqtn0{{LinearTerm{1,1}, LinearTerm{3,1}}, 4.0};
    ConstraintEQ eqtn1{{LinearTerm{1,1}, LinearTerm{2,-3}}, 0.0};
    ConstraintEQ correct{{LinearTerm{3,1}, LinearTerm{2,3}}, 4.0};
    subs_test(eqtn0, eqtn1, correct);
}

TEST_CASE("SubstituteWithTermsNegation", "[constraint]") 
{
    ConstraintEQ eqtn0{{LinearTerm{1,-1}, LinearTerm{3,1}}, 4.0};
    ConstraintEQ eqtn1{{LinearTerm{1,1}, LinearTerm{2,-3}}, 0.0};
    ConstraintEQ correct{{LinearTerm{3,1}, LinearTerm{2,-3}}, 4.0};
    subs_test(eqtn0, eqtn1, correct);
}

TEST_CASE("SubstituteWithTermsAddToPreexistingTerm", "[constraint]") 
{
    ConstraintEQ eqtn0{{LinearTerm{1,1}, LinearTerm{3,1}}, 4.0};
    ConstraintEQ eqtn1{{LinearTerm{1,1}, LinearTerm{3,-3}}, 0.0};
    ConstraintEQ correct{{LinearTerm{3,4}}, 4.0};
    subs_test(eqtn0, eqtn1, correct);
}

TEST_CASE("ShiftConstraints", "[constraint]") 
{
    ConstraintEQ eqtn0{{LinearTerm{1,1}, LinearTerm{3,1}}, 4.0};
    auto cs = shift_constraints({eqtn0}, 1);
    REQUIRE(cs[0].terms[0].dof == 2);
    REQUIRE(cs[0].terms[1].dof == 4);
}

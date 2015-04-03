#include "UnitTest++.h"
#include "constraint.h"

using namespace tbem;

TEST(ConstraintEquality) {
    ConstraintEQ eqtn0{{LinearTerm{0,3}, LinearTerm{1,-1}, LinearTerm{2,4}}, 13.7};
    ConstraintEQ eqtn1{{LinearTerm{1,-1}}, 1.0};
    CHECK(eqtn0 == eqtn0);
    CHECK(eqtn1 == eqtn1);
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

int main(int, char const *[])
{
    return UnitTest::RunAllTests();
}

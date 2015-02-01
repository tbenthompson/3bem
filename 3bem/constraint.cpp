#include "constraint.h"

namespace tbem {

RearrangedConstraintEQ isolate_term_on_lhs(const ConstraintEQ& c, 
                                                  size_t constrained_index) {
    assert(constrained_index < c.terms.size());
    const auto& constrained_term = c.terms[constrained_index]; 

    std::vector<LinearTerm> divided_negated_terms;
    for (const auto& t: c.terms) {
        if (t == constrained_term) {
            continue;
        }
        LinearTerm divided_negated{t.dof, -t.weight / constrained_term.weight};
        divided_negated_terms.push_back(divided_negated);
    }

    double divided_rhs = c.rhs / constrained_term.weight;

    return RearrangedConstraintEQ{
        constrained_term.dof,
        divided_negated_terms,
        divided_rhs
    };
}


size_t find_last_dof_index(const ConstraintEQ& c) {
    auto highest_dof_term = std::max_element(c.terms.begin(), c.terms.end(),
        [&] (const LinearTerm& a, const LinearTerm& b) {
            return a.dof < b.dof;
        });
    return std::distance(c.terms.begin(), highest_dof_term);
}

typedef std::vector<LinearTerm>::const_iterator LinearTermIterator;
LinearTermIterator find_term_with_dof(const std::vector<LinearTerm>& terms, size_t dof) {
    return std::find_if(terms.begin(), terms.end(),
        [&] (const LinearTerm& sub_lt) {
            return sub_lt.dof == dof;
        }
    );
}

template <typename T>
bool none_found(const typename T::const_iterator& it, const T& vec) {
    return it == std::end(vec);
}

template <typename T>
void remove(std::vector<T>& vec, size_t pos)
{
    typename std::vector<T>::iterator it = vec.begin();
    std::advance(it, pos);
    vec.erase(it);
}

ConstraintEQ substitute(const ConstraintEQ& c, size_t constrained_dof_index,
                        const RearrangedConstraintEQ& subs_in) 
{
    const auto& constrained_term = c.terms[constrained_dof_index];
    assert(constrained_term.dof == subs_in.constrained_dof);
    double multiplicative_factor = constrained_term.weight;

    std::vector<LinearTerm> out_terms;
    auto which_subs_terms_unused = range(subs_in.terms.size());
    for (size_t i = 0; i < c.terms.size(); i++) {
        if (i == constrained_dof_index) {
            continue;
        }
        const auto& term = c.terms[i];
        const auto& subs_term = find_term_with_dof(subs_in.terms, term.dof);

        if (none_found(subs_term, subs_in.terms)) {
            out_terms.push_back(term); 
        } else {
            double additive_weight = subs_term->weight * multiplicative_factor;
            out_terms.push_back(LinearTerm{term.dof, term.weight + additive_weight});
            size_t subs_term_index = std::distance(subs_in.terms.begin(), subs_term);
            assert(subs_term_index < subs_in.terms.size());
            remove(which_subs_terms_unused, subs_term_index);
        } 
    }

    for (size_t subs_term_index: which_subs_terms_unused) {
        const auto& subs_term = subs_in.terms[subs_term_index];
        double additive_weight = subs_term.weight * multiplicative_factor;
        out_terms.push_back(LinearTerm{subs_term.dof, additive_weight});
    }

    double out_rhs = c.rhs - multiplicative_factor * subs_in.rhs;
    return {out_terms, out_rhs};
}

ConstraintEQ filter_zero_terms(const ConstraintEQ& c, double eps) 
{
    std::vector<LinearTerm> out_terms;
    for(const auto& t: c.terms) {
        if (std::fabs(t.weight) > eps) {
            out_terms.push_back(t);
        }
    }
    return {out_terms, c.rhs};
}

} //END namespace tbem

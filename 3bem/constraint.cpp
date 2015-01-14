#include "constraint.h"

namespace tbem {

RearrangedConstraintEQ isolate_term_on_lhs(const ConstraintEQ& c, 
                                                  int constrained_index) {
    assert(constrained_index < (int)c.terms.size());
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

bool is_constrained(const ConstraintMapT& dof_constraint_map,
                           int dof) {
    const auto& it = dof_constraint_map.find(dof);
    if (it == dof_constraint_map.end()) {
        return false;
    }
    return true;
}

int find_last_dof_index(const ConstraintEQ& c) {
    auto highest_dof_term = std::max_element(c.terms.begin(), c.terms.end(),
        [&] (const LinearTerm& a, const LinearTerm& b) {
            return a.dof < b.dof;
        });
    return std::distance(c.terms.begin(), highest_dof_term);
}

typedef std::vector<LinearTerm>::const_iterator LinearTermIterator;
LinearTermIterator find_term_with_dof(const std::vector<LinearTerm>& terms, int dof) {
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

ConstraintEQ substitute(const ConstraintEQ& c, int constrained_dof_index,
                        const RearrangedConstraintEQ& subs_in) {

    const auto& constrained_term = c.terms[constrained_dof_index];
    assert(constrained_term.dof == subs_in.constrained_dof);
    double multiplicative_factor = constrained_term.weight;

    std::vector<LinearTerm> out_terms;
    auto which_subs_terms_unused = range(subs_in.terms.size());
    for (size_t i = 0; i < c.terms.size(); i++) {
        if ((int)i == constrained_dof_index) {
            continue;
        }
        const auto& term = c.terms[i];
        const auto& subs_term = find_term_with_dof(subs_in.terms, term.dof);

        if (none_found(subs_term, subs_in.terms)) {
            out_terms.push_back(term); 
        } else {
            double additive_weight = subs_term->weight * multiplicative_factor;
            out_terms.push_back(LinearTerm{term.dof, term.weight + additive_weight});
            int subs_term_index = std::distance(subs_in.terms.begin(), subs_term);
            assert(subs_term_index < (int)subs_in.terms.size());
            remove(which_subs_terms_unused, subs_term_index);
        } 
    }

    for (int subs_term_index: which_subs_terms_unused) {
        const auto& subs_term = subs_in.terms[subs_term_index];
        double additive_weight = subs_term.weight * multiplicative_factor;
        out_terms.push_back(LinearTerm{subs_term.dof, additive_weight});
    }

    double out_rhs = c.rhs - multiplicative_factor * subs_in.rhs;
    return {out_terms, out_rhs};
}

ConstraintEQ filter_zero_terms(const ConstraintEQ& c, double eps) {
    std::vector<LinearTerm> out_terms;
    for(const auto& t: c.terms) {
        if (std::fabs(t.weight) > eps) {
            out_terms.push_back(t);
        }
    }
    return {out_terms, c.rhs};
}

RearrangedConstraintEQ make_lower_triangular(const ConstraintEQ& c,
                                             const ConstraintMapT& map) {
    if (c.terms.size() == 0) {
        std::string msg = "Function: make_lower_triangular has found either an empty constraint or a cyclic set of constraints.";
        throw std::invalid_argument(msg);
    }

    int last_dof_index = find_last_dof_index(c);
    int last_dof = c.terms[last_dof_index].dof;
    if (is_constrained(map, last_dof)) {
        const auto& last_dof_constraint = map.find(last_dof)->second;
        ConstraintEQ c_subs = substitute(c, last_dof_index, last_dof_constraint);
        ConstraintEQ c_subs_filtered = filter_zero_terms(c_subs);
        return make_lower_triangular(c_subs_filtered, map);
    }
    return isolate_term_on_lhs(c, last_dof_index);
}

ConstraintMatrix ConstraintMatrix::from_constraints(
        const std::vector<ConstraintEQ>& constraints) {
    ConstraintMapT new_map;

    for (size_t i = 0; i < constraints.size(); i++) {
        const auto& c = constraints[i];
        try {
            auto lower_tri_constraint = make_lower_triangular(c, new_map);
            new_map[lower_tri_constraint.constrained_dof] =
                std::move(lower_tri_constraint);
        } catch (const std::invalid_argument& e) {
            continue;
        }
    }

    return ConstraintMatrix{new_map};
};

} //END namespace tbem

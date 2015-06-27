#include "constraint.h"
#include <iostream>
#include <algorithm>
#include <cmath>

namespace tbem {

bool LinearTerm::operator==(const LinearTerm& other) const 
{
    return dof == other.dof && weight == other.weight;
}

std::ostream& operator<<(std::ostream& os, const LinearTerm& lt) 
{
    os << "(" << lt.dof << ", " << lt.weight << ")";
    return os;
}
    
bool ConstraintEQ::operator==(const ConstraintEQ& other) const 
{
    bool equal = other.terms.size() == terms.size();
    if (!equal) {
        return false;
    }
    for (size_t i = 0; i < terms.size(); i++) {
        equal = equal && terms[i] == other.terms[i];
    }
    equal = equal && (rhs == other.rhs);
    return equal;
}

std::ostream& operator<<(std::ostream& os, const ConstraintEQ& c) 
{
    os << "ConstraintEQ[[";
    os << "(rhs, " << c.rhs << "), ";
    for (auto t: c.terms) {
        os << t << ", ";
    }
    os << "]]";
    return os;
}

std::ostream& operator<<(std::ostream& os, 
    const RearrangedConstraintEQ& c) 
{
    os << "RearrangedConstraintEQ[[";
    os << "(constrained_dof=" << c.constrained_dof << ", 1), ";
    os << "(rhs, " << c.rhs << "), ";
    for (auto t: c.terms) {
        os << t << ", ";
    }
    os << "]]";
    return os;
}

ConstraintEQ boundary_condition(size_t dof, double value) 
{
    return {{LinearTerm{dof, 1.0}}, value};
}

ConstraintEQ continuity_constraint(size_t dof1, size_t dof2) 
{
    return {
        {LinearTerm{dof1, 1.0}, LinearTerm{dof2, -1.0}},
        0.0
    };
}

RearrangedConstraintEQ isolate_term_on_lhs(const ConstraintEQ& c, 
    size_t constrained_index) 
{
    assert(constrained_index < c.terms.size());
    auto constrained_term = c.terms[constrained_index]; 

    std::vector<LinearTerm> divided_negated_terms;
    for (auto t: c.terms) {
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


size_t find_last_dof_index(const ConstraintEQ& c) 
{
    auto highest_dof_term = std::max_element(c.terms.begin(), c.terms.end(),
        [&] (const LinearTerm& a, const LinearTerm& b) {
            return a.dof < b.dof;
        });
    return std::distance(c.terms.begin(), highest_dof_term);
}

typedef std::vector<LinearTerm>::const_iterator LinearTermIterator;
LinearTermIterator find_term_with_dof(const std::vector<LinearTerm>& terms, size_t dof) 
{
    return std::find_if(terms.begin(), terms.end(),
        [&] (const LinearTerm& sub_lt) {
            return sub_lt.dof == dof;
        }
    );
}

template <typename T>
bool none_found(const typename T::const_iterator& it, const T& vec) 
{
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
    auto constrained_term = c.terms[constrained_dof_index];
    assert(constrained_term.dof == subs_in.constrained_dof);
    double multiplicative_factor = constrained_term.weight;

    std::vector<LinearTerm> out_terms;
    auto which_subs_terms_unused = range(subs_in.terms.size());
    for (size_t i = 0; i < c.terms.size(); i++) {
        if (i == constrained_dof_index) {
            continue;
        }
        auto term = c.terms[i];
        auto subs_term = find_term_with_dof(subs_in.terms, term.dof);

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
        auto subs_term = subs_in.terms[subs_term_index];
        double additive_weight = subs_term.weight * multiplicative_factor;
        out_terms.push_back(LinearTerm{subs_term.dof, additive_weight});
    }

    double out_rhs = c.rhs - multiplicative_factor * subs_in.rhs;
    return {out_terms, out_rhs};
}

ConstraintEQ filter_zero_terms(const ConstraintEQ& c, double eps) 
{
    std::vector<LinearTerm> out_terms;
    for(auto t: c.terms) {
        if (std::fabs(t.weight) > eps) {
            out_terms.push_back(t);
        }
    }
    return {out_terms, c.rhs};
}

std::vector<ConstraintEQ> 
shift_constraints(const std::vector<ConstraintEQ>& constraints, size_t shift_dof) {
    std::vector<ConstraintEQ> out;
    for (size_t i = 0; i < constraints.size(); i++) {
        auto n_terms = constraints[i].terms.size();
        std::vector<LinearTerm> out_terms; 
        for (size_t t_idx = 0; t_idx < n_terms; t_idx++) {
            out_terms.push_back({
                constraints[i].terms[t_idx].dof + shift_dof,
                constraints[i].terms[t_idx].weight
            });
        }
        out.push_back({out_terms, constraints[i].rhs});
    }

    return out;
}

} //END namespace tbem

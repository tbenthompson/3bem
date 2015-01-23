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

bool is_constrained(const ConstraintMatrix& dof_constraint_map,
                           size_t dof) {
    const auto& it = dof_constraint_map.find(dof);
    if (it == dof_constraint_map.end()) {
        return false;
    }
    return true;
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
                        const RearrangedConstraintEQ& subs_in) {

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
                                             const ConstraintMatrix& matrix) {
    if (c.terms.size() == 0) {
        std::string msg = "Function: make_lower_triangular has found either an empty constraint or a cyclic set of constraints.";
        throw std::invalid_argument(msg);
    }

    auto last_dof_index = find_last_dof_index(c);
    auto last_dof = c.terms[last_dof_index].dof;
    if (is_constrained(matrix, last_dof)) {
        const auto& last_dof_constraint = matrix.find(last_dof)->second;
        ConstraintEQ c_subs = substitute(c, last_dof_index, last_dof_constraint);
        ConstraintEQ c_subs_filtered = filter_zero_terms(c_subs);
        return make_lower_triangular(c_subs_filtered, matrix);
    }
    return isolate_term_on_lhs(c, last_dof_index);
}

ConstraintMatrix from_constraints(
        const std::vector<ConstraintEQ>& constraints) {
    ConstraintMatrix new_mat;

    for (size_t i = 0; i < constraints.size(); i++) {
        const auto& c = constraints[i];
        try {
            auto lower_tri_constraint = make_lower_triangular(c, new_mat);
            new_mat[lower_tri_constraint.constrained_dof] =
                std::move(lower_tri_constraint);
        } catch (const std::invalid_argument& e) {
            continue;
        }
    }

    return ConstraintMatrix{new_mat};
};

std::vector<double>
distribute_vector(const ConstraintMatrix& matrix, const std::vector<double>& in,
                  size_t total_dofs) 
{
    std::vector<double> out(total_dofs); 

    size_t next_reduced_dof = 0;

    for (size_t dof_index = 0; dof_index < total_dofs; dof_index++) {
        if(is_constrained(matrix, dof_index)) {
            continue;
        }
        out[dof_index] = in[next_reduced_dof];
        next_reduced_dof++;
    }

    for (size_t dof_index = 0; dof_index < total_dofs; dof_index++) {
        if(!is_constrained(matrix, dof_index)) {
            continue;
        }
        auto constraint = matrix.find(dof_index)->second;
        auto val = constraint.rhs;
        for (size_t j = 0; j < constraint.terms.size(); j++) {
            val += constraint.terms[j].weight * out[constraint.terms[j].dof];
        }
        out[dof_index] = val;
    }
    return out;
}

std::vector<LinearTerm>
add_term_with_constraints(const ConstraintMatrix& matrix, const LinearTerm& entry) {
    if (!is_constrained(matrix, entry.dof)) {
        return {entry};
    }

    const auto& constraint = matrix.find(entry.dof)->second;
    const auto& terms = constraint.terms;
    std::vector<LinearTerm> out_terms;
    for (size_t i = 0; i < terms.size(); i++) {
        double recurse_weight = terms[i].weight * entry.weight;
        LinearTerm new_entry{terms[i].dof, recurse_weight};
        const auto& terms = add_term_with_constraints(matrix, new_entry);
        for (const auto& t: terms) {
            out_terms.push_back(std::move(t));
        }
    }
    return out_terms;
}

std::vector<double> condense_vector(const ConstraintMatrix& matrix,
    const std::vector<double>& all) 
{
    std::vector<double> condensed_dofs(all.size(), 0.0);
    for (size_t dof_idx = 0; dof_idx < all.size(); dof_idx++) {
        LinearTerm term_to_add{dof_idx, all[dof_idx]};
        auto expanded_term = add_term_with_constraints(matrix, term_to_add);
        for (const auto& t: expanded_term) {
            condensed_dofs[t.dof] += t.weight;
        }
    }

    std::vector<double> out;
    for (size_t dof_idx = 0; dof_idx < all.size(); dof_idx++) {
        if (is_constrained(matrix, dof_idx)) {
            continue;
        }
        out.push_back(condensed_dofs[dof_idx]);
    }

    return out;
}

std::vector<MatrixEntry>
entries_from_terms(const std::vector<LinearTerm>& terms, size_t axis, 
        MatrixEntry original_entry)
{
    size_t other_axis = (axis + 1) % 2;
    std::vector<MatrixEntry> entries;
    for (const auto& t: terms) {
        size_t loc[2];
        loc[axis] = t.dof;
        loc[other_axis] = original_entry.loc[other_axis];
        entries.push_back(MatrixEntry{{loc[0], loc[1]}, t.weight});
    }
    return entries;
}

std::vector<MatrixEntry>
add_entry_with_constraints(const ConstraintMatrix& constraint_matrix,
    const MatrixEntry& entry) 
{
    int constrained_axis;
    if (!is_constrained(constraint_matrix, entry.loc[0])) {
        if (!is_constrained(constraint_matrix, entry.loc[1])) {
            return {entry};
        }
        constrained_axis = 1;
    } else {
        constrained_axis = 0;
    }

    auto applied = add_term_with_constraints(constraint_matrix,
        LinearTerm{entry.loc[constrained_axis], entry.value});
    auto entries = entries_from_terms(applied, constrained_axis, entry);

    std::vector<MatrixEntry> out;
    for (const auto& e: entries) {
        auto recurse_entries = add_entry_with_constraints(constraint_matrix, e);
        for (const auto& e_recurse: recurse_entries) {
            out.push_back(e_recurse);
        }
    }
    
    return out;
}

Matrix
condense_matrix(const ConstraintMatrix& constraint_matrix, const Matrix& matrix)
{
    std::vector<double> condensed(matrix.n_rows * matrix.n_cols, 0.0);
    for (size_t row_idx = 0; row_idx < matrix.n_rows; row_idx++) {
        for (size_t col_idx = 0; col_idx < matrix.n_cols; col_idx++) {
            MatrixEntry entry_to_add{
                {row_idx, col_idx},
                matrix.data[row_idx * matrix.n_cols + col_idx]
            };

            auto expanded_entry =
                add_entry_with_constraints(constraint_matrix, entry_to_add);
            for (const auto& e: expanded_entry) {
                condensed[e.loc[0] * matrix.n_cols + e.loc[1]] += e.value;
            }
        }
    }

    std::vector<double> out;
    size_t n_rows = 0;
    size_t n_cols = 0;
    for (size_t row_idx = 0; row_idx < matrix.n_rows; row_idx++) {
        if (is_constrained(constraint_matrix, row_idx)) {
            continue;
        }
        n_rows++;
        for (size_t col_idx = 0; col_idx < matrix.n_cols; col_idx++) {
            if (is_constrained(constraint_matrix, col_idx)) {
                continue;
            }
            if (n_rows == 1) {
                n_cols++;
            }
            out.push_back(condensed[row_idx * matrix.n_cols + col_idx]);
        }
    }

    return {n_rows, n_cols, out};
}

} //END namespace tbem

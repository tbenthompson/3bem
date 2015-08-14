#include "constraint_matrix.h"
#include <algorithm>
#include <map>
#include <set>

namespace tbem {

bool is_constrained(const ConstraintMatrixData& cm, size_t dof) 
{
    auto it = cm.find(dof);
    if (it == cm.end()) {
        return false;
    }
    return true;
}

ConstraintEQ make_lower_triangular(const ConstraintEQ& c,
    const ConstraintMatrixData& matrix) 
{
    if (c.terms.size() == 0) {
        return c;
    }
    auto last_dof_index = find_last_dof_index(c);
    auto last_dof = c.terms[last_dof_index].dof;
    if (is_constrained(matrix, last_dof)) {
        auto last_dof_constraint = matrix.find(last_dof)->second;
        ConstraintEQ c_subs = substitute(c, last_dof_index, last_dof_constraint);
        ConstraintEQ c_subs_filtered = filter_zero_terms(c_subs);
        return make_lower_triangular(c_subs_filtered, matrix);
    }
    return c;
}


ConstraintMatrix from_constraints(const std::vector<ConstraintEQ>& constraints) 
{
    ConstraintMatrixData out;

    for (size_t i = 0; i < constraints.size(); i++) {
        auto combined = combine_terms(constraints[i]);
        auto zeros_filtered = filter_zero_terms(combined);
        auto lower_tri_constraint = make_lower_triangular(zeros_filtered, out);
        if (lower_tri_constraint.terms.size() == 0) {
            continue;
        }
        auto last_dof_index = find_last_dof_index(lower_tri_constraint);
        auto separated = isolate_term_on_lhs(lower_tri_constraint, last_dof_index);
        out.insert(std::make_pair(separated.constrained_dof, separated));
    }

    return {out};
};

std::vector<double> distribute_vector(const ConstraintMatrix& matrix,
    const std::vector<double>& in, size_t n_total_dofs) 
{
    std::vector<double> out(n_total_dofs); 

    size_t next_reduced_dof = 0;

    for (size_t dof_index = 0; dof_index < n_total_dofs; dof_index++) {
        if(is_constrained(matrix.map, dof_index)) {
            continue;
        }
        out[dof_index] = in[next_reduced_dof];
        next_reduced_dof++;
    }

    for (size_t dof_index = 0; dof_index < n_total_dofs; dof_index++) {
        if(!is_constrained(matrix.map, dof_index)) {
            continue;
        }
        auto constraint = matrix.map.find(dof_index)->second;

        auto val = constraint.rhs;
        for (size_t j = 0; j < constraint.terms.size(); j++) {
            assert(constraint.terms[j].dof < dof_index);
            val += constraint.terms[j].weight * out[constraint.terms[j].dof];
        }
        out[dof_index] = val;
    }
    return out;
}

void add_term_with_constraints(const ConstraintMatrixData& matrix,
    std::vector<double>& modifiable_vec, const LinearTerm& entry) 
{
    if (!is_constrained(matrix, entry.dof)) {
        modifiable_vec[entry.dof] += entry.weight;
        return;
    }

    auto constraint = matrix.find(entry.dof)->second;
    auto terms = constraint.terms;
    std::vector<LinearTerm> out_terms;
    for (auto t: terms) {
        double recurse_weight = t.weight * entry.weight;
        modifiable_vec[t.dof] += recurse_weight;
    }
}

std::vector<double> condense_vector(const ConstraintMatrix& matrix,
    const std::vector<double>& all) 
{
    std::vector<double> condensed_dofs(all.size(), 0.0);
    for (int dof_idx = all.size() - 1; dof_idx >= 0; dof_idx--) {
        double condensed_value = condensed_dofs[dof_idx];
        condensed_dofs[dof_idx] = 0.0;
        LinearTerm term_to_add{(size_t)dof_idx, all[dof_idx] + condensed_value};
        add_term_with_constraints(matrix.map, condensed_dofs, term_to_add);
    }

    std::vector<double> out(all.size() - matrix.size());
    size_t next_out_idx = 0;
    for (size_t dof_idx = 0; dof_idx < all.size(); dof_idx++) {
        if (is_constrained(matrix.map, dof_idx)) {
            continue;
        }
        assert(next_out_idx < out.size());
        out[next_out_idx] = condensed_dofs[dof_idx];
        next_out_idx++;
    }

    return out;
}

template <typename T>
void add_entry_with_constraints(const ConstraintMatrixData& row_cm, 
    const ConstraintMatrixData& col_cm, size_t n_rows, size_t n_cols,
    T& modifiable_matrix, const MatrixEntry& entry) 
{
    int constrained_axis;
    int other_axis;
    const ConstraintMatrixData* cm;
    if (!is_constrained(row_cm, entry.loc[0])) {
        if (!is_constrained(col_cm, entry.loc[1])) {
            size_t entry_idx = entry.loc[0] * n_cols + entry.loc[1];
            // modifiable_matrix[entry_idx] += entry.value;
            modifiable_matrix[entry_idx] += entry.value;
            return;
        }
        cm = &col_cm;
        constrained_axis = 1;
        other_axis = 0;
    } else {
        cm = &row_cm;
        constrained_axis = 0;
        other_axis = 1;
    }

    auto constrained_dof = entry.loc[constrained_axis];
    auto constraint = cm->find(constrained_dof)->second;
    auto terms = constraint.terms;

    for (auto t: terms) {
        double recurse_weight = t.weight * entry.value;
        size_t loc[2];
        loc[constrained_axis] = t.dof;
        loc[other_axis] = entry.loc[other_axis];
        assert(loc[constrained_axis] < constrained_dof);
        add_entry_with_constraints(
            row_cm, col_cm, n_rows, n_cols, modifiable_matrix, 
            {loc[0], loc[1], recurse_weight}
        );
    }
}


DenseOperator remove_constrained(const ConstraintMatrixData& row_cm,
    const ConstraintMatrixData& col_cm, const DenseOperator& matrix) 
{
    assert(matrix.n_rows() >= row_cm.size());
    assert(matrix.n_cols() >= col_cm.size());

    auto n_rows_out = matrix.n_rows() - row_cm.size();
    auto n_cols_out = matrix.n_cols() - col_cm.size();
    auto out = DenseOperator(n_rows_out, n_cols_out);

    size_t out_row_idx = 0;
    for (size_t in_row_idx = 0; in_row_idx < matrix.n_rows(); in_row_idx++) {
        if (is_constrained(row_cm, in_row_idx)) {
            continue;
        }
        size_t out_col_idx = 0;
        for (size_t in_col_idx = 0; in_col_idx < matrix.n_cols(); in_col_idx++) {
            if (is_constrained(col_cm, in_col_idx)) {
                continue;
            }
            auto in_entry_idx = in_row_idx * matrix.n_cols() + in_col_idx;
            auto out_entry_idx = out_row_idx * n_cols_out + out_col_idx;
            out[out_entry_idx] = matrix[in_entry_idx];

            out_col_idx++;
        }
        out_row_idx++;
    }

    return out;
}

DenseOperator condense_matrix(const ConstraintMatrix& row_cm,
    const ConstraintMatrix& col_cm, const DenseOperator& matrix)
{
    auto condensed = DenseOperator(matrix.n_rows(), matrix.n_cols(), 0.0);

    for (int row_idx = matrix.n_rows() - 1; row_idx >= 0; --row_idx) {
        for (int col_idx = matrix.n_cols() - 1; col_idx >= 0; --col_idx) {
            auto entry_idx = row_idx * matrix.n_cols() + col_idx;
            MatrixEntry entry_to_add{
                (size_t)row_idx, (size_t)col_idx,
                matrix[entry_idx]
            };

            add_entry_with_constraints(
                row_cm.map, col_cm.map, matrix.n_rows(), matrix.n_cols(),
                condensed, entry_to_add
            );
        }
    }
    
    auto out = remove_constrained(row_cm.map, col_cm.map, condensed);
    return out;
}

//TODO: This could be used in dense condensation and vector condensation too.
std::map<size_t,size_t> constraint_dof_map(const ConstraintMatrixData& cm, size_t n_dofs) 
{
    std::map<size_t,size_t> dof_map;
    size_t out_idx = 0;
    for (size_t in_idx = 0; in_idx < n_dofs; in_idx++) {
        if (is_constrained(cm, in_idx)) {
            continue;
        }
        dof_map[in_idx] = out_idx;
        out_idx++;
    }
    return dof_map;
}

std::vector<MatrixEntry> remove_constrained(const ConstraintMatrixData& row_cm,
    const ConstraintMatrixData& col_cm, size_t n_rows, size_t n_cols, 
    const std::map<size_t,double>& matrix) 
{
    std::vector<MatrixEntry> out;
    auto row_map = constraint_dof_map(row_cm, n_rows);
    auto col_map = constraint_dof_map(col_cm, n_cols);
    for (auto it = matrix.begin(); it != matrix.end(); ++it) {
        //TODO: General purpose DOK-->COO sparse matrix function
        auto matrix_idx = it->first;
        auto in_col = matrix_idx % n_cols;
        auto in_row = (matrix_idx - in_col) / n_cols;

        auto out_row = row_map[in_row];
        auto out_col = col_map[in_col];
        out.push_back({out_row, out_col, it->second});
    }
    return out;
}



SparseOperator condense_matrix(const ConstraintMatrix& row_cm,
    const ConstraintMatrix& col_cm, const SparseOperator& matrix)
{
    std::map<size_t,double> condensed;

    auto n_in_rows = matrix.n_rows();
    auto n_in_cols = matrix.n_cols();

    size_t count = 0;
    for (int row = n_in_rows - 1; row >= 0; --row) {
        int first_col = static_cast<int>(matrix.row_ptrs[row]);
        int last_col = static_cast<int>(matrix.row_ptrs[row + 1]);
        for (int col_idx = last_col - 1; col_idx >= first_col; --col_idx) {
            auto col = matrix.column_indices[col_idx];
            count++;
            auto val = matrix.values[col_idx];
            MatrixEntry entry_to_add{(size_t)row, (size_t)col, val};
            add_entry_with_constraints(
                row_cm.map, col_cm.map, n_in_rows, n_in_cols,
                condensed, entry_to_add
            );
        }
    }
    
    auto entries = remove_constrained(
        row_cm.map, col_cm.map, n_in_rows, n_in_cols, condensed
    );
    auto out_rows = n_in_rows - row_cm.map.size();
    auto out_cols = n_in_cols - col_cm.map.size();
    auto out = SparseOperator::csr_from_coo(out_rows, out_cols, entries);
    return out;
}

ConstraintMatrix homogenize_constraints(const ConstraintMatrix& cm)
{
    ConstraintMatrixData out;
    for (auto it = cm.map.begin(); it != cm.map.end(); ++it) {
        out.insert(std::make_pair(it->first, RearrangedConstraintEQ{
            it->first, it->second.terms, 0.0
        }));
    }
    return {out};
}

bool is_ignorable(const std::set<size_t>& ignored, const RearrangedConstraintEQ& c)
{
    for (auto t: c.terms) {
        if (ignored.count(t.dof) == 0) {
            return false; 
        }
    }
    return true;
}

size_t max_constrained_dof(const ConstraintMatrixData& cm) 
{
    size_t max_dof = 0;
    for (auto it = cm.begin(); it != cm.end(); ++it) {
        max_dof = std::max(it->first, max_dof);
    }
    return max_dof;
}

std::set<size_t> _identify_ignored_dofs_set(const ConstraintMatrix& cm)
{
    std::set<size_t> ignored;

    for (size_t dof_index = 0; dof_index <= max_constrained_dof(cm.map); dof_index++) {
        if(!is_constrained(cm.map, dof_index)) {
            continue;
        }
        auto constraint = cm.map.find(dof_index)->second;
        if (is_ignorable(ignored, constraint)) {
            ignored.insert(constraint.constrained_dof);
        }
    }
    return ignored;
}

std::vector<size_t> identify_ignored_dofs(const ConstraintMatrix& cm)
{
    auto ignored = _identify_ignored_dofs_set(cm);
    std::vector<size_t> out(ignored.begin(), ignored.end());
    return out;
}

DenseOperator distribute_row_zeros(const DenseOperator& matrix, 
    const ConstraintMatrix& cm)
{
    auto ignored_rows = _identify_ignored_dofs_set(cm);
    auto n_total_rows = matrix.n_rows() + ignored_rows.size();
    auto next_in_row = 0;
    std::vector<double> out_data(n_total_rows * matrix.n_cols());
    for (size_t i = 0; i < n_total_rows; i++) {
        if (ignored_rows.count(i) > 0) {
            for (size_t j = 0; j < matrix.n_cols(); j++) {
                out_data[i * matrix.n_cols() + j] = 0.0;
            }
        } else {
            for (size_t j = 0; j < matrix.n_cols(); j++) {
                auto value = matrix[next_in_row * matrix.n_cols() + j];
                out_data[i * matrix.n_cols() + j] = value;
            }
            next_in_row++;
        }
    }
    return DenseOperator(n_total_rows, matrix.n_cols(), out_data);
}

} // END namespace tbem

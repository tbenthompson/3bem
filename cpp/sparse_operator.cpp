#include "sparse_operator.h"
#include "dense_operator.h"

namespace tbem {

SparseOperator::SparseOperator(size_t n_rows, size_t n_cols,
    const std::vector<double>& values,
    const std::vector<size_t>& column_indices,
    const std::vector<size_t>& row_ptrs):
    shape{n_rows, n_cols},
    values(values),
    column_indices(column_indices),
    row_ptrs(row_ptrs)
{
    assert(row_ptrs.size() == n_rows + 1);
}

std::vector<double> SparseOperator::apply(const std::vector<double>& x) const 
{
    assert(x.size() == n_cols());
    std::vector<double> out(n_rows(), 0.0);
    for (size_t i = 0; i < row_ptrs.size() - 1; i++) {
        for (size_t c_idx = row_ptrs[i]; c_idx < row_ptrs[i + 1]; c_idx++) {
            out[i] += values[c_idx] * x[column_indices[c_idx]];
        }
    }
    return out;
}


std::vector<double> SparseOperator::to_dense() const
{
    std::vector<double> dense(n_cols() * n_rows(), 0.0);
    for (size_t i = 0; i < row_ptrs.size() - 1; i++) {
       for (size_t c_idx = row_ptrs[i]; c_idx < row_ptrs[i + 1]; c_idx++) {
           dense[i * n_cols() + column_indices[c_idx]] += values[c_idx];
       }
    }
    return dense;
}

DenseOperator SparseOperator::left_multiply_with_dense(const DenseOperator& other) const
{
    assert(n_rows() == other.n_cols());

    std::vector<double> out_matrix(other.n_rows() * n_cols(), 0.0);
    for (size_t i = 0; i < other.n_rows(); i++) {
        for (size_t j = 0; j < n_rows(); j++) {
            for (size_t c_idx = row_ptrs[j]; c_idx < row_ptrs[j + 1]; c_idx++) {
                auto out_col = column_indices[c_idx];
                auto out_idx = i * n_cols() + out_col;
                auto other_idx = i * other.n_cols() + j;
                out_matrix[out_idx] += other[other_idx] * values[c_idx];
            }
        }
    }
    return DenseOperator(other.n_rows(), n_cols(), out_matrix);
}

SparseOperator SparseOperator::csr_from_coo(size_t n_rows, size_t n_cols,
    const std::vector<MatrixEntry>& entries)
{
    //compute number of non-zero entries per row of A 
    std::vector<size_t> nonzeros_per_row(n_rows, 0);

    for (size_t i = 0; i < entries.size(); i++) {
        assert(entries[i].row() < n_rows);
        assert(entries[i].col() < n_cols);
        nonzeros_per_row[entries[i].row()]++;
    }

    std::vector<size_t> row_ptrs(n_rows + 1);
    size_t cumulative = 0;
    for (size_t i = 0; i < row_ptrs.size() - 1; i++) {
        row_ptrs[i] = cumulative;
        cumulative += nonzeros_per_row[i];
    }
    row_ptrs[n_rows] = entries.size();

    std::vector<size_t> in_row_already(n_rows, 0);
    std::vector<size_t> column_indices(entries.size());
    std::vector<double> values(entries.size());
    for (size_t i = 0; i < entries.size(); i++) {
        auto row = entries[i].row();
        auto row_next = row_ptrs[row] + in_row_already[row];
        column_indices[row_next] = entries[i].col(); 
        values[row_next] = entries[i].value;
        in_row_already[row]++;
    }

    return SparseOperator(n_rows, n_cols, values, column_indices, row_ptrs);
}

} //end namespace tbem

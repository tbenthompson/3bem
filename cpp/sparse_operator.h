#ifndef __1231231231898972_SPARSE_OPERATOR_H
#define __1231231231898972_SPARSE_OPERATOR_H
#include <vector>
#include <assert.h>
#include <iostream>
#include "operator.h"

namespace tbem {

struct MatrixEntry 
{
    size_t loc[2];
    double value;

    size_t row() const {return loc[0];}
    size_t col() const {return loc[1];}
};

struct SparseOperator: public OperatorI
{
    //TODO: get rid of the two types of shape. one type of shape!
    const OperatorShape shape;
    const std::vector<double> values;
    const std::vector<size_t> column_indices;
    const std::vector<size_t> row_ptrs;

    SparseOperator(size_t n_rows, size_t n_cols,
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

    size_t nnz() const {return row_ptrs.back();}

    std::vector<double> to_dense() {
        std::vector<double> dense(n_cols() * n_rows(), 0.0);
        for (size_t i = 0; i < row_ptrs.size() - 1; i++) {
           for (size_t c_idx = row_ptrs[i]; c_idx < row_ptrs[i + 1]; c_idx++) {
               dense[i * n_cols() + column_indices[c_idx]] += values[c_idx];
           }
        }
        return dense;
    }

    virtual size_t n_rows() const {return shape.n_rows;}
    virtual size_t n_cols() const {return shape.n_cols;}

    virtual std::vector<double> apply(const std::vector<double>& x) const 
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

    static SparseOperator csr_from_coo(
        size_t n_rows, size_t n_cols,
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
};


} // end namespace tbem

#endif

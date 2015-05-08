#ifndef __1231231231898972_SPARSE_OPERATOR_H
#define __1231231231898972_SPARSE_OPERATOR_H
#include <vector>
#include "block_operator.h"

namespace tbem {

struct MatrixEntry 
{
    size_t loc[2];
    double value;

    size_t row() const {return loc[0];}
    size_t col() const {return loc[1];}
};

struct BlockSparseOperator: public BlockOperatorI
{
    const OperatorShape component_shape;
    const OperatorShape block_shape;
    const std::vector<double> values;
    const std::vector<size_t> column_indices;
    const std::vector<size_t> row_ptrs;

    BlockSparseOperator(size_t n_component_rows, size_t n_component_cols,
        size_t n_block_rows, size_t n_block_cols,
        const std::vector<double>& values,
        const std::vector<size_t>& column_indices,
        const std::vector<size_t>& row_ptrs):
        component_shape{n_component_rows, n_component_cols},
        block_shape{n_block_rows, n_block_cols},
        values(values),
        column_indices(column_indices),
        row_ptrs(row_ptrs)
    {
        assert(row_ptrs.size() == n_total_rows() + 1);
    }

    size_t nnz() const {return row_ptrs.back();}

    virtual size_t n_block_rows() const {return block_shape.n_rows;}
    virtual size_t n_block_cols() const {return block_shape.n_cols;}
    virtual size_t n_total_rows() const 
    {
        return block_shape.n_rows * component_shape.n_rows;
    }
    virtual size_t n_total_cols() const
    {
        return block_shape.n_cols * component_shape.n_cols;
    }

    virtual std::vector<double> apply(const std::vector<double>& x) const 
    {
        assert(x.size() == n_total_cols());
        std::vector<double> out(n_total_rows(), 0.0);
        for (size_t i = 0; i < row_ptrs.size() - 1; i++) {
           for (size_t c_idx = row_ptrs[i]; c_idx < row_ptrs[i + 1]; c_idx++) {
                out[i] += values[c_idx] * x[column_indices[c_idx]];
           }
        }
        return out;
    }

    static BlockSparseOperator csr_from_coo(
        size_t n_component_rows, size_t n_component_cols,
        size_t n_block_rows, size_t n_block_cols,
        const std::vector<MatrixEntry>& entries)
    {
        //compute number of non-zero entries per row of A 
        auto n_total_rows = n_block_rows * n_component_rows;
        std::vector<size_t> nonzeros_per_row(n_total_rows, 0.0);

        for (size_t i = 0; i < entries.size(); i++) {
            nonzeros_per_row[entries[i].row()]++;
        }

        std::vector<size_t> row_ptrs(n_total_rows + 1);
        size_t cumulative = 0;
        for (size_t i = 0; i < row_ptrs.size(); i++) {
            row_ptrs[i] = cumulative;
            cumulative += nonzeros_per_row[i];
        }
        row_ptrs[n_total_rows] = entries.size();

        std::vector<size_t> in_row_already(n_total_rows, 0);
        std::vector<size_t> column_indices(entries.size());
        std::vector<double> values(entries.size());
        for (size_t i = 0; i < entries.size(); i++) {
            auto row = entries[i].row();
            auto row_next = row_ptrs[row] + in_row_already[row];
            column_indices[row_next] = entries[i].col(); 
            values[row_next] = entries[i].value;
            in_row_already[row]++;
        }
        return BlockSparseOperator(
            n_component_rows, n_component_cols, n_block_rows, n_block_cols,
            values, column_indices, row_ptrs
        );
    }
};


} // end namespace tbem

#endif

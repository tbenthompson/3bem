#ifndef __1231231231898972_SPARSE_OPERATOR_H
#define __1231231231898972_SPARSE_OPERATOR_H
#include <vector>
#include "block_operator.h"

namespace tbem {

struct SparseMatrixEntry 
{
    const size_t row;
    const size_t col;
    const double value;
};

struct BlockSparseOperator: public BlockOperatorI
{
    const OperatorShape component_shape;
    const OperatorShape block_shape;
    const std::vector<SparseMatrixEntry> storage;

    BlockSparseOperator(size_t n_component_rows, size_t n_component_cols,
        size_t n_block_rows, size_t n_block_cols,
        const std::vector<SparseMatrixEntry>& entries):
        component_shape{n_component_rows, n_component_cols},
        block_shape{n_block_rows, n_block_cols},
        storage(entries)
    {}

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
        for (size_t i = 0; i < storage.size(); i++) {
            out[storage[i].row] += storage[i].value * x[storage[i].col];
        }
        return out;
    }
};

} // end namespace tbem

#endif

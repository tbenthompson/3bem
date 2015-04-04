#ifndef __1231231231898972_SPARSE_OPERATOR_H
#define __1231231231898972_SPARSE_OPERATOR_H
#include <vector>
#include "vectorx.h"
#include "operator.h"

namespace tbem {

struct SparseMatrixEntry 
{
    const size_t loc[2];
    const double value;
};

struct SparseOperator: public OperatorI
{
    typedef std::vector<SparseMatrixEntry> DataT;

    const OperatorShape shape;
    const DataT storage;

    SparseOperator(size_t n_rows, size_t n_cols,
        const std::vector<SparseMatrixEntry>& entries):
        shape{n_rows, n_cols},
        storage(entries)
    {}

    virtual size_t n_rows() const 
    {
        return shape.n_rows; 
    }

    virtual size_t n_cols() const
    {
        return shape.n_cols;
    }

    virtual VectorX apply(const VectorX& x) const {
        VectorX out(shape.n_rows, 0.0);
        for (size_t i = 0; i < storage.size(); i++) {
            out[storage[i].loc[0]] += storage[i].value * x[storage[i].loc[1]];
        }
        return out;
    }
};

template <typename T>
struct BlockOperator;
typedef BlockOperator<SparseOperator> BlockSparseOperator;

} // end namespace tbem

#endif

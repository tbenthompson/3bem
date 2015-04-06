#ifndef __1234567876161_BLOCK_OPERATOR_H
#define __1234567876161_BLOCK_OPERATOR_H

#include <vector>
#include <cassert>
#include "vectorx.h"
#include "operator.h"

namespace tbem {

struct BlockOperatorI {
    virtual size_t n_block_rows() const = 0;
    virtual size_t n_block_cols() const = 0;
    virtual size_t n_total_rows() const = 0;
    virtual size_t n_total_cols() const = 0;
    virtual BlockVectorX apply(const BlockVectorX& x) const = 0;

    VectorX apply_scalar(const VectorX& x) const {
        assert(n_block_cols() == 1);
        assert(n_block_rows() == 1);
        return apply(BlockVectorX{x})[0]; 
    }
    virtual VectorX apply(const VectorX& x) const {
        return apply_scalar(x);
    }
};

template <typename T>
struct BlockOperator: public BlockOperatorI {
    const OperatorShape shape;
    std::vector<T> ops;

    BlockOperator(size_t n_block_rows, size_t n_block_cols, std::vector<T> ops):
        shape{n_block_rows, n_block_cols},
        ops(std::move(ops))
    {}

    virtual size_t n_block_rows() const {return shape.n_rows;}
    virtual size_t n_block_cols() const {return shape.n_cols;}
   
    virtual size_t n_total_rows() const {
        size_t n_rows = 0;
        for (size_t d1 = 0; d1 < n_block_rows(); d1++) {
            n_rows += ops[d1 * n_block_cols()].n_rows();
        }
        return n_rows;
    }

    virtual size_t n_total_cols() const {
        size_t n_cols = 0;
        for (size_t d2 = 0; d2 < n_block_cols(); d2++) {
            n_cols += ops[d2].n_cols();
        }
        return n_cols;
    }

    virtual BlockVectorX apply(const BlockVectorX& x) const {
        assert(n_block_rows() * x.size() == ops.size());
        assert(n_block_cols() == x.size());

        BlockVectorX res(n_block_rows(), VectorX(ops[0].n_rows(), 0.0));
        for (size_t d1 = 0; d1 < n_block_rows(); d1++) {
            for (size_t d2 = 0; d2 < n_block_cols(); d2++) {
                size_t comp_idx = d1 * n_block_cols() + d2;

                const auto& op = ops[comp_idx];
                res[d1] += op.apply(x[d2]);
            }
        }
        return res;
    }
    
    const T& get_block(size_t row, size_t col) const {
        return ops[row * n_block_cols() + col];
    }
};


}//end namespace tbem

#endif

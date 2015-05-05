#ifndef __1234567876161_BLOCK_OPERATOR_H
#define __1234567876161_BLOCK_OPERATOR_H

#include <vector>
#include <cassert>
#include "operator.h"

namespace tbem {

struct BlockOperatorI {
    virtual size_t n_block_rows() const = 0;
    virtual size_t n_block_cols() const = 0;
    virtual size_t n_total_rows() const = 0;
    virtual size_t n_total_cols() const = 0;
    virtual std::vector<double> apply(const std::vector<double>& x) const = 0;
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

    virtual std::vector<double> apply(const std::vector<double>& x) const {
        assert(x.size() == n_total_cols());

        std::vector<double> res(n_block_rows() * ops[0].n_rows(), 0.0);
        for (size_t d1 = 0; d1 < n_block_rows(); d1++) {
            for (size_t d2 = 0; d2 < n_block_cols(); d2++) {
                size_t comp_idx = d1 * n_block_cols() + d2;

                const auto& op = ops[comp_idx];
                std::vector<double> input(ops[0].n_cols());
                for (size_t i = 0; i < ops[0].n_cols(); i++) {
                    input[i] = x[d2 * ops[0].n_cols() + i];
                }
                auto applied = op.apply(input);
                assert(applied.size() == ops[0].n_rows());
                for (size_t i = 0; i < ops[0].n_rows(); i++) {
                    res[d1 * ops[0].n_rows() + i] += applied[i];
                }
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

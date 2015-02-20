#include "dense_operator.h"

namespace tbem {

VectorX Operator::apply(const VectorX& x) const {
    assert(x.size() == n_cols());
   
    VectorX res(n_rows(), 0.0);
#pragma omp parallel for
    for (size_t i = 0; i < n_rows(); i++) {
        for (size_t j = 0; j < n_cols(); j++) {
            size_t matrix_idx = i * n_cols() + j;
            res[i] += (*storage)[matrix_idx] * x[j];
        }
    }

    return res;
}

size_t BlockOperator::n_block_rows() const {
    return shape.n_rows;
}

size_t BlockOperator::n_block_cols() const {
    return shape.n_cols;
}

size_t BlockOperator::n_total_rows() const {
    size_t n_rows = 0;
    for (size_t d1 = 0; d1 < n_block_rows(); d1++) {
        n_rows += ops[d1 * n_block_cols()].n_rows();
    }
    return n_rows;
}

size_t BlockOperator::n_total_cols() const {
    size_t n_cols = 0;
    for (size_t d2 = 0; d2 < n_block_cols(); d2++) {
        n_cols += ops[d2].n_cols();
    }
    return n_cols;
}

    
BlockVectorX BlockOperator::apply(const BlockVectorX& x) const {
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

Operator BlockOperator::combine_components() {
    Operator out(n_total_rows(), n_total_cols());

    size_t n_rows_so_far = 0;
    for (size_t d1 = 0; d1 < n_block_rows(); d1++) {
        auto n_this_comp_rows = ops[d1 * n_block_cols()].n_rows();
        for (size_t row_idx = 0; row_idx < n_this_comp_rows; row_idx++) {
            auto out_row_idx = n_rows_so_far + row_idx;

            size_t n_cols_so_far = 0;
            for (size_t d2 = 0; d2 < n_block_cols(); d2++) {
                auto in_comp = d1 * n_block_cols() + d2;

                const auto& op = ops[in_comp];
                for (size_t col_idx = 0; col_idx < op.n_cols(); col_idx++) {
                    auto out_col_idx = n_cols_so_far + col_idx;

                    auto in_element = row_idx * op.n_cols() + col_idx;
                    auto out_element = out_row_idx * n_total_cols() + out_col_idx;
                    out[out_element] = op[in_element];
                }
                n_cols_so_far += op.n_cols();
            }
        }
        n_rows_so_far += n_this_comp_rows;
    }

    return out;
}


} //end namespace tbem

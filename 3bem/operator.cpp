#include "operator.h"

namespace tbem {

VectorX Operator::apply(const VectorX& x) const {
    assert(x.size() == n_cols());
   
    VectorX res(n_rows(), 0.0);
#pragma omp parallel for
    for (size_t i = 0; i < n_rows(); i++) {
        for (size_t j = 0; j < n_cols(); j++) {
            size_t matrix_idx = i * n_cols() + j;
            res[i] += (*_data)[matrix_idx] * x[j];
        }
    }

    return res;
}

size_t BlockOperator::n_cols() const {
    size_t n_cols = 0;
    for (size_t d2 = 0; d2 < n_comp_cols; d2++) {
        n_cols += ops[d2].n_cols();
    }
    return n_cols;
}

size_t BlockOperator::n_rows() const {
    size_t n_rows = 0;
    for (size_t d1 = 0; d1 < n_comp_rows; d1++) {
        n_rows += ops[d1 * n_comp_cols].n_rows();
    }
    return n_rows;
}
    
BlockVectorX BlockOperator::apply(const BlockVectorX& x) const {
    assert(n_comp_rows * x.size() == ops.size());
    assert(n_comp_cols == x.size());

    BlockVectorX res(n_comp_rows, VectorX(ops[0].n_rows(), 0.0));
    for (size_t d1 = 0; d1 < n_comp_rows; d1++) {
        for (size_t d2 = 0; d2 < n_comp_cols; d2++) {
            size_t comp_idx = d1 * n_comp_cols + d2;

            const auto& op = ops[comp_idx];
            assert(op.n_rows() * x[d2].size() == op.data().size());
            assert(op.n_cols() == x[d2].size());

            res[d1] += op.apply(x[d2]);
        }
    }
    return res;
}

Operator combine_components(const BlockOperator& block_op) {
    auto n_cols = block_op.n_cols();
    auto n_rows = block_op.n_rows();

    Operator out(n_rows, n_cols);

    size_t n_rows_so_far = 0;
    for (size_t d1 = 0; d1 < block_op.n_comp_rows; d1++) {
        auto n_this_comp_rows = block_op.ops[d1 * block_op.n_comp_cols].n_rows();
        for (size_t row_idx = 0; row_idx < n_this_comp_rows; row_idx++) {
            auto out_row_idx = n_rows_so_far + row_idx;

            size_t n_cols_so_far = 0;
            for (size_t d2 = 0; d2 < block_op.n_comp_cols; d2++) {
                auto in_comp = d1 * block_op.n_comp_cols + d2;

                const auto& op = block_op.ops[in_comp];
                for (size_t col_idx = 0; col_idx < op.n_cols(); col_idx++) {
                    auto out_col_idx = n_cols_so_far + col_idx;

                    auto in_element = row_idx * op.n_cols() + col_idx;
                    auto out_element = out_row_idx * n_cols + out_col_idx;
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

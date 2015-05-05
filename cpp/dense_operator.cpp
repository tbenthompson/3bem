#include "dense_operator.h"
#include "block_operator.h"
#include <cassert>

namespace tbem {

DenseOperator::DenseOperator(size_t n_rows, size_t n_cols,
        const std::vector<double>& data):
    shape{n_rows, n_cols},
    storage(std::make_shared<std::vector<double>>(data))
{}

DenseOperator::DenseOperator(size_t n_rows, size_t n_cols):
    shape{n_rows, n_cols},
    storage(std::make_shared<std::vector<double>>(n_rows * n_cols))
{}

DenseOperator::DenseOperator(size_t n_rows, size_t n_cols, double val):
    shape{n_rows, n_cols},
    storage(std::make_shared<std::vector<double>>(n_rows * n_cols, val))
{}

size_t DenseOperator::n_rows() const {
    return shape.n_rows;
}

size_t DenseOperator::n_cols() const {
    return shape.n_cols;
}

size_t DenseOperator::n_elements() const {
    return shape.n_rows * shape.n_cols;
}

std::vector<double> DenseOperator::apply(const std::vector<double>& x) const {
    assert(x.size() == n_cols());
   
    std::vector<double> res(n_rows(), 0.0);
#pragma omp parallel for
    for (size_t i = 0; i < n_rows(); i++) {
        for (size_t j = 0; j < n_cols(); j++) {
            size_t matrix_idx = i * n_cols() + j;
            res[i] += (*storage)[matrix_idx] * x[j];
        }
    }

    return res;
}

const DenseOperator::DataT& DenseOperator::data() const {
    return *storage;
}

double& DenseOperator::operator[](size_t idx) {
    return (*storage)[idx];
}

const double& DenseOperator::operator[](size_t idx) const {
    return (*storage)[idx];
}

void copy_matrix(const DenseOperator& in_op, size_t out_start_row,
        size_t out_start_col, DenseOperator& out_op) 
{
    
    for (size_t row_idx = 0; row_idx < in_op.n_rows(); row_idx++) {
        auto out_row_idx = out_start_row + row_idx;

        for (size_t col_idx = 0; col_idx < in_op.n_cols(); col_idx++) {
            auto out_col_idx = out_start_col + col_idx;

            auto in_element = row_idx * in_op.n_cols() + col_idx;
            auto out_element = out_row_idx * out_op.n_cols() + out_col_idx;
            out_op[out_element] = in_op[in_element];
        }
    }
}

DenseOperator combine_components(const BlockDenseOperator& block) 
{
    DenseOperator out(block.n_total_rows(), block.n_total_cols());

    size_t n_rows_so_far = 0;
    for (size_t d1 = 0; d1 < block.n_block_rows(); d1++) {
        auto n_this_comp_rows = block.ops[d1 * block.n_block_cols()].n_rows();

        size_t n_cols_so_far = 0;
        for (size_t d2 = 0; d2 < block.n_block_cols(); d2++) {
            auto in_comp = d1 * block.n_block_cols() + d2;
            const auto& op = block.ops[in_comp];

            copy_matrix(op, n_rows_so_far, n_cols_so_far, out);

            n_cols_so_far += op.n_cols();
        }
        n_rows_so_far += n_this_comp_rows;
    }

    return out;
}


}//end namespace tbem

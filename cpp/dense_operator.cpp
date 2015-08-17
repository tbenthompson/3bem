#include "dense_operator.h"
#include <cassert>
#include "blas_wrapper.h"

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

size_t DenseOperator::n_rows() const 
{
    return shape.n_rows;
}

size_t DenseOperator::n_cols() const 
{
    return shape.n_cols;
}

size_t DenseOperator::n_elements() const 
{
    return shape.n_rows * shape.n_cols;
}

std::vector<double> DenseOperator::apply(const std::vector<double>& x) const 
{
    assert(x.size() == n_cols());
    if (n_cols() == 0) {
        return std::vector<double>(n_rows(), 0.0);
    }
    auto out = matrix_vector_product(*storage, x);
    return out;
}

const DenseOperator::DataT& DenseOperator::data() const 
{
    return *storage;
}

double& DenseOperator::operator[](size_t idx) 
{
    return (*storage)[idx];
}

const double& DenseOperator::operator[](size_t idx) const 
{
    return (*storage)[idx];
}

DenseOperator DenseOperator::add(const DenseOperator& B) 
{
    assert(n_rows() == B.n_rows());
    assert(n_cols() == B.n_cols());
    std::vector<double> out_data = data();
    for (size_t i = 0; i < out_data.size(); i++) {
        out_data[i] += B.data()[i];
    }
    return DenseOperator(n_rows(), n_cols(), out_data);
}

std::unique_ptr<OperatorI> DenseOperator::clone() const 
{
    return std::unique_ptr<DenseOperator>(new DenseOperator(
        shape.n_rows, shape.n_cols, *storage
    ));
}

DenseOperator compose_dense_ops(const std::vector<DenseOperator>& ops,
    const std::vector<size_t>& start_rows, const std::vector<size_t>& start_cols,
    const std::vector<double>& multipliers, size_t n_out_rows, size_t n_out_cols)
{
    DenseOperator result(n_out_rows, n_out_cols, 0.0);
    for (size_t i = 0; i < ops.size(); i++) {
        auto s_row = start_rows[i];
        auto s_col = start_cols[i];
        auto n_rows = ops[i].n_rows();
        auto n_cols = ops[i].n_cols();
        auto mult = multipliers[i];
        for (size_t row = 0; row < n_rows; row++) {
            for (size_t col = 0; col < n_cols; col++) {
                auto out_entry = (row + s_row) * n_out_cols + (col + s_col);
                auto in_entry = row * n_cols + col;
                result[out_entry] += ops[i][in_entry] * mult;
            }
        }
    }
    return result;
}

}//end namespace tbem

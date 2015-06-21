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

}//end namespace tbem

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
    auto out = matrix_vector_product(*storage, x);
    return out;
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

}//end namespace tbem

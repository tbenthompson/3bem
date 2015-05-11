#include "dense_operator.h"
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

}//end namespace tbem

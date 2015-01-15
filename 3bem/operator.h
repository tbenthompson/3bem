#ifndef __123123123789798798_OPERATOR_H
#define __123123123789798798_OPERATOR_H

#include <cassert>
#include <vector>
#include "numbers.h"
//TODO: Try to remove this dependency on vec_ops
#include "vec_ops.h"

namespace tbem {
    
template <typename I, typename O, typename Op>
struct MatrixOperator 
{
    typedef I InType;
    typedef O OutType;
    typedef Op OperatorType;

    const size_t rows;
    const size_t cols;
    const std::vector<OperatorType> data;
};

double apply_operator(double A, double x) {
    return A * x;
}

Vec<double,2> apply_operator(const Vec<Vec<double,2>,2>& A, const Vec<double,2>& x) {
    return {dot_product(A[0], x), dot_product(A[1], x)};
}

Vec<double,3> apply_operator(const Vec<Vec<double,3>,3>& A, const Vec<double,3>& x) {
    return {dot_product(A[0], x), dot_product(A[1], x), dot_product(A[2], x)};
}

template <typename InType, typename OutType, typename OperatorType>
std::vector<OutType>
apply_operator(const MatrixOperator<InType,OutType,OperatorType>& A,
    const std::vector<InType>& x) 
{
    assert(A.rows * x.size() == A.data.size());
    std::vector<OutType> res(A.rows, zeros<OutType>::make());
#pragma omp parallel for
    for (int i = 0; i < A.rows; i++) {
        for (size_t j = 0; j < x.size(); j++) {
            res[i] += apply_operator(A.data[i * x.size() + j], x[j]);
        }
    }
    return res;
}

} // end namespace tbem

#endif

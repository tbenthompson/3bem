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

template <typename T, size_t dim>
T apply_operator(const Vec<T,dim>& A, const Vec<double,dim>& x) {
    return dot_product(x, A);
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
// template <size_t dim, typename InType, typename OutType, typename OperatorType>
// Function<OutType,dim>
// apply_operator(const MatrixOperator<InType,OutType,OperatorType>& A,
//                const Function<InType,dim> x) 
// {
//     assert(A.rows * x.n_dofs() == A.data.size());
//     std::vector<
//     return Function<OutType,dim>{};
// //     assert(A.rows * x.size() == A.data.size());
// //     std::vector<OutType> res(A.rows, zeros<OutType>::make());
// // #pragma omp parallel for
// //     for (int i = 0; i < A.rows; i++) {
// //         for (size_t j = 0; j < x.size(); j++) {
// //             res[i] += apply_operator(A.data[i * x.size() + j], x[j]);
// //         }
// //     }
// //     return res;
// }

} // end namespace tbem

#endif

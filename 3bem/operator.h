#ifndef __123123123789798798_OPERATOR_H
#define __123123123789798798_OPERATOR_H

namespace tbem {
    
template <typename I, typename O, typename Op>
struct MatrixOperator 
{
    typedef I InType;
    typedef O OutType;
    typedef Op OperatorType;

    const size_t rows;
    const size_t cols;
    std::vector<OperatorType> data;
};

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
            res[i] += dot_product(x[j], A.data[i * x.size() + j]);
        }
    }
    return res;
}

} // end namespace tbem

#endif

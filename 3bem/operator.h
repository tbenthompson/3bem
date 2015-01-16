#ifndef __123123123789798798_OPERATOR_H
#define __123123123789798798_OPERATOR_H

#include <cassert>
#include <vector>
#include "numbers.h"
//TODO: Try to remove this dependency on vec_ops
#include "vec_ops.h"

namespace tbem {

struct MatrixOperator 
{
    //TODO: rename to n_rows
    const size_t rows;
    //TODO: rename to n_cols
    const size_t cols;
    const size_t n_comp_rows;
    const size_t n_comp_cols;
    const std::vector<std::vector<double>> data;
};

MatrixOperator 
reshape_to_operator(const size_t rows, const size_t cols, const std::vector<double>& A);

template <size_t dim>
MatrixOperator reshape_to_operator(const size_t rows, const size_t cols, 
    std::vector<Vec<Vec<double,dim>,dim>> A);

std::vector<std::vector<double>>
apply_operator(const MatrixOperator& A, const std::vector<std::vector<double>>& x);

std::vector<double> 
apply_operator(const MatrixOperator& A, const std::vector<double>& x); 

} // end namespace tbem

#endif

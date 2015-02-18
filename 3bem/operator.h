#ifndef __123123123789798798_OPERATOR_H
#define __123123123789798798_OPERATOR_H

#include <cstdlib>
#include <cassert>
#include <vector>
#include <memory>
#include "vec.h"
#include "function.h"

namespace tbem {

struct OperatorI {
    virtual size_t n_rows() const = 0;
    virtual size_t n_cols() const = 0;
    virtual Function apply_operator(const Function& x) const = 0;
};

struct DenseMatrixOperator: public OperatorI;
struct MatrixFreeFarfieldOperator: public OperatorI;
struct SparseMatrixFMMOperator: public OperatorI;
struct MatrixFreeFMMOperator: public OperatorI;
struct BlockOperatorI;

struct Operator {
    const size_t n_rows;
    const size_t n_cols;

    typedef std::shared_ptr<std::vector<double>> OperatorDataPtr;
    OperatorDataPtr internal_data;

    const std::vector<double>& data() const {
        return *internal_data;
    }

    double& operator[] (size_t idx) {
        return (*internal_data)[idx];
    }

    const double& operator[] (size_t idx) const {
        return (*internal_data)[idx];
    }
};

struct BlockOperator 
{
    const size_t n_comp_rows;
    const size_t n_comp_cols;
    std::vector<Operator> ops;
};

Operator make_operator(size_t n_rows, size_t n_cols, const std::vector<double>& data);

BlockOperator 
reshape_to_operator(const size_t rows, const size_t cols, const std::vector<double>& A);

template <size_t dim>
BlockOperator reshape_to_operator(const size_t rows, const size_t cols, 
    const std::vector<Vec<Vec<double,dim>,dim>>& A);

BlockFunction apply_operator(const BlockOperator& A, const BlockFunction& x);

Function apply_operator(const BlockOperator& A, const Function& x); 

BlockOperator combine_components(const BlockOperator& op);

} // end namespace tbem

#endif

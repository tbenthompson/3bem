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

// struct DenseMatrixOperator: public OperatorI;
// struct MatrixFreeFarfieldOperator: public OperatorI;
// struct SparseMatrixFMMOperator: public OperatorI;
// struct MatrixFreeFMMOperator: public OperatorI;
// struct BlockOperatorI;

struct Operator {
    const size_t n_rows;
    const size_t n_cols;

    typedef std::vector<double> DataT;
    //TODO: Consider changing to unique_ptr
    typedef std::shared_ptr<DataT> OperatorDataPtr;
    OperatorDataPtr internal_data;

    const DataT& data() const {
        return *internal_data;
    }

    double& operator[] (size_t idx) {
        return (*internal_data)[idx];
    }

    const double& operator[] (size_t idx) const {
        return (*internal_data)[idx];
    }

    static Operator empty(size_t n_rows, size_t n_cols) {
        return {n_rows, n_cols, std::make_shared<DataT>(n_rows * n_cols)};
    }

    static Operator constant(size_t n_rows, size_t n_cols, double val) {
        return {n_rows, n_cols, std::make_shared<DataT>(n_rows * n_cols, val)};
    }
};

struct BlockOperator 
{
    const size_t n_comp_rows;
    const size_t n_comp_cols;
    std::vector<Operator> ops;
};

//TODO: Make this a static member 
Operator make_operator(size_t n_rows, size_t n_cols, const std::vector<double>& data);

BlockFunction apply_operator(const BlockOperator& A, const BlockFunction& x);

Function apply_operator(const BlockOperator& A, const Function& x); 

BlockOperator combine_components(const BlockOperator& op);

} // end namespace tbem

#endif

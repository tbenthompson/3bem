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
    virtual Function apply(const Function& x) const = 0;
};

// struct DenseMatrixOperator: public OperatorI;
// struct MatrixFreeFarfieldOperator: public OperatorI;
// struct SparseMatrixFMMOperator: public OperatorI;
// struct MatrixFreeFMMOperator: public OperatorI;
// struct BlockOperatorI;

struct Operator: public OperatorI {
    typedef std::vector<double> DataT;
    //TODO: Consider changing to unique_ptr
    typedef std::shared_ptr<DataT> OperatorDataPtr;

    const size_t _n_rows;
    const size_t _n_cols;
    OperatorDataPtr _data;

    Operator(size_t n_rows, size_t n_cols, const std::vector<double>& data):
        _n_rows(n_rows), _n_cols(n_cols), 
        _data(std::make_shared<std::vector<double>>(data))
    {}

    Operator(size_t n_rows, size_t n_cols):
        _n_rows(n_rows), _n_cols(n_cols), 
        _data(std::make_shared<std::vector<double>>(n_rows * n_cols))
    {}

    Operator(size_t n_rows, size_t n_cols, double val):
        _n_rows(n_rows), _n_cols(n_cols), 
        _data(std::make_shared<std::vector<double>>(n_rows * n_cols, val))
    {}

    virtual size_t n_rows() const {return _n_rows;}
    virtual size_t n_cols() const {return _n_cols;}
    size_t n_elements() const {return _n_rows * _n_cols;}

    virtual Function apply(const Function& x) const;

    const DataT& data() const {
        return *_data;
    }

    double& operator[] (size_t idx) {
        return (*_data)[idx];
    }

    const double& operator[] (size_t idx) const {
        return (*_data)[idx];
    }
};

struct BlockOperator 
{
    const size_t n_comp_rows;
    const size_t n_comp_cols;
    std::vector<Operator> ops;
};

BlockFunction apply_operator(const BlockOperator& A, const BlockFunction& x);

Function apply_operator(const BlockOperator& A, const Function& x); 

BlockOperator combine_components(const BlockOperator& op);

} // end namespace tbem

#endif

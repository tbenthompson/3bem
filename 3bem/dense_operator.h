#ifndef __123123123789798798_DENSE_OPERATOR_H
#define __123123123789798798_DENSE_OPERATOR_H

#include <cstdlib>
#include <cassert>
#include <vector>
#include <memory>
#include "vec.h"
#include "function.h"
#include "operator.h"

namespace tbem {

struct OperatorShape {
    size_t n_rows;
    size_t n_cols;
};

struct Operator: public OperatorI {
    typedef std::vector<double> DataT;
    //TODO: Consider changing to unique_ptr
    typedef std::shared_ptr<DataT> OperatorDataPtr;

    OperatorShape shape;
    OperatorDataPtr storage;

    Operator(size_t n_rows, size_t n_cols, const std::vector<double>& data):
        shape{n_rows, n_cols},
        storage(std::make_shared<std::vector<double>>(data))
    {}

    Operator(size_t n_rows, size_t n_cols):
        shape{n_rows, n_cols},
        storage(std::make_shared<std::vector<double>>(n_rows * n_cols))
    {}

    Operator(size_t n_rows, size_t n_cols, double val):
        shape{n_rows, n_cols},
        storage(std::make_shared<std::vector<double>>(n_rows * n_cols, val))
    {}

    virtual size_t n_rows() const {return shape.n_rows;}
    virtual size_t n_cols() const {return shape.n_cols;}
    size_t n_elements() const {return shape.n_rows * shape.n_cols;}

    virtual VectorX apply(const VectorX& x) const;

    const DataT& data() const {
        return *storage;
    }

    double& operator[] (size_t idx) {
        return (*storage)[idx];
    }

    const double& operator[] (size_t idx) const {
        return (*storage)[idx];
    }
};

struct BlockOperator: public BlockOperatorI
{
    const OperatorShape shape;
    std::vector<Operator> ops;

    BlockOperator(size_t n_block_rows, size_t n_block_cols, 
            const std::vector<Operator>& ops):
        shape{n_block_rows, n_block_cols},
        ops(ops)
    {}

    Operator combine_components();

    virtual size_t n_block_rows() const;
    virtual size_t n_block_cols() const;
    virtual size_t n_total_rows() const;
    virtual size_t n_total_cols() const;
    virtual BlockVectorX apply(const BlockVectorX& x) const;
};


// struct DenseMatrixOperator: public OperatorI;
// struct MatrixFreeFarfieldOperator: public OperatorI;
// struct SparseMatrixFMMOperator: public OperatorI;
// struct MatrixFreeFMMOperator: public OperatorI;
// struct BlockOperatorI;


} // end namespace tbem

#endif

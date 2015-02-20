#ifndef __123123123789798798_DENSE_OPERATOR_H
#define __123123123789798798_DENSE_OPERATOR_H

#include <vector>
#include <memory>
#include "fwd_vectorx.h"
#include "operator.h"

namespace tbem {

struct OperatorShape {
    size_t n_rows;
    size_t n_cols;
};

struct DenseOperator: public OperatorI {
    typedef std::vector<double> DataT;
    //TODO: Consider changing to unique_ptr
    typedef std::shared_ptr<DataT> OperatorDataPtr;

    OperatorShape shape;
    OperatorDataPtr storage;

    DenseOperator(size_t n_rows, size_t n_cols, const std::vector<double>& data);
    DenseOperator(size_t n_rows, size_t n_cols);
    DenseOperator(size_t n_rows, size_t n_cols, double val);

    virtual size_t n_rows() const;
    virtual size_t n_cols() const;
    size_t n_elements() const; 

    virtual VectorX apply(const VectorX& x) const;

    const DataT& data() const;
    double& operator[] (size_t idx); 
    const double& operator[] (size_t idx) const;
};

struct BlockDenseOperator: public BlockOperatorI
{
    const OperatorShape shape;
    std::vector<DenseOperator> ops;

    BlockDenseOperator(size_t n_block_rows, size_t n_block_cols, 
                       const std::vector<DenseOperator>& ops);

    DenseOperator combine_components();

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

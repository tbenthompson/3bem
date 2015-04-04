#ifndef __123123123789798798_DENSE_OPERATOR_H
#define __123123123789798798_DENSE_OPERATOR_H

#include <vector>
#include <memory>
#include "fwd_vectorx.h"
#include "operator.h"

namespace tbem {


struct DenseOperator: public OperatorI {
    typedef std::vector<double> DataT;
    //TODO: Consider changing to unique_ptr
    typedef std::shared_ptr<DataT> OperatorDataPtr;

    const OperatorShape shape;
    const OperatorDataPtr storage;

    DenseOperator(size_t n_rows, size_t n_cols, const std::vector<double>& data);
    DenseOperator(size_t n_rows, size_t n_cols);
    DenseOperator(size_t n_rows, size_t n_cols, double val);

    virtual size_t n_rows() const override;
    virtual size_t n_cols() const override;
    size_t n_elements() const; 

    virtual VectorX apply(const VectorX& x) const override;

    const DataT& data() const;
    double& operator[] (size_t idx); 
    const double& operator[] (size_t idx) const;
};

template <typename T>
struct BlockOperator;
typedef BlockOperator<DenseOperator> BlockDenseOperator;

DenseOperator combine_components(const BlockDenseOperator& block);

} // end namespace tbem

#endif

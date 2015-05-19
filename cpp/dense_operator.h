#ifndef __123123123789798798_DENSE_OPERATOR_H
#define __123123123789798798_DENSE_OPERATOR_H

#include <vector>
#include <memory>
#include "operator.h"

namespace tbem {


struct DenseOperator: public OperatorI {
    typedef std::vector<double> DataT;
    typedef std::shared_ptr<DataT> OperatorDataPtr;

    const OperatorShape shape;
    const OperatorDataPtr storage;

    DenseOperator(size_t n_rows, size_t n_cols, const std::vector<double>& data);
    DenseOperator(size_t n_rows, size_t n_cols);
    DenseOperator(size_t n_rows, size_t n_cols, double val);

    virtual size_t n_rows() const;
    virtual size_t n_cols() const;
    size_t n_elements() const; 

    virtual std::vector<double> apply(const std::vector<double>& x) const override;

    const DataT& data() const;
    double& operator[] (size_t idx); 
    const double& operator[] (size_t idx) const;
};

} // end namespace tbem

#endif

#ifndef TBEM123123123789798798_DENSE_OPERATOR_H
#define TBEM123123123789798798_DENSE_OPERATOR_H

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

    virtual size_t n_rows() const override;
    virtual size_t n_cols() const override;
    size_t n_elements() const; 

    virtual std::vector<double> apply(const std::vector<double>& x) const override;

    const DataT& data() const;
    double& operator[] (size_t idx); 
    const double& operator[] (size_t idx) const;

    DenseOperator add(const DenseOperator& B);
};

DenseOperator compose_dense_ops(const std::vector<DenseOperator>& ops,
    const std::vector<size_t>& start_rows, const std::vector<size_t>& start_cols,
    const std::vector<double>& multipliers, size_t n_out_rows, size_t n_out_cols);

} // end namespace tbem

#endif

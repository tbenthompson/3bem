#ifndef TBEM_ROW_ZERO_DISTRIBUTOR_H
#define TBEM_ROW_ZERO_DISTRIBUTOR_H

#include "operator.h"

#include <set>

namespace tbem {

struct ConstraintMatrix;
struct DenseOperator;

std::vector<size_t> identify_ignored_dofs(const ConstraintMatrix& cm);

DenseOperator distribute_row_zeros(const DenseOperator& matrix, 
    const ConstraintMatrix& cm); 


struct RowZeroDistributor: public OperatorI
{
    const std::set<size_t> ignored_rows;
    const std::unique_ptr<OperatorI> wrapped_op;

    RowZeroDistributor(const ConstraintMatrix& cm, const OperatorI& op);
    RowZeroDistributor(const std::set<size_t>& cm, const OperatorI& op);

    virtual size_t n_rows() const;
    virtual size_t n_cols() const;
    virtual std::vector<double> apply(const std::vector<double>& x) const;
    virtual std::unique_ptr<OperatorI> clone() const;
};

} // end namespace tbem

#endif

#include "row_zero_distributor.h"
#include "dense_operator.h"
#include "constraint_matrix.h"

#include <set>

namespace tbem {


bool is_ignorable(const std::set<size_t>& ignored, const RearrangedConstraintEQ& c)
{
    for (auto t: c.terms) {
        if (ignored.count(t.dof) == 0) {
            return false; 
        }
    }
    return true;
}

size_t max_constrained_dof(const ConstraintMatrixData& cm) 
{
    size_t max_dof = 0;
    for (auto it = cm.begin(); it != cm.end(); ++it) {
        max_dof = std::max(it->first, max_dof);
    }
    return max_dof;
}

std::set<size_t> _identify_ignored_dofs_set(const ConstraintMatrix& cm)
{
    std::set<size_t> ignored;

    for (size_t dof_index = 0; dof_index <= max_constrained_dof(cm.map); dof_index++) {
        if(!is_constrained(cm.map, dof_index)) {
            continue;
        }
        auto constraint = cm.map.find(dof_index)->second;
        if (is_ignorable(ignored, constraint)) {
            ignored.insert(constraint.constrained_dof);
        }
    }
    return ignored;
}

std::vector<size_t> identify_ignored_dofs(const ConstraintMatrix& cm)
{
    auto ignored = _identify_ignored_dofs_set(cm);
    std::vector<size_t> out(ignored.begin(), ignored.end());
    return out;
}

DenseOperator distribute_row_zeros(const DenseOperator& matrix, 
    const ConstraintMatrix& cm)
{
    auto ignored_rows = _identify_ignored_dofs_set(cm);
    auto n_total_rows = matrix.n_rows() + ignored_rows.size();
    auto next_in_row = 0;
    std::vector<double> out_data(n_total_rows * matrix.n_cols());
    for (size_t i = 0; i < n_total_rows; i++) {
        if (ignored_rows.count(i) > 0) {
            for (size_t j = 0; j < matrix.n_cols(); j++) {
                out_data[i * matrix.n_cols() + j] = 0.0;
            }
        } else {
            for (size_t j = 0; j < matrix.n_cols(); j++) {
                auto value = matrix[next_in_row * matrix.n_cols() + j];
                out_data[i * matrix.n_cols() + j] = value;
            }
            next_in_row++;
        }
    }
    return DenseOperator(n_total_rows, matrix.n_cols(), out_data);
}

RowZeroDistributor::RowZeroDistributor(const ConstraintMatrix& cm, const OperatorI& op):
    ignored_rows(_identify_ignored_dofs_set(cm)),
    wrapped_op(op.clone())
{}

RowZeroDistributor::RowZeroDistributor(const std::set<size_t>& ignored_rows,
        const OperatorI& op):
    ignored_rows(ignored_rows),
    wrapped_op(op.clone())
{}

size_t RowZeroDistributor::n_rows() const 
{
    return wrapped_op->n_rows() + ignored_rows.size();
}

size_t RowZeroDistributor::n_cols() const
{
    return wrapped_op->n_cols();
}

std::vector<double> RowZeroDistributor::apply(const std::vector<double>& x) const
{
    auto intermediate = wrapped_op->apply(x);

    std::vector<double> out(n_rows());
    size_t next_in_row = 0;
    for (size_t i = 0; i < n_rows(); i++) {
        if (ignored_rows.count(i) > 0) {
            out[i] = 0.0;
        } else {
            out[i] = intermediate[next_in_row];
            next_in_row++;
        }
    }

    return out;
}

std::unique_ptr<OperatorI> RowZeroDistributor::clone() const
{
    return std::unique_ptr<OperatorI>(new RowZeroDistributor(
        ignored_rows, *wrapped_op
    ));
}

} // end namespace tbem

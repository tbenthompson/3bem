#ifndef __JJJJJJJJJJJJJJ_PETSC_INTERFACE_H
#define __JJJJJJJJJJJJJJ_PETSC_INTERFACE_H

#include <functional>
#include <memory>
#include <vector>
#include "operator.h"
#include "fwd_vectorx.h"

struct _p_Mat;
typedef struct _p_Mat* Mat;

namespace tbem {

struct MatrixEntry;

struct SparseOperator: public OperatorI {
    SparseOperator(size_t n_rows, size_t n_cols,
        const std::vector<MatrixEntry>& entries);

    virtual size_t n_rows() const;
    virtual size_t n_cols() const;
    virtual VectorX apply(const VectorX& v) const;

    Mat internal_mat;
};
template <typename T>
struct BlockOperator;
typedef BlockOperator<SparseOperator> BlockSparseOperator;

typedef std::function<void(std::vector<double>& x, std::vector<double>& y)> MatVecFnc;

std::vector<double> 
solve_system(const double* rhs, int n_dofs, double tolerance, MatVecFnc fnc);

std::vector<double> solve_system(const VectorX& rhs, double tolerance, MatVecFnc fnc);

} // END namespace tbem
#endif

#ifndef __JJJJJJJJJJJJJJ_PETSC_INTERFACE_H
#define __JJJJJJJJJJJJJJ_PETSC_INTERFACE_H

#include <functional>
#include <memory>
#include <vector>
#include "fwd_vectorx.h"

struct _p_Mat;
typedef struct _p_Mat* Mat;

namespace tbem {

struct MatrixEntry;

struct PETScSparseMatWrapper {
    PETScSparseMatWrapper(size_t n_rows, size_t n_cols,
        const std::vector<MatrixEntry>& entries);

    size_t n_rows() const;
    size_t n_cols() const;
    VectorX mat_vec_prod(const VectorX& v) const;

    Mat internal_mat;
};

typedef std::function<void(std::vector<double>& x, std::vector<double>& y)> MatVecFnc;

std::vector<double> 
solve_system(const double* rhs, int n_dofs, double tolerance, MatVecFnc fnc);

std::vector<double> solve_system(const VectorX& rhs, double tolerance, MatVecFnc fnc);

} // END namespace tbem
#endif

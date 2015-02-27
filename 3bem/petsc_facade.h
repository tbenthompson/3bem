#ifndef __JJJJJJJJJJJJJJ_PETSC_INTERFACE_H
#define __JJJJJJJJJJJJJJ_PETSC_INTERFACE_H

#include <functional>
#include <memory>
#include <vector>
#include "vectorx.h"

struct _p_Mat;
typedef struct _p_Mat* Mat;

namespace tbem {

class PETScMatWrapper {
     Mat internal;
};

typedef std::function<void(std::vector<double>& x, std::vector<double>& y)> MatVecFnc;

std::vector<double> 
solve_system(const double* rhs, int n_dofs, double tolerance, MatVecFnc fnc);

std::vector<double> solve_system(const VectorX& rhs, double tolerance, MatVecFnc fnc);

} // END namespace tbem
#endif

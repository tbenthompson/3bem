#ifndef __JJJJJJJJJJJJJJ_PETSC_INTERFACE_H
#define __JJJJJJJJJJJJJJ_PETSC_INTERFACE_H

#include <functional>
#include <vector>

namespace tbem {

typedef std::function<void(std::vector<double>& x, std::vector<double>& y)> MatVecFnc;

std::vector<double> solve_system(const double* rhs, int n_dofs,
                                 double tolerance, MatVecFnc fnc);

std::vector<double> solve_system(const std::vector<double>& rhs,
                                 double tolerance, MatVecFnc fnc);

} // END namespace tbem
#endif

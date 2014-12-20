#ifndef __JJJJJJJJJJJJJJ_PETSC_INTERFACE_H
#define __JJJJJJJJJJJJJJ_PETSC_INTERFACE_H

#include <functional>
#include <vector>

namespace tbem {

typedef std::function<void(std::vector<double>& x, std::vector<double>& y)> MatVecFnc;

std::vector<double> solve_system(std::vector<double> rhs,
                                 double tolerance,
                                 MatVecFnc fnc);

} // END namespace tbem
#endif

#ifndef __NWQPOISJMNLJSDROIT_BASIS_H
#define __NWQPOISJMNLJSDROIT_BASIS_H
#include "mesh.h"
#include "constraint_matrix.h"
#include "vectorx.h"

namespace tbem {

template <size_t dim, typename Fnc> 
VectorX constrained_interpolate(const Mesh<dim>& mesh,
                                            const Fnc& fnc,
                                            const ConstraintMatrix& matrix) {
    int n_dofs = dim * mesh.facets.size();
    std::vector<double> res(mesh.n_dofs());
#pragma omp parallel for
    for (unsigned int i = 0; i < mesh.facets.size(); i++) {
        for (int d = 0; d < dim; d++) {
            int dof = dim * i + d;
            if(is_constrained(matrix, dof)) {
                continue;
            }
            res[i * dim + d] = fnc(mesh.facets[i][d]);
        }
    }
    return distribute_vector(matrix, res, n_dofs);
}

/* Interpolates a function onto the linear basis defined by the specified
 * mesh.
 */
template <size_t dim, typename Fnc>
VectorX interpolate(const Mesh<dim>& mesh, const Fnc& fnc) {
    return constrained_interpolate<dim,Fnc>(mesh, fnc, ConstraintMatrix{});
}


} //END NAMESPACE tbem

#endif

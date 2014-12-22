#ifndef __NWQPOISJMNLJSDROIT_BASIS_H
#define __NWQPOISJMNLJSDROIT_BASIS_H
#include "mesh.h"
#include "constraint.h"

//TODO: higher order basis!
//TODO: move some of the stuff from numerics here.
//TODO: interpolation with constraints

namespace tbem {

template <int dim, typename Fnc> 
std::vector<double> constrained_interpolate(const Mesh<dim>& mesh,
                                            const Fnc& fnc,
                                            const ConstraintMatrix& c_mat) {
    int n_dofs = dim * mesh.facets.size();
    std::vector<double> res;
    for (unsigned int i = 0; i < mesh.facets.size(); i++) {
        for (int d = 0; d < dim; d++) {
            int dof = dim * i + d;
            auto dof_constraint = c_mat.c_map.find(dof);
            if (dof_constraint != c_mat.c_map.end()) {
                continue;
            }
            res.push_back(fnc(mesh.facets[i].vertices[d]));
        }
    }
    return c_mat.get_all(res, n_dofs);
}

/* Interpolates a function onto the linear basis defined by the specified
 * mesh.
 */
template <int dim, typename Fnc>
std::vector<double> interpolate(const Mesh<dim>& mesh,
                                const Fnc& fnc) {
    return constrained_interpolate<dim,Fnc>(mesh, fnc, ConstraintMatrix{});
}


} //END NAMESPACE tbem

#endif

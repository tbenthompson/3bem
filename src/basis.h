#ifndef __NWQPOISJMNLJSDROIT_BASIS_H
#define __NWQPOISJMNLJSDROIT_BASIS_H
#include "mesh.h"

//TODO: higher order basis!
//TODO: move some of the stuff from numerics here.

template <int dim, typename Fnc>
std::vector<double> interpolate(const Mesh<dim>& mesh,
                                const Fnc& fnc) {
    int n_dofs = dim * mesh.facets.size();
    std::vector<double> res(n_dofs);
    for (unsigned int i = 0; i < mesh.facets.size(); i++) {
        for (int d = 0; d < dim; d++) {
            res[dim * i + d] = fnc(mesh.facets[i].vertices[d]);
        }
    }
    return res;
}

#endif

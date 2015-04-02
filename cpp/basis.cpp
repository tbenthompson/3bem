#include "basis.h"
#include "vectorx.h"
#include "mesh.h"

namespace tbem {

template <size_t dim>
VectorX interpolate(const Mesh<dim>& mesh, 
    const std::function<double(const Vec<double,dim>&)>& fnc) 
{
    std::vector<double> res(mesh.n_dofs());
// #pragma omp parallel for
    for (unsigned int i = 0; i < mesh.facets.size(); i++) {
        for (int d = 0; d < dim; d++) {
            int dof = dim * i + d;
            res[dof] = fnc(mesh.facets[i][d]);
        }
    }
    return res;
}

template VectorX interpolate(const Mesh<2>& mesh, 
    const std::function<double(const Vec<double,2>&)>& fnc);
template VectorX interpolate(const Mesh<3>& mesh, 
    const std::function<double(const Vec<double,3>&)>& fnc);

} //end namespace tbem

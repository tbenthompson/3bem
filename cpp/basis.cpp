#include "basis.h"
#include "mesh.h"

namespace tbem {

template <size_t dim>
std::vector<double> interpolate(const Mesh<dim>& mesh, 
    const std::function<double(const Vec<double,dim>&)>& fnc) 
{
    std::vector<double> res(mesh.n_dofs());
// #pragma omp parallel for
    for (size_t i = 0; i < mesh.facets.size(); i++) {
        for (size_t d = 0; d < dim; d++) {
            auto dof = dim * i + d;
            res[dof] = fnc(mesh.facets[i][d]);
        }
    }
    return res;
}

template std::vector<double> interpolate(const Mesh<2>& mesh, 
    const std::function<double(const Vec<double,2>&)>& fnc);
template std::vector<double> interpolate(const Mesh<3>& mesh, 
    const std::function<double(const Vec<double,3>&)>& fnc);

} //end namespace tbem

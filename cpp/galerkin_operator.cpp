#include "galerkin_operator.h"
#include "mesh.h"
#include "gauss_quad.h"
#include "sparse_operator.h"
#include "numerics.h"

namespace tbem {

template <size_t dim>
SparseOperator make_galerkin_operator(size_t n_components,
    const Mesh<dim>& obs_mesh, const QuadRule<dim-1>& obs_quad)
{
    auto n_obs_dofs = obs_mesh.n_dofs();
    auto n_quad_pts = obs_mesh.n_facets() * obs_quad.size();
    std::vector<MatrixEntry> entries;
    for (size_t obs_idx = 0; obs_idx < obs_mesh.n_facets(); obs_idx++) 
    {
        auto obs_jacobian = facet_jacobian<dim>(obs_mesh.facets[obs_idx]);
        for (size_t obs_q = 0; obs_q < obs_quad.size(); obs_q++) 
        {
            auto quad_idx = obs_idx * obs_quad.size() + obs_q;
            auto basis = linear_basis(obs_quad[obs_q].x_hat);
            auto weight = obs_jacobian * obs_quad[obs_q].w;

            for (size_t obs_basis_idx = 0; obs_basis_idx < dim; obs_basis_idx++) 
            {
                auto obs_dof = dim * obs_idx + obs_basis_idx;
                auto entry_value = basis[obs_basis_idx] * weight;
                for (size_t d = 0; d < n_components; d++) {
                    entries.push_back({
                        d * n_obs_dofs + obs_dof,
                        d * n_quad_pts + quad_idx,
                        entry_value
                    });
                }
            }
        }
    }
    return SparseOperator::csr_from_coo(
        n_components * obs_mesh.n_dofs(),
        n_components * obs_mesh.n_facets() * obs_quad.size(),
        entries
    );
}

template 
SparseOperator make_galerkin_operator(size_t n_components,
    const Mesh<2>& obs_mesh, const QuadRule<1>& obs_quad);
template 
SparseOperator make_galerkin_operator(size_t n_components,
    const Mesh<3>& obs_mesh, const QuadRule<2>& obs_quad);

}//end namespace tbem

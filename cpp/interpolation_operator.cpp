#include "interpolation_operator.h"
#include "mesh.h"
#include "quadrature.h"
#include "sparse_operator.h"
#include "numerics.h"

namespace tbem {

template <size_t dim>
SparseOperator make_interpolation_operator(size_t n_components,
    const Mesh<dim>& src_mesh, const QuadRule<dim-1>& src_quad)
{
    auto n_quadrature_pts = src_mesh.n_facets() * src_quad.size();
    std::vector<MatrixEntry> entries;
    for (size_t idx = 0; idx < src_mesh.facets.size(); idx++) 
    {
        for (size_t q = 0; q < src_quad.size(); q++) 
        {
            auto basis = linear_basis(src_quad[q].x_hat);
            auto interp_dof = idx * src_quad.size() + q;
            for (size_t src_basis_idx = 0; src_basis_idx < dim; src_basis_idx++) 
            {
                auto src_dof = idx * dim + src_basis_idx;
                for (size_t d = 0; d < n_components; d++) {
                    entries.push_back({
                        d * n_quadrature_pts + interp_dof,
                        d * src_mesh.n_dofs() + src_dof,
                        basis[src_basis_idx]
                    });
                }
            }
        }
    }
    return SparseOperator::csr_from_coo(
        n_components * n_quadrature_pts,
        n_components * src_mesh.n_dofs(),
        entries
    );
}

template 
SparseOperator make_interpolation_operator(size_t n_components,
    const Mesh<2>& src_mesh, const QuadRule<1>& src_quad);
template 
SparseOperator make_interpolation_operator(size_t n_components,
    const Mesh<3>& src_mesh, const QuadRule<2>& src_quad);

} // end namespace tbem

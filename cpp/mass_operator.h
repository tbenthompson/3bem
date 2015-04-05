#ifndef __123786172627689192_MASS_OPERATOR_H
#define __123786172627689192_MASS_OPERATOR_H

#include "galerkin_operator.h"

namespace tbem {
template <size_t dim>
struct BlockMassOperator: public BlockOperatorI
{
    const OperatorShape shape;
    const Mesh<dim> obs_mesh;
    const QuadRule<dim-1> obs_quad;
    const BlockGalerkinOperator<dim> galerkin;
    BlockMassOperator(const OperatorShape& shape,
        const Mesh<dim>& obs_mesh, const QuadRule<dim-1>& obs_quad):
        shape(shape), obs_mesh(obs_mesh), obs_quad(obs_quad),
        galerkin{shape, obs_mesh, obs_quad}
    {}

    virtual size_t n_block_rows() const {return shape.n_rows;}
    virtual size_t n_block_cols() const {return shape.n_cols;}
    virtual size_t n_total_rows() const {return shape.n_rows * obs_mesh.n_dofs();} 
    virtual size_t n_total_cols() const {return shape.n_cols * obs_mesh.n_dofs();}

    virtual BlockVectorX apply(const BlockVectorX& x) const {
        auto n_quadrature_pts = obs_mesh.n_facets() * obs_quad.size();
        BlockVectorX interpolated(shape.n_cols, VectorX(n_quadrature_pts, 0.0));
#pragma omp parallel for
        for (size_t obs_idx = 0; obs_idx < obs_mesh.facets.size(); obs_idx++) 
        {
            for (size_t obs_q = 0; obs_q < obs_quad.size(); obs_q++) 
            {
                auto basis = linear_basis(obs_quad[obs_q].x_hat);
                auto interp_dof = obs_idx * obs_quad.size() + obs_q;
                for (size_t src_basis_idx = 0; src_basis_idx < dim; src_basis_idx++) 
                {
                    auto src_dof = obs_idx * dim + src_basis_idx;
                    for (size_t d = 0; d < shape.n_rows; d++) {
                        interpolated[d][interp_dof] += 
                            basis[src_basis_idx] * x[d][src_dof];
                    }
                }
            }
        }
        return galerkin.apply(interpolated);
    }
};

template <size_t dim, size_t R, size_t C>
BlockMassOperator<dim>
mass_operator(const Mesh<dim>& obs_mesh, size_t n_q)
{
    return BlockMassOperator<dim>{{R,C}, obs_mesh, gauss_facet<dim>(n_q)};
}

} // END NAMESPACE tbem
#endif

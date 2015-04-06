#ifndef __LKJLKJK12312111_GALERKIN_OPERATOR_H
#define __LKJLKJK12312111_GALERKIN_OPERATOR_H

#include "numerics.h"
#include "block_operator.h"
#include "mesh.h"
#include "quadrature.h"

namespace tbem {

template <size_t dim>
struct BlockGalerkinOperator: public BlockOperatorI
{
    const OperatorShape shape;
    const Mesh<dim> obs_mesh;
    const QuadRule<dim-1> obs_quad;

    BlockGalerkinOperator(const OperatorShape& shape,
        const Mesh<dim>& obs_mesh, const QuadRule<dim-1>& obs_quad):
        shape(shape), obs_mesh(obs_mesh), obs_quad(obs_quad)
    {}

    virtual size_t n_block_rows() const {return shape.n_rows;}
    virtual size_t n_block_cols() const {return shape.n_cols;}

    virtual size_t n_total_rows() const {
        return shape.n_rows * obs_mesh.n_dofs();
    }

    virtual size_t n_total_cols() const {
        return shape.n_cols * obs_mesh.n_facets() * obs_quad.size();
    }

    virtual BlockVectorX apply(const BlockVectorX& x) const {
        auto n_obs_dofs = obs_mesh.n_dofs();

        BlockVectorX out(shape.n_rows, VectorX(n_obs_dofs, 0.0));
#pragma omp parallel for
        for (size_t obs_idx = 0; obs_idx < obs_mesh.facets.size(); obs_idx++) 
        {
            auto obs_jacobian = facet_jacobian<dim>(obs_mesh.facets[obs_idx]);
            for (size_t obs_q = 0; obs_q < obs_quad.size(); obs_q++) 
            {
                auto interp_dof = obs_idx * obs_quad.size() + obs_q;
                auto basis = linear_basis(obs_quad[obs_q].x_hat);
                auto weight = obs_jacobian * obs_quad[obs_q].w;

                for (size_t obs_basis_idx = 0; obs_basis_idx < dim; obs_basis_idx++) 
                {
                    auto obs_dof = dim * obs_idx + obs_basis_idx;
                    auto entry_value = basis[obs_basis_idx] * weight;
                    for (size_t d = 0; d < shape.n_rows; d++) {
                        out[d][obs_dof] += entry_value * x[d][interp_dof];
                    }
                }
            }
        }

        return out;
             
    }
};

} // end namespace tbem

#endif

#ifndef __KLJLKJLKJNNBMBMNV_INTERPOLATION_OPERATOR_H
#define __KLJLKJLKJNNBMBMNV_INTERPOLATION_OPERATOR_H

#include "block_operator.h"
#include "quadrature.h"

namespace tbem {

template <size_t dim>
struct BlockInterpolationOperator: public BlockOperatorI
{
    const OperatorShape shape;
    const Mesh<dim> mesh;
    const QuadRule<dim-1> quad;

    BlockInterpolationOperator(const OperatorShape& shape,
        const Mesh<dim>& mesh, const QuadRule<dim-1>& quad):
        shape(shape), mesh(mesh), quad(quad)
    {}

    virtual size_t n_block_rows() const {return shape.n_rows;}
    virtual size_t n_block_cols() const {return shape.n_cols;}
    virtual size_t n_total_rows() const {return shape.n_rows * mesh.n_dofs();} 
    virtual size_t n_total_cols() const {return shape.n_cols * mesh.n_dofs();}

    virtual BlockVectorX apply(const BlockVectorX& x) const {
        auto n_quadrature_pts = mesh.n_facets() * quad.size();
        BlockVectorX interpolated(shape.n_cols, VectorX(n_quadrature_pts, 0.0));
#pragma omp parallel for
        for (size_t idx = 0; idx < mesh.facets.size(); idx++) 
        {
            for (size_t q = 0; q < quad.size(); q++) 
            {
                auto basis = linear_basis(quad[q].x_hat);
                auto interp_dof = idx * quad.size() + q;
                for (size_t src_basis_idx = 0; src_basis_idx < dim; src_basis_idx++) 
                {
                    auto src_dof = idx * dim + src_basis_idx;
                    for (size_t d = 0; d < shape.n_rows; d++) {
                        interpolated[d][interp_dof] += 
                            basis[src_basis_idx] * x[d][src_dof];
                    }
                }
            }
        }
        return interpolated;
    }
};

} //end namespace tbem

#endif

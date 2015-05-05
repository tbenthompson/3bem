#ifndef __123786172627689192_MASS_OPERATOR_H
#define __123786172627689192_MASS_OPERATOR_H

#include "galerkin_operator.h"
#include "interpolation_operator.h"

namespace tbem {


template <size_t dim>
struct BlockMassOperator: public BlockOperatorI
{
    const OperatorShape shape;
    const BlockInterpolationOperator<dim> interp;
    const BlockGalerkinOperator<dim> galerkin;

    BlockMassOperator(const OperatorShape& shape,
        const Mesh<dim>& obs_mesh, const QuadRule<dim-1>& obs_quad):
        shape(shape),
        interp(shape, obs_mesh, obs_quad),
        galerkin(shape, obs_mesh, obs_quad)
    {}

    virtual size_t n_block_rows() const {return shape.n_rows;}
    virtual size_t n_block_cols() const {return shape.n_cols;}
    virtual size_t n_total_rows() const {return galerkin.n_total_rows();} 
    virtual size_t n_total_cols() const {return interp.n_total_cols();}

    virtual std::vector<double> apply(const std::vector<double>& x) const {
        return galerkin.apply(interp.apply(x));
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

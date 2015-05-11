#ifndef __123786172627689192_MASS_OPERATOR_H
#define __123786172627689192_MASS_OPERATOR_H

#include "galerkin_operator.h"
#include "interpolation_operator.h"

namespace tbem {


template <size_t dim>
struct MassOperator: public OperatorI
{
    const OperatorShape shape;
    const InterpolationOperator<dim> interp;
    const GalerkinOperator<dim> galerkin;

    MassOperator(const OperatorShape& shape,
        const Mesh<dim>& obs_mesh, const QuadRule<dim-1>& obs_quad):
        shape(shape),
        interp(shape, obs_mesh, obs_quad),
        galerkin(shape, obs_mesh, obs_quad)
    {}

    virtual size_t n_rows() const {return galerkin.n_rows();} 
    virtual size_t n_cols() const {return interp.n_cols();}

    virtual std::vector<double> apply(const std::vector<double>& x) const {
        return galerkin.apply(interp.apply(x));
    }
};

template <size_t dim, size_t R, size_t C>
MassOperator<dim>
mass_operator(const Mesh<dim>& obs_mesh, size_t n_q)
{
    return MassOperator<dim>{{R,C}, obs_mesh, gauss_facet<dim>(n_q)};
}

} // END NAMESPACE tbem
#endif

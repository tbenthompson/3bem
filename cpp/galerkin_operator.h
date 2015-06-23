#ifndef TBEMLKJLKJK12312111_GALERKIN_OPERATOR_H
#define TBEMLKJLKJK12312111_GALERKIN_OPERATOR_H

#include <cstdlib>
#include "quad_rule.h"

namespace tbem {

struct SparseOperator;
struct OperatorShape;
template <size_t dim> struct Mesh;

template <size_t dim>
SparseOperator make_galerkin_operator(size_t n_components,
    const Mesh<dim>& obs_mesh, const QuadRule<dim-1>& obs_quad);

} // end namespace tbem

#endif

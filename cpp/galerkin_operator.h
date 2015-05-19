#ifndef __LKJLKJK12312111_GALERKIN_OPERATOR_H
#define __LKJLKJK12312111_GALERKIN_OPERATOR_H

#include <cstdlib>
#include <vector>

namespace tbem {

struct SparseOperator;
struct OperatorShape;
template <size_t dim> struct Mesh;
template <size_t dim>
struct QuadPt;
template <size_t dim>
using QuadRule = std::vector<QuadPt<dim>>;

template <size_t dim>
SparseOperator make_galerkin_operator(size_t n_components,
    const Mesh<dim>& obs_mesh, const QuadRule<dim-1>& obs_quad);

} // end namespace tbem

#endif

#ifndef TBEMKLJLKJLKJNNBMBMNV_INTERPOLATION_OPERATOR_H
#define TBEMKLJLKJLKJNNBMBMNV_INTERPOLATION_OPERATOR_H

#include <cstdlib>
#include <vector>

namespace tbem {

struct SparseOperator;
template <size_t dim> struct Mesh;
template <size_t dim>
struct QuadPt;
template <size_t dim>
using QuadRule = std::vector<QuadPt<dim>>;

template <size_t dim>
SparseOperator make_interpolation_operator(size_t n_components,
    const Mesh<dim>& src_mesh, const QuadRule<dim-1>& src_quad);

} //end namespace tbem

#endif

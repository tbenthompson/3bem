#ifndef __NWQPOISJMNLJSDROIT_BASIS_H
#define __NWQPOISJMNLJSDROIT_BASIS_H
#include <functional>
#include "fwd_vectorx.h"
#include "vec.h"

namespace tbem {

template <size_t dim>
struct Mesh;

/* Interpolates a function onto the linear basis defined by the specified
 * mesh.
 */
template <size_t dim>
VectorX interpolate(const Mesh<dim>& mesh, 
    const std::function<double(const Vec<double,dim>&)>& fnc);

} //END NAMESPACE tbem

#endif

#ifndef __GGggGGggTTFDSSDf_MESH_H
#define __GGggGGggTTFDSSDf_MESH_H
#include <vector>
#include "function.h"

namespace tbem {

/* There's a pretty interesting correspondence between the interpolation
 * of a field like position between the topological corners of a triangle
 * and the interpolation of a field like displacment, traction, or heat flux
 * between the same corners. I have generalized the Mesh structure into a 
 * Function structure. Imagine that a 2D mesh is defined by 2 functions. The
 * first gives the value along the x-axis for each vertex. The second gives 
 * the value along the y-axis for each vertex.
 */

template <size_t dim>
using Facet = FacetFunction<Vec<double,dim>,dim>;

template <size_t dim>
using Mesh = Function<Vec<double,dim>,dim>;

} //END NAMESPACE tbem

#endif

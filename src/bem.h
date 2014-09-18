#ifndef __BEM_H
#define __BEM_H

#include <array>
#include <vector>

void mesh_to_subsegs(std::vector<std::array<double, 2>>& vertices,
                     std::vector<std::array<int, 2>>& segments,
                     int order) 
{
     
}

// chunks to write:
// refine func
// summation func
// clenshaw curtis quadrature + tanh-sinh quadrature
// richardson extrapolation quadrature
// mappings from reference to real space -- only linear?
// constraints --> boundary conditions, non-singular traction BCs
// which kernels to use for the different boundary integral equations
// the kernels
// evaluate solution on the surface after calculation -- just use the 
//     integral equation again
// interior triangulation and evaluation

#endif

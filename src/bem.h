#ifndef __BEM_H
#define __BEM_H

#include <array>
#include <vector>
#include "numerics.h"
/* [NOTE - What quadrature to use]
 * Source integrals are always smooth because of the numerical limit taken by the Richardson Extrapolation Quadrature method. Gauss quadrature works well for smooth integrals. The order of the quadrature will need to be increased for the nearfield integrals, especially the ones where the limiting parameter is small. For higher order, p, computing a Gauss quadrature rule costs time O(p^2) or requires complex O(p) algorithms. Double exponential quadrature is very simple for any order, so for high order quadrature using a double exponential quadrature formula may be better.  
 * For observation integrals, the outer integral of a galerkin method may be weakly singular at its endpoints. As a result, Gauss quadrature is not a good choice. A simple alternative that can handle endpoint singularities is the Double Exponential (also known as Tanh-Sinh) quadrature rule.
 */

/* A relatively brainless mesh class. There are many operations that modify a mesh to be in a more friendly state. 
 * Better data structures are possible. For example a quad-edge or winged-edge structure has better cache-locality and makes it mucher easier to refine. However, this is simple and sufficient!
 */
class Mesh {
public:
    std::vector<std::array<double, 2>> vertices;
    std::vector<std::array<int, 2>> segments;
};


/* Refine a 2D mesh by simply placing vertices in the middle of the requested edges. The list of requested edges is pass by value and not referenced, because it is internally sorted and thus will be modified if passed by reference. Note that this function entirely reconstructs the mesh -- the original mesh is not modified.
 */
Mesh refine_mesh(const Mesh& m, std::vector<int> refine_these);


/* Data for the evaluation of nearfield integrals. */
class NearEval {
public:
    NearEval(int n_steps, int n_obs);

    void zero_nears(int i);

    static constexpr double initial_dist = 1.0;

    const int n_steps;
    std::vector<std::vector<double>> steps;
    std::vector<QuadratureRule> quad;
    std::vector<double> dist;
};

class SegmentInfo {
public:
    std::array<double, 2> v0;
    std::array<double, 2> v1;
    double length;
    double v0_val;
    double v1_val;
};

/* Compute the minimum distance between the vertices of two segments. This does not strictly compute the distance between arbitrary segments, because if the segments intersect at a non-vertex location, the distance should be 0 and will be non-zero. However, for meshes input to the direct_interact function, there should be no vertices that intersect anywhere except the vertices. */
double appx_segment_distance(std::array<double, 2> v00,
                             std::array<double, 2> v01, 
                             std::array<double, 2> v10, 
                             std::array<double, 2> v11);


inline double laplace_single(double r, double dx, double dy) {
    return (1.0 / (4 * M_PI * r));
}

inline double one(double, double, double) {
    return 1.0;
}

typedef std::function<double (double, double, double)> KernelFnc;
std::vector<double> direct_interact(Mesh& src_mesh,
                                    Mesh& obs_mesh,
                                    QuadratureRule src_quad,
                                    QuadratureRule obs_quad,
                                    KernelFnc kernel,
                                    std::vector<double>& src_strength,
                                    int n_steps); 
// chunks to write:
// mesh cleaning and region determination
// only inputs are vertices, segments and boundary conditions on those segments
// subsegmentation (DONE)
// refine func (DONE) make it cache-friendly
// summation func
// richardson extrapolation quadrature
// mappings from reference to real space
// --linear (DONE)
// --polynomial
// constraints --> boundary conditions, non-singular traction BCs
// which kernels to use for the different boundary integral equations
// the kernels
// evaluate solution on the surface after calculation -- just use the 
//     integral equation again
// interior meshing and evaluation

#endif

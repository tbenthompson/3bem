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
    NearEval(int n_steps);

    static constexpr double initial_dist = 1.0;

    const int n_steps;
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
double appx_segment_distance(std::array<double, 2> pt,
                             std::array<double, 2> v10, 
                             std::array<double, 2> v11);


//TODO: Autogenerate kernels using some python sympy script
inline double laplace_single(double r, double dx, double dy,
                             double nx, double ny) {
    return -1.0 / (2 * M_PI) * log(r);
}

inline double laplace_double(double r, double dx, double dy,
                             double nx, double ny) {
    return (dx * nx + dy * ny) / (2 * M_PI * r * r);
}

inline double one(double, double, double, double, double) {
    return 1.0;
}

typedef std::function<double (double, double, double, double, double)> KernelFnc;
std::vector<double> direct_interact(Mesh& src_mesh,
                                    Mesh& obs_mesh,
                                    QuadratureRule src_quad,
                                    QuadratureRule obs_quad,
                                    KernelFnc kernel,
                                    std::vector<double>& src_strength,
                                    int n_steps); 

double eval_integral_equation(Mesh& src_mesh,
                              QuadratureRule& src_quad,
                              NearEval& near_eval, 
                              std::array<double, 2> obs_pt,
                              std::array<double, 2> obs_normal,
                              KernelFnc& kernel,
                              std::vector<double>& src_strength);

double integral(QuadratureRule& quad_rule,
                KernelFnc& kernel,
                const std::array<double, 2>& src_v0,
                const std::array<double, 2>& src_v1,
                const double src_length,
                const double v0_val, 
                const double v1_val,
                double obs_x, 
                double obs_y);
// chunks to write:
// mesh cleaning and region determination 
// only inputs are pairs/triplets of vertices boundary conditions on those segments
// subsegmentation (NOT YET)
// refine func (DONE)
// summation func (DONE)
// richardson extrapolation quadrature (PARTIAL)
// mappings from reference to real space
// --linear (DONE)
// --polynomial (NOT YET)
// basis 
// --linear (DONE)
// --polynomial (NOT YET)
// constraints --> boundary conditions, non-singular traction BCs
// which kernels to use for the different boundary integral equations
// the kernels
// evaluate solution on the surface after calculation -- just use the 
//     integral equation again
// interior meshing and evaluation
// adaptivity?
// that old list of good things to do!
// transcribe all those old sheets of thoughts
// look into UFL from Fenics as the 

#endif

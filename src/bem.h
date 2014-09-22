#ifndef __BEM_H
#define __BEM_H

#include <array>
#include <vector>
#include "numerics.h"
/* [NOTE - What quadrature to use]
 * Since sources are meant the represent the far-field where kernels are smooth, gauss quadrature is a good choice. Double exponential quadrature would also work, though it will be less efficient for low orders. For higher order, p, computing a Gauss quadrature rule costs time O(p^2) or requires complex O(p) algorithms. Double exponential quadrature is very simple for any order.  
 * For observation points, the outer integral of a galerkin method may be weakly singular at its endpoints. As a result, Gauss quadrature is not a good choice. A simple alternative that can handle endpoint singularities is the Double Exponential (also known as Tanh-Sinh) quadrature rule.
 */

/* A relatively brainless mesh class. There are many operations that modify a mesh to be in a more friendly state. 
 * Better data structures are possible. For example a quad-edge or winged-edge structure
 * has better cache-locality and makes it mucher easier to refine. However, this is simple
 * and sufficient!
 */
class Mesh {
public:
    std::vector<std::array<double, 2>> vertices;
    std::vector<std::array<int, 2>> segments;
};


/* Refine a 2D mesh by simply placing vertices in the middle of the requested edges. The list of requested edges is pass by value and not referenced, because it is internally sorted and thus will be modified if passed by reference. Note that this function entirely reconstructs the mesh -- the original mesh is not modified.
 */
Mesh refine_mesh(const Mesh& m, std::vector<int> refine_these);


class Subsegments {
public:
    std::vector<double> ref_left;
    std::vector<double> ref_center;
    std::vector<double> ref_right;
    std::vector<double> ref_weight;
    std::vector<std::array<double,2>> left;
    std::vector<std::array<double,2>> center;
    std::vector<std::array<double,2>> right;
    std::vector<int> owner;
};

class NearEval {
public:
    NearEval(int n_near_steps, int n_obs);

    void zero_nears(int i);

    static constexpr double initial_dist = 1.0;

    const int n_near_steps;
    std::vector<std::vector<double>> near_steps;
    std::vector<QuadratureRule> near_quad;
    std::vector<double> near_dist;
};

/* Use a quadrature rule to convert a mesh into a set of point sources or observation points with identifying information.  Point sources are useful for treating the BEM problem as an N-body problem.  Observation points are used for the outer integral in a galerkin boundary element method.  
 *
 * The information in each subsegment is provided in order to locate the surrounding region of the point. Near-field evaluations may require higher quadrature order and thus the subsegments will need to be effectively "refined".
 */
Subsegments get_src_obs(Mesh& m, const QuadratureRule& quad_rule);


std::vector<double> direct_interact(Mesh& src_mesh,
                                    Subsegments& src,
                                    Subsegments& obs, 
                                    std::vector<double> src_strength,
                                    int n_near_steps);
// chunks to write:
// mesh cleaning and region determination
// only inputs are vertices, segments and boundary conditions on those segments
// subsegmentation (DONE)
// refine func (DONE)
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

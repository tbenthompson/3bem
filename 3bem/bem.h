#ifndef __AAAAAAAAA_BEM_H
#define __AAAAAAAAA_BEM_H

#include <cassert>
#include <functional>
#include <array>
#include <vector>
#include "vec.h"
#include "numerics.h"
#include "kernels.h"
#include "mesh.h"
#include "adaptive_quad.h"
#include "quadrature.h"


namespace tbem {

template <int dim>
struct Problem {
    const Mesh<dim>& src_mesh;
    const Mesh<dim>& obs_mesh;
    const Kernel<dim>& K;
    const std::vector<double>& src_strength;
};

template <int dim>
class FacetInfo {
public:
    //The responsibility is on the user to maintain the lifetime of the facet.
    const Facet<dim>& face;

    const double area_scale;
    const double length_scale;
    const double jacobian;
    const Vec<double,dim> normal;

    static const double inv_ref_facet_area;
    static FacetInfo<dim> build(const Facet<dim>& facet);
};

template <> const double FacetInfo<3>::inv_ref_facet_area = 2.0;
template <> const double FacetInfo<2>::inv_ref_facet_area = 0.5;

template <int dim>
struct ObsPt {
    static ObsPt<dim> from_face(const Vec<double,dim-1>& ref_loc,
                                const FacetInfo<dim>& obs_face);

    const double len_scale;
    const Vec<double,dim> loc;
    const Vec<double,dim> normal;
    const Vec<double,dim> richardson_dir;
};

template <int dim, typename K>
Vec<double,dim> eval_quad_pt(const Vec<double,dim-1>& x_hat,
                          const K& kernel,
                          const FacetInfo<dim>& face,
                          const Vec<double,dim>& obs_loc,
                          const Vec<double,dim>& obs_n);

std::vector<double> bem_mat_mult(const std::vector<double>& A, 
                                 int n_obs_dofs,
                                 const std::vector<double>& x);

/* Data transfer object for computing integral terms. 
 * Values are stored by reference for efficiency's sake.
 * This means that the responsibility for maintaining their lifetime is on the
 * user.
 */
template <int dim, typename KT>
struct IntegralTerm {
    const QuadStrategy<dim>& qs;
    const KT& k;
    const ObsPt<dim>& obs;
    const FacetInfo<dim>& src_face;
    const double appx_pt_face_dist_squared;
};

// TODO: This template specialization stuff is kind of ugly...
template <int dim>
struct UnitFacetAdaptiveIntegrator {
    template <typename KT>
    Vec<double,dim> operator()(const IntegralTerm<dim,KT>& term, 
                 const Vec<double,dim>& nf_obs_pt);
};


template <int dim, typename KT>
Vec<double,dim> compute_term(const IntegralTerm<dim, KT>& term);

/* It may be that the exact distance to an element is not required,
 * but a low precision approximate distance is useful.
 * This function simply approximates the distance by the distance from
 * the given point to each vertex of the element.
 */
template <int dim>
double appx_face_dist2(const Vec<double,dim>& pt,
                       const std::array<Vec<double,dim>,dim>& vs);

/* Given an unknown function defined in terms of a polynomial basis:
 * u(x) = \sum_i U_i \phi_i(x)
 * Determine the coefficients for the expansion of an integral term in
 * terms of U_i.
 * In other words, this function calculates the vector Q_i where
 * Q_i = \int_{S_{src}} K(x,y) \phi_i(y) dy
 */
template <int dim>
std::vector<double> integral_equation_vector(const Problem<dim>& p,
                                             const QuadStrategy<dim>& qs,
                                             const ObsPt<dim>& obs);

/* Evaluate the integral equation for a specific observation point:
 * \int_{S_{src}} K(x,y) u(y) dy
 * y is given by the ObsPt<dim> obs.
 * The caller provides a quadrature strategy that specifies the order and 
 * tolerance for evaluating the integral equation.
 */
template <int dim>
double eval_integral_equation(const Problem<dim>& p, const QuadStrategy<dim>& qs,
                              const ObsPt<dim>& obs);

/* Given a kernel function and two meshes this function calculates the
 * Galerkin boundary element matrix representing the operator 
 * \int_{S_{obs}} \phi_i(x) \int_{S_{src}} K(x,y) \phi_j(y) dy dx
 * where S_{obs} is the observation mesh, S_{src} is the source mesh,
 * K(x,y) is the kernel function and \phi_i(x) is a basis function.
 */
template <int dim>
std::vector<double> interact_matrix(const Problem<dim>& p,
                                    const QuadStrategy<dim>& qs);

template <int dim>
std::vector<double> direct_interact(const Problem<dim>& p,
                                    const QuadStrategy<dim>& qs);

/* In many integral equations, one of the functions of interest appears
 * outside of an integration, possibly as the result of integrating against
 * a delta function. 
 * In a galerkin boundary element formulation, these free terms look like:
 * \int_S \phi_i(x) u(x) dx.
 * This function calculates such integrals using gaussian quadrature.
 */
template <int dim>
std::vector<double> mass_term(const Problem<dim>& p,
                              const QuadStrategy<dim>& qs);

template <int dim>
double get_len_scale(Mesh<dim>& mesh, int which_face, int q);

#include "bem_impl.h"

} // END NAMESPACE tbem
#endif

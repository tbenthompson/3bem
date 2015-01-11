#ifndef __AAAAAAAAA_BEM_H
#define __AAAAAAAAA_BEM_H

#include <cassert>
#include "vec.h"
#include "numerics.h"
#include "mesh.h"
#include "adaptive_quad.h"
#include "quadrature.h"

namespace tbem {

template <size_t dim>
struct IdentityScalar {
    typedef double OutType;
    typedef double InType;
    typedef double OperatorType;

    OperatorType call_with_no_params() const {
        return 1.0;
    }

    OperatorType operator()(const double& r2, const Vec<double,dim>& delta,
        const Vec<double,dim>& nsrc, const Vec<double,dim>& nobs) const 
    {
        return call_with_no_params();
    }
};

template <size_t dim, size_t n_rows, size_t n_cols>
struct IdentityTensor {
    typedef Vec<double,n_rows> OutType;
    typedef Vec<double,n_cols> InType;
    typedef Vec<Vec<double,n_cols>,n_rows> OperatorType;

    OperatorType call_with_no_params() const {
        auto out = zeros<OperatorType>::make();
        for (size_t i = 0; i < n_rows; i++) {
            for (size_t j = 0; j < n_cols; j++) {
                if (i == j) {
                    out[i][j] = 1.0;
                }
            }
        }
        return out;
    }

    OperatorType operator()(const double& r2, const Vec<double,dim>& delta,
        const Vec<double,dim>& nsrc, const Vec<double,dim>& nobs) const 
    {
        return call_with_no_params();
    }
};

template <size_t dim>
struct ObsPt;
template <size_t dim>
struct FacetInfo;

/* Data transfer object for computing integral terms. 
 * Values are stored by reference for efficiency's sake.
 * This means that the responsibility for maintaining their lifetime is on the
 * user.
 */
template <size_t dim, typename KT>
struct IntegralTerm {
    const QuadStrategy<dim>& qs;
    const KT& k;
    const ObsPt<dim>& obs;
    const FacetInfo<dim>& src_face;
    const double appx_pt_face_dist_squared;
};

template <size_t dim, typename KT>
IntegralTerm<dim,KT> make_integral_term(const QuadStrategy<dim>& qs,
        const KT& k, const ObsPt<dim>& obs, const FacetInfo<dim>& src_face,
        const double appx_pt_face_dist_squared) {
    return IntegralTerm<dim,KT>{qs, k, obs, src_face, appx_pt_face_dist_squared};
}

/* Given observation point information and source face information, the
 * fucntion evaluates the influence of a single source quadrature point
 * on the observation point.
 */
template <size_t dim, typename KT>
Vec<typename KT::OperatorType,dim> eval_point_influence(
    const Vec<double,dim-1>& x_hat, const KT& kernel, const FacetInfo<dim>& face,
    const Vec<double,dim>& obs_loc, const Vec<double,dim>& obs_n);

template <size_t dim>
struct UnitFacetAdaptiveIntegrator {
    template <typename KT>
    Vec<typename KT::OperatorType,dim> 
    operator()(const IntegralTerm<dim,KT>& term, 
               const Vec<double,dim>& nf_obs_pt);
};


/* Compute the full influence of a source facet on an observation point, given
 * a kernel function/Green's function
 */
template <size_t dim, typename KT>
Vec<typename KT::OperatorType,dim> compute_term(const IntegralTerm<dim,KT>& term);

template <size_t dim, typename KT>
struct Problem {
    const Mesh<dim>& src_mesh;
    const Mesh<dim>& obs_mesh;
    const KT& K;
    const std::vector<typename KT::InType>& src_strength;
};

template <size_t dim, typename KT>
Problem<dim,KT> make_problem(const Mesh<dim>& src_mesh,
                             const Mesh<dim>& obs_mesh, const KT& k,
                             const std::vector<typename KT::InType>& src_strength) {
    return {src_mesh, obs_mesh, k, src_strength};
}

template <size_t dim>
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

template <size_t dim>
struct ObsPt {
    static ObsPt<dim> from_face(const Vec<double,dim-1>& ref_loc,
                                const FacetInfo<dim>& obs_face);

    const double len_scale;
    const Vec<double,dim> loc;
    const Vec<double,dim> normal;
    const Vec<double,dim> richardson_dir;
};


/* Given an unknown function defined in terms of a polynomial basis:
 * u(x) = \sum_i U_i \phi_i(x)
 * Determine the coefficients for the expansion of an integral term in
 * terms of U_i.
 * In other words, this function calculates the vector Q_i where
 * Q_i = \int_{S_{src}} K(x,y) \phi_i(y) dy
 */
template <size_t dim, typename KT>
std::vector<typename KT::OperatorType> 
integral_equation_vector(const Problem<dim,KT>& p, const QuadStrategy<dim>& qs,
                         const ObsPt<dim>& obs);

/* Evaluate the integral equation for a specific observation point:
 * \int_{S_{src}} K(x,y) u(y) dy
 * y is given by the ObsPt<dim> obs.
 * The caller provides a quadrature strategy that specifies the order and 
 * tolerance for evaluating the integral equation.
 */
template <size_t dim, typename KT>
double eval_integral_equation(const Problem<dim,KT>& p, const QuadStrategy<dim>& qs,
                              const ObsPt<dim>& obs);

/* Given a kernel function and two meshes this function calculates the
 * Galerkin boundary element matrix representing the operator 
 * \int_{S_{obs}} \phi_i(x) \int_{S_{src}} K(x,y) \phi_j(y) dy dx
 * where S_{obs} is the observation mesh, S_{src} is the source mesh,
 * K(x,y) is the kernel function and \phi_i(x) is a basis function.
 */
template <size_t dim, typename KT>
std::vector<typename KT::OperatorType> 
interact_matrix(const Problem<dim,KT>& p, const QuadStrategy<dim>& qs);

template <size_t dim, typename KT>
std::vector<typename KT::OutType> direct_interact(const Problem<dim,KT>& p,
                                                  const QuadStrategy<dim>& qs);


template <typename KT>
std::vector<typename KT::OutType>
bem_mat_mult(const std::vector<typename KT::OperatorType>& A, const KT& k,
             int n_obs_dofs, const std::vector<typename KT::InType>& x);

/* In many integral equations, one of the functions of interest appears
 * outside of an integration, possibly as the result of integrating against
 * a delta function. 
 * In a galerkin boundary element formulation, these free terms look like:
 * \int_S \phi_i(x) u(x) dx.
 * This function calculates such integrals using gaussian quadrature.
 */
template <size_t dim, typename KT>
std::vector<typename KT::OutType> mass_term(const Problem<dim,KT>& p,
    const QuadStrategy<dim>& qs);

template <size_t dim>
double get_len_scale(Mesh<dim>& mesh, int which_face, int q);

#include "bem_impl.h"

} // END NAMESPACE tbem
#endif

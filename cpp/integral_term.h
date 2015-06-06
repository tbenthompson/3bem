#ifndef TBEMPPPPPPPPKPKPKPKP_INTEGRAL_TERM_H
#define TBEMPPPPPPPPKPKPKPKP_INTEGRAL_TERM_H

#include "vec.h"
#include "kernel.h"
#include "facet_info.h"
#include "quadrature.h"

namespace tbem {

template <size_t dim>
struct ObsPt {
    const Vec<double,dim> loc;
    const Vec<double,dim> normal;
    const Vec<double,dim> richardson_dir;
};

/* Data transfer object for computing integral terms. 
 * Values are stored by reference for efficiency's sake.
 * This means that the responsibility for maintaining their lifetime is on the
 * user.
 */
template <size_t dim, size_t R, size_t C>
struct IntegralTerm {
    const ObsPt<dim>& obs;
    const FacetInfo<dim>& src_face;

    /* Given observation point information and source face information, the
     * fucntion evaluates the influence of a single source quadrature point
     * on the observation point.
     */
    Vec<Vec<Vec<double,C>,R>,dim> eval_point_influence(const Kernel<dim,R,C>& k,
        const Vec<double,dim-1>& x_hat, const Vec<double,dim>& moved_obs_loc) const; 

    Vec<Vec<Vec<double,C>,R>,dim> eval_point_influence(const Kernel<dim,R,C>& k,
        const Vec<double,dim-1>& x_hat) const; 
};

enum class FarNearType {
    Singular,
    Nearfield,
    Farfield 
};

template <size_t dim>
struct NearestPoint {
    const Vec<double,dim-1> ref_pt;
    const Vec<double,dim> pt;
    const double distance;
    const FarNearType type;
};

template <size_t dim>
struct FarNearLogic {
    double far_threshold;
    double singular_threshold;
    
    NearestPoint<dim> decide(const Vec<double,dim>& pt, const FacetInfo<dim>& facet);
};

template <size_t dim,size_t R, size_t C>
struct NearfieldIntegratorI {
    virtual Vec<Vec<Vec<double,C>,R>,dim> 
    compute_nearfield(const Kernel<dim,R,C>&, const IntegralTerm<dim,R,C>&,
        const NearestPoint<dim>&) const = 0;
};

template <size_t dim,size_t R, size_t C>
struct IntegrationStrategy {
    std::shared_ptr<Kernel<dim,R,C>> K;
    QuadRule<dim-1> src_far_quad;
    QuadRule<dim-1> obs_quad;
    std::vector<double> singular_steps;
    double far_threshold;
    std::unique_ptr<NearfieldIntegratorI<dim,R,C>> nearfield_integrator;

    Vec<Vec<Vec<double,C>,R>,dim>
    compute_singular(const IntegralTerm<dim,R,C>&, const NearestPoint<dim>&) const;

    Vec<Vec<Vec<double,C>,R>,dim> 
    compute_farfield(const IntegralTerm<dim,R,C>&, const NearestPoint<dim>&) const;

    /* Compute the full influence of a source facet on an observation point, given
     * a kernel function/Green's function
     */
    Vec<Vec<Vec<double,C>,R>,dim> compute_term(const IntegralTerm<dim,R,C>& term,
        const NearestPoint<dim>& nearest_pt) const;
};

template <size_t dim,size_t R, size_t C>
struct AdaptiveIntegrator: public NearfieldIntegratorI<dim,R,C> {
    double near_tol;

    AdaptiveIntegrator(double near_tol):
        near_tol(near_tol)
    {}

    virtual Vec<Vec<Vec<double,C>,R>,dim> 
    compute_nearfield(const Kernel<dim,R,C>& K, const IntegralTerm<dim,R,C>&,
        const NearestPoint<dim>&) const;
};

template <size_t dim,size_t R, size_t C>
IntegrationStrategy<dim,R,C> make_adaptive_integrator(
    const QuadStrategy<dim>& qs, const Kernel<dim,R,C>& K) 
{
    return {
        K.clone(), 
        qs.src_far_quad,
        qs.obs_quad,
        qs.singular_steps,
        qs.far_threshold,
        std::unique_ptr<AdaptiveIntegrator<dim,R,C>>(
            new AdaptiveIntegrator<dim,R,C>(qs.near_tol)
        )
    };
}

template <size_t dim,size_t R, size_t C>
struct SinhIntegrator: public NearfieldIntegratorI<dim,R,C> {
    size_t sinh_order;

    SinhIntegrator(size_t sinh_order):
        sinh_order(sinh_order)
    {}

    virtual Vec<Vec<Vec<double,C>,R>,dim> 
    compute_nearfield(const Kernel<dim,R,C>& K, const IntegralTerm<dim,R,C>&,
        const NearestPoint<dim>&) const;
};

template <size_t dim,size_t R, size_t C>
IntegrationStrategy<dim,R,C> make_sinh_integrator(size_t sinh_order,
    const QuadStrategy<dim>& qs, const Kernel<dim,R,C>& K) 
{
    return {
        K.clone(), 
        qs.src_far_quad,
        qs.obs_quad,
        qs.singular_steps,
        qs.far_threshold,
        std::unique_ptr<SinhIntegrator<dim,R,C>>(
            new SinhIntegrator<dim,R,C>(sinh_order)
        )
    };
}


} //end namespace tbem

#endif

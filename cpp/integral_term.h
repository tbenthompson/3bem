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
struct IntegrationMethodI {
    virtual Vec<Vec<Vec<double,C>,R>,dim>
    compute_singular(const IntegralTerm<dim,R,C>&, const NearestPoint<dim>&) const = 0;

    virtual Vec<Vec<Vec<double,C>,R>,dim> 
    compute_nearfield(const IntegralTerm<dim,R,C>&, const NearestPoint<dim>&) const = 0;

    virtual Vec<Vec<Vec<double,C>,R>,dim> 
    compute_farfield(const IntegralTerm<dim,R,C>&, const NearestPoint<dim>&) const = 0;

    virtual double far_threshold() const = 0;
    virtual const Kernel<dim,R,C>& get_kernel() const = 0;
    virtual QuadRule<dim-1> get_src_quad() const = 0;
    virtual QuadRule<dim-1> get_obs_quad() const = 0;

    /* Compute the full influence of a source facet on an observation point, given
     * a kernel function/Green's function
     */
    Vec<Vec<Vec<double,C>,R>,dim> compute_term(const IntegralTerm<dim,R,C>& term,
        const NearestPoint<dim>& nearest_pt) const;
};

template <size_t dim,size_t R, size_t C>
struct AdaptiveIntegrationMethod: public IntegrationMethodI<dim,R,C> {
    QuadStrategy<dim> qs;
    std::shared_ptr<Kernel<dim,R,C>> K;

    AdaptiveIntegrationMethod(const QuadStrategy<dim>& qs, const Kernel<dim,R,C>& K):
        qs(qs), K(K.clone())
    {}

    virtual Vec<Vec<Vec<double,C>,R>,dim>
    compute_singular(const IntegralTerm<dim,R,C>&, const NearestPoint<dim>&) const;

    virtual Vec<Vec<Vec<double,C>,R>,dim> 
    compute_nearfield(const IntegralTerm<dim,R,C>&, const NearestPoint<dim>&) const;

    virtual Vec<Vec<Vec<double,C>,R>,dim> 
    compute_farfield(const IntegralTerm<dim,R,C>&, const NearestPoint<dim>&) const;

    virtual double far_threshold() const {return qs.far_threshold;}
    virtual const Kernel<dim,R,C>& get_kernel() const {return *K;}
    virtual QuadRule<dim-1> get_src_quad() const {return qs.src_far_quad;}
    virtual QuadRule<dim-1> get_obs_quad() const {return qs.obs_quad;}
};

template <size_t dim,size_t R, size_t C>
AdaptiveIntegrationMethod<dim,R,C> make_adaptive_integration_mthd(
    const QuadStrategy<dim>& qs, const Kernel<dim,R,C>& K) 
{
    return AdaptiveIntegrationMethod<dim,R,C>(qs, K);
}

template <size_t dim,size_t R, size_t C>
struct SinhIntegrationMethod: public IntegrationMethodI<dim,R,C> {
    size_t sinh_order;
    QuadStrategy<dim> qs;
    AdaptiveIntegrationMethod<dim,R,C> adaptive;
    std::shared_ptr<Kernel<dim,R,C>> K;

    SinhIntegrationMethod(size_t sinh_order, const QuadStrategy<dim>& qs, 
        const Kernel<dim,R,C>& K):
        sinh_order(sinh_order), qs(qs), adaptive(qs, K), K(K.clone())
    {}

    virtual Vec<Vec<Vec<double,C>,R>,dim>
    compute_singular(const IntegralTerm<dim,R,C>&, const NearestPoint<dim>&) const;

    virtual Vec<Vec<Vec<double,C>,R>,dim> 
    compute_nearfield(const IntegralTerm<dim,R,C>&, const NearestPoint<dim>&) const;

    virtual Vec<Vec<Vec<double,C>,R>,dim> 
    compute_farfield(const IntegralTerm<dim,R,C>&, const NearestPoint<dim>&) const;

    virtual double far_threshold() const {return qs.far_threshold;}
    virtual const Kernel<dim,R,C>& get_kernel() const {return *K;}
    virtual QuadRule<dim-1> get_src_quad() const {return qs.src_far_quad;}
    virtual QuadRule<dim-1> get_obs_quad() const {return qs.obs_quad;}
};

template <size_t dim,size_t R, size_t C>
SinhIntegrationMethod<dim,R,C> make_sinh_integration_mthd(size_t sinh_order,
    const QuadStrategy<dim>& qs, const Kernel<dim,R,C>& K) 
{
    return SinhIntegrationMethod<dim,R,C>(sinh_order, qs, K);
}


} //end namespace tbem

#endif

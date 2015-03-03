#ifndef __PPPPPPPPKPKPKPKP_INTEGRAL_TERM_H
#define __PPPPPPPPKPKPKPKP_INTEGRAL_TERM_H

#include <cassert>
#include <cstdlib>
#include <iomanip> 
#include "vec.h"
#include "quadrature.h"
#include "adaptive_quad.h"
#include "numerics.h"
#include "obs_pt.h"
#include "facet_info.h"

namespace tbem {

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
};

template <size_t dim, typename KT>
IntegralTerm<dim,KT> make_integral_term(const QuadStrategy<dim>& qs,
        const KT& k, const ObsPt<dim>& obs, const FacetInfo<dim>& src_face)
{
    return IntegralTerm<dim,KT>{qs, k, obs, src_face};
}

/* Given observation point information and source face information, the
 * fucntion evaluates the influence of a single source quadrature point
 * on the observation point.
 */
template <size_t dim, typename KT>
Vec<typename KT::OperatorType,dim> eval_point_influence(const Vec<double,dim-1>& x_hat,
                                  const KT& kernel,
                                  const FacetInfo<dim>& face,
                                  const Vec<double,dim>& obs_loc,
                                  const Vec<double,dim>& obs_n) {
    const auto src_pt = ref_to_real(x_hat, face.face);
    const auto d = src_pt - obs_loc;
    const auto r2 = dot_product(d, d);
    const auto kernel_val = kernel(r2, d, face.normal, obs_n);
    return outer_product(linear_basis(x_hat), kernel_val * face.jacobian);
}

template <size_t dim>
struct UnitFacetAdaptiveIntegrator {
    template <typename KT>
    Vec<typename KT::OperatorType,dim> 
    operator()(const IntegralTerm<dim,KT>& term, 
               const Vec<double,dim>& nf_obs_pt);
};


template <>
struct UnitFacetAdaptiveIntegrator<2> {
    template <typename KT>
    Vec<typename KT::OperatorType,2> operator()(const IntegralTerm<2,KT>& term, 
                 const Vec<double,2>& nf_obs_pt) {
        return adaptive_integrate<Vec<typename KT::OperatorType,2>>(
            [&] (double x_hat) {
                Vec<double,1> q_pt = {x_hat};
                return eval_point_influence<2>(q_pt, term.k, term.src_face,
                                       nf_obs_pt, term.obs.normal);
            }, -1.0, 1.0, term.qs.near_tol);
    }
};

template <>
struct UnitFacetAdaptiveIntegrator<3> {
    template <typename KT>
    Vec<typename KT::OperatorType,3> operator()(const IntegralTerm<3,KT>& term, 
                                                const Vec<double,3>& nf_obs_pt) {
        return adaptive_integrate<Vec<typename KT::OperatorType,3>>(
            [&] (double x_hat) {
                if (x_hat == 1.0) {
                    return zeros<Vec<typename KT::OperatorType,3>>::make();
                }
                return adaptive_integrate<Vec<typename KT::OperatorType,3>>(
                    [&] (double y_hat) {
                        Vec<double,2> q_pt = {x_hat, y_hat};
                        return eval_point_influence<3>(q_pt, term.k, term.src_face,
                                                 nf_obs_pt, term.obs.normal);
                    }, 0.0, 1 - x_hat, term.qs.near_tol);
            }, 0.0, 1.0, term.qs.near_tol);
    }
};

template <size_t dim, typename KT> 
Vec<typename KT::OperatorType,dim> 
compute_adaptively(const IntegralTerm<dim, KT>& term,
                   const Vec<double,dim>& nf_obs_pt) {
    UnitFacetAdaptiveIntegrator<dim> integrator;
    return integrator(term, nf_obs_pt);
}

template <size_t dim, typename KT>
Vec<double,dim> get_step_loc(const IntegralTerm<dim,KT>& term, int step_idx) {
    const double safe_dist_ratio = 5.0;
    double step_distance = safe_dist_ratio * term.obs.len_scale * 
                           term.qs.singular_steps[step_idx];
    return term.obs.loc + step_distance * term.obs.richardson_dir;
}

template <typename T>
T richardson_limit(const std::vector<T>& values) {
    assert(values.size() > 1);

    int n_steps = values.size();
    std::vector<T> last_level = values;
    std::vector<T> this_level;
    for (int m = 1; m < n_steps; m++) {
        this_level.resize(n_steps - m);
        for (int i = 0; i < n_steps - m; i++) {
            const double mult = pow(2, m);
            const double factor = 1.0 / (mult - 1.0);
            T low = last_level[i];
            T high = last_level[i + 1];
            T moreacc = factor * (mult * high - low);
            this_level[i] = moreacc;
        }
        last_level = this_level;
    }
    return this_level[0];
}

template <size_t dim, typename KT> 
Vec<typename KT::OperatorType,dim> compute_as_limit(const IntegralTerm<dim, KT>& term) {
    std::vector<Vec<typename KT::OperatorType,dim>> 
        near_steps(term.qs.n_singular_steps);

    for (int step_idx = 0; step_idx < term.qs.n_singular_steps; step_idx++) {
        auto step_loc = get_step_loc(term, step_idx);
        near_steps[step_idx] = compute_adaptively<dim>(term, step_loc);
    }

    return richardson_limit(near_steps);
}

template <size_t dim, typename KT>
Vec<typename KT::OperatorType,dim> compute_far_term(const IntegralTerm<dim, KT>& term) {
    auto integrals = zeros<Vec<typename KT::OperatorType,dim>>::make();
    for (size_t i = 0; i < term.qs.src_far_quad.size(); i++) {
        integrals += eval_point_influence<dim>(
                        term.qs.src_far_quad[i].x_hat,
                        term.k, term.src_face,
                        term.obs.loc, term.obs.normal
                    ) * term.qs.src_far_quad[i].w;
    }
    return integrals;
}

/* It may be that the exact distance to an element is not required,
 * but a low precision approximate distance is useful.
 * This function simply approximates the distance by the distance from
 * the given point to each vertex of the element.
 * A better approximation to the distance to a face might include
 * the centroid (see appx_face_dist2)
 */
template <size_t dim>
double appx_face_dist2(const Vec<double,dim>& pt,
                       const std::array<Vec<double,dim>,dim>& vs) {
    double res = dist2(pt, vs[0]);
    for (int d = 1; d < dim; d++) {
        res = std::min(res, dist2(pt, vs[d])); 
    }
    return res;
}

enum class FarNearType {
    Singular,
    Nearfield,
    Farfield 
};

struct FarNearLogic {
    double far_threshold;
    double singular_threshold;
    
    template <size_t dim>
    FarNearType decide(const ObsPt<dim>& obs, const FacetInfo<dim>& facet) {
        auto appx_dist2 = appx_face_dist2(obs.loc, facet.face);
        bool nearfield = appx_dist2 < pow(far_threshold, 2) * facet.area_scale;
        if (nearfield) {
            bool singular = appx_dist2 < singular_threshold * facet.area_scale;
            if (singular) { 
                return FarNearType::Singular;
            } else {
                return FarNearType::Nearfield;
            }
        } else {
            return FarNearType::Farfield;
        }
    }
};

/* Compute the full influence of a source facet on an observation point, given
 * a kernel function/Green's function
 */
template <size_t dim, typename KT>
Vec<typename KT::OperatorType,dim> compute_term(const IntegralTerm<dim,KT>& term) {
    FarNearLogic far_near_logic{term.qs.far_threshold, 3.0};
    switch (far_near_logic.decide(term.obs, term.src_face)) {
        case FarNearType::Singular:
            return compute_as_limit(term);
            break;
        case FarNearType::Nearfield:
            return compute_adaptively<dim>(term, term.obs.loc);
            break;
        case FarNearType::Farfield:
            return compute_far_term(term);
            break;
    }
    throw std::exception();
}



} //end namespace tbem

#endif

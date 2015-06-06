#include <assert.h>
#include "integral_term.h"
#include "richardson.h"
#include "adaptive_quad.h"
#include "numerics.h"
#include "closest_pt_on_facet.h"

namespace tbem {

template <size_t dim>
NearestPoint<dim> FarNearLogic<dim>::decide(const Vec<double,dim>& pt,
    const FacetInfo<dim>& facet) 
{
    auto near_ref_pt = closest_pt_facet(pt, facet.face);
    auto near_pt = ref_to_real(near_ref_pt, facet.face);
    auto exact_dist2 = dist2(near_pt, pt);
    auto appx_dist2 = exact_dist2;//appx_face_dist2(obs.loc, facet.face);
    bool nearfield = appx_dist2 < pow(far_threshold, 2) * facet.area_scale;
    if (nearfield) {
        bool singular = appx_dist2 < singular_threshold * facet.area_scale;
        if (singular) { 
            return {
                near_ref_pt, near_pt, std::sqrt(exact_dist2), FarNearType::Singular
            };
        } else {
            return {
                near_ref_pt, near_pt, std::sqrt(exact_dist2), FarNearType::Nearfield
            };
        }
    } else {
        return {near_ref_pt, near_pt, std::sqrt(exact_dist2), FarNearType::Farfield};
    }
}

template struct FarNearLogic<2>;
template struct FarNearLogic<3>;

template <size_t dim, size_t R, size_t C>
Vec<Vec<Vec<double,C>,R>,dim> IntegralTerm<dim,R,C>::eval_point_influence(
    const Kernel<dim,R,C>& k, const Vec<double,dim-1>& x_hat,
    const Vec<double,dim>& moved_obs_loc) const 
{
    const auto src_pt = ref_to_real(x_hat, src_face.face);
    auto kernel_val = k(moved_obs_loc, src_pt, obs.normal, src_face.normal);
    return outer_product(linear_basis(x_hat), kernel_val * src_face.jacobian);
}

template <size_t dim, size_t R, size_t C>
Vec<Vec<Vec<double,C>,R>,dim> IntegralTerm<dim,R,C>::eval_point_influence(
    const Kernel<dim,R,C>& k, const Vec<double,dim-1>& x_hat) const 
{
    return eval_point_influence(k, x_hat, obs.loc);
}

template struct IntegralTerm<2,1,1>;
template struct IntegralTerm<2,2,2>;
template struct IntegralTerm<3,1,1>;
template struct IntegralTerm<3,3,3>;

template <size_t dim, size_t R, size_t C>
Vec<Vec<Vec<double,C>,R>,dim> IntegrationStrategy<dim,R,C>::compute_term(
    const IntegralTerm<dim,R,C>& term, const NearestPoint<dim>& nearest_pt) const
{
    switch (nearest_pt.type) {
        case FarNearType::Singular:
            return compute_singular(term, nearest_pt);
            break;
        case FarNearType::Nearfield:
            return nearfield_integrator->compute_nearfield(*K, term, nearest_pt);
            break;
        case FarNearType::Farfield:
            return compute_farfield(term, nearest_pt);
            break;
    }
    throw std::exception();
}

template <size_t dim, size_t R, size_t C>
Vec<Vec<Vec<double,C>,R>,dim> 
IntegrationStrategy<dim,R,C>::compute_singular(const IntegralTerm<dim,R,C>& term,
    const NearestPoint<dim>& nearest_pt) const 
{
    (void)nearest_pt;
    std::vector<Vec<Vec<Vec<double,C>,R>,dim>> steps(singular_steps.size());

    for (int step_idx = 0; step_idx < singular_steps.size(); step_idx++) {
        auto step_loc = get_step_loc(term.obs, singular_steps[step_idx]);
        auto shifted_near_pt = FarNearLogic<dim>{far_threshold, 1.0}
            .decide(step_loc, term.src_face);
        steps[step_idx] = nearfield_integrator->compute_nearfield(
            *K,
            {
                {step_loc, term.obs.normal, term.obs.richardson_dir},
                term.src_face
            },
            shifted_near_pt
        );
    }

    return richardson_limit(2, steps);
}

template <size_t dim, size_t R, size_t C>
Vec<Vec<Vec<double,C>,R>,dim> 
IntegrationStrategy<dim,R,C>::compute_farfield(const IntegralTerm<dim,R,C>& term,
    const NearestPoint<dim>& nearest_pt) const 
{
    (void)nearest_pt;
    auto integrals = zeros<Vec<Vec<Vec<double,C>,R>,dim>>::make();
    for (size_t i = 0; i < src_far_quad.size(); i++) {
        integrals += term.eval_point_influence(*K, src_far_quad[i].x_hat) *
                     src_far_quad[i].w;
    }
    return integrals;
}

template struct IntegrationStrategy<2,1,1>;
template struct IntegrationStrategy<2,2,2>;
template struct IntegrationStrategy<3,1,1>;
template struct IntegrationStrategy<3,3,3>;

template <size_t dim>
struct UnitFacetAdaptiveIntegrator {
    template <size_t R, size_t C>
    Vec<Vec<Vec<double,C>,R>,2> operator()(const Kernel<dim,R,C>& k, double tolerance, 
        const IntegralTerm<dim,R,C>& term, const Vec<double,dim>& nf_obs_pt);
};

template <>
struct UnitFacetAdaptiveIntegrator<2> {
    template <size_t R, size_t C>
    Vec<Vec<Vec<double,C>,R>,2> operator()(const Kernel<2,R,C>& k, double tolerance, 
        const IntegralTerm<2,R,C>& term, const Vec<double,2>& nf_obs_pt) 
    {
        return adaptive_integrate<Vec<Vec<Vec<double,C>,R>,2>>(
            [&] (double x_hat) {
                Vec<double,1> q_pt = {x_hat};
                return term.eval_point_influence(k, q_pt, nf_obs_pt);
            }, -1.0, 1.0, tolerance);
    }
};

template <>
struct UnitFacetAdaptiveIntegrator<3> {
    template <size_t R, size_t C>
    Vec<Vec<Vec<double,C>,R>,3> operator()(const Kernel<3,R,C>& k, double tolerance, 
        const IntegralTerm<3,R,C>& term, const Vec<double,3>& nf_obs_pt) 
    {
        return adaptive_integrate<Vec<Vec<Vec<double,C>,R>,3>>(
            [&] (double x_hat) {
                if (x_hat == 1.0) {
                    return zeros<Vec<Vec<Vec<double,C>,R>,3>>::make();
                }
                return adaptive_integrate<Vec<Vec<Vec<double,C>,R>,3>>(
                    [&] (double y_hat) {
                        Vec<double,2> q_pt = {x_hat, y_hat};
                        return term.eval_point_influence(k, q_pt, nf_obs_pt);
                    }, 0.0, 1 - x_hat, tolerance);
            }, 0.0, 1.0, tolerance);
    }
};

template <size_t dim>
Vec<double,dim> get_step_loc(const ObsPt<dim>& obs, double step_size) {
    const double safe_dist_ratio = 1.0;
    double step_distance = safe_dist_ratio * step_size;
    return obs.loc + step_distance * obs.richardson_dir;
}


template <size_t dim, size_t R, size_t C>
Vec<Vec<Vec<double,C>,R>,dim> 
AdaptiveIntegrator<dim,R,C>::compute_nearfield(const Kernel<dim,R,C>& K, 
    const IntegralTerm<dim,R,C>& term,
    const NearestPoint<dim>& nearest_pt) const 
{
    (void)nearest_pt;
    return UnitFacetAdaptiveIntegrator<dim>()(K, near_tol, term, term.obs.loc);
}

template struct AdaptiveIntegrator<2,1,1>;
template struct AdaptiveIntegrator<2,2,2>;
template struct AdaptiveIntegrator<3,1,1>;
template struct AdaptiveIntegrator<3,3,3>;

template <size_t dim>
QuadRule<dim-1> choose_sinh_quad(size_t sinh_order,
    double S, double l, Vec<double,dim-1> singular_pt);

template <>
inline QuadRule<1> choose_sinh_quad<2>(size_t sinh_order, 
    double S, double l, Vec<double,1> singular_pt) 
{
    assert(l > 0);
    assert(S > 0);
    if ((l / S) > 1) {
        return gauss(10);
    }
    else {
        size_t n = static_cast<size_t>(sinh_order * (1 + std::log(S / l)));
        return sinh_transform(gauss(n), singular_pt[0], l, false);
    }
}

template <>
inline QuadRule<2> choose_sinh_quad<3>(size_t sinh_order,
    double S, double l, Vec<double,2> singular_pt) 
{
    assert(l > 0);
    assert(S > 0);
    if ((l / S) > 0.5) {
        return tri_gauss(10); 
    }
    else {
        size_t n = static_cast<size_t>(sinh_order * (1 + std::log(S / l)));
        auto q = sinh_sigmoidal_transform(gauss(2.0 * n), gauss(n), 
            singular_pt[0], singular_pt[1], l, false);
        return q;
    }
}

template <size_t dim, size_t R, size_t C>
Vec<Vec<Vec<double,C>,R>,dim> sinh_integrate(size_t sinh_order, 
        const Kernel<dim,R,C>& k, const IntegralTerm<dim,R,C>& term,
        const NearestPoint<dim>& near_pt, const Vec<double,dim>& nf_obs_pt) 
{
    assert(near_pt.distance > 0);
    auto S = term.src_face.length_scale;
    auto l = near_pt.distance;
    auto q = choose_sinh_quad<dim>(sinh_order, S, l, near_pt.ref_pt);
    auto integrals = zeros<Vec<Vec<Vec<double,C>,R>,dim>>::make();
    for (size_t i = 0; i < q.size(); i++) {
        integrals += term.eval_point_influence(k, q[i].x_hat, nf_obs_pt) * q[i].w;
    }
    return integrals;
}

template <size_t dim, size_t R, size_t C>
Vec<Vec<Vec<double,C>,R>,dim> 
SinhIntegrator<dim,R,C>::compute_nearfield(const Kernel<dim,R,C>& K, 
    const IntegralTerm<dim,R,C>& term,
    const NearestPoint<dim>& nearest_pt) const 
{
    return sinh_integrate(sinh_order, K, term, nearest_pt, term.obs.loc);
}

template struct SinhIntegrator<2,1,1>;
template struct SinhIntegrator<2,2,2>;
template struct SinhIntegrator<3,1,1>;
template struct SinhIntegrator<3,3,3>;

} //end namespace tbem

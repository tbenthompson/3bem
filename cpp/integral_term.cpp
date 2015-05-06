#include <assert.h>
#include "integral_term.h"
#include "richardson.h"
#include "adaptive_quad.h"
#include "numerics.h"
#include "closest_pt.h"

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
Vec<Vec<Vec<double,C>,R>,dim> IntegrationMethodI<dim,R,C>::compute_term(
    const IntegralTerm<dim,R,C>& term, const NearestPoint<dim>& nearest_pt) const
{
    switch (nearest_pt.type) {
        case FarNearType::Singular:
            return compute_singular(term, nearest_pt);
            break;
        case FarNearType::Nearfield:
            return compute_nearfield(term, nearest_pt);
            break;
        case FarNearType::Farfield:
            return compute_farfield(term, nearest_pt);
            break;
    }
    throw std::exception();
}

template struct IntegrationMethodI<2,1,1>;
template struct IntegrationMethodI<2,2,2>;
template struct IntegrationMethodI<3,1,1>;
template struct IntegrationMethodI<3,3,3>;

template <size_t dim>
struct UnitFacetAdaptiveIntegrator {
    template <size_t R, size_t C>
    Vec<Vec<Vec<double,C>,R>,2> operator()(const Kernel<dim,R,C>& k, double tolerance, 
        const IntegralTerm<dim,R,C>& term, const NearestPoint<dim>& near_pt,
        const Vec<double,dim>& nf_obs_pt);
};

template <>
struct UnitFacetAdaptiveIntegrator<2> {
    template <size_t R, size_t C>
    Vec<Vec<Vec<double,C>,R>,2> operator()(const Kernel<2,R,C>& k, double tolerance, 
        const IntegralTerm<2,R,C>& term, const NearestPoint<2>& near_pt,
        const Vec<double,2>& nf_obs_pt) 
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
        const IntegralTerm<3,R,C>& term, const NearestPoint<3>& near_pt,
        const Vec<double,3>& nf_obs_pt) 
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
    const double safe_dist_ratio = 5.0;
    double step_distance = safe_dist_ratio * obs.len_scale * step_size;
    return obs.loc + step_distance * obs.richardson_dir;
}

template <size_t dim, size_t R, size_t C>
Vec<Vec<Vec<double,C>,R>,dim> 
AdaptiveIntegrationMethod<dim,R,C>::compute_singular(const IntegralTerm<dim,R,C>& term,
    const NearestPoint<dim>& nearest_pt) const 
{
    std::vector<Vec<Vec<Vec<double,C>,R>,dim>> near_steps(qs.n_singular_steps);

    for (int step_idx = 0; step_idx < qs.n_singular_steps; step_idx++) {
        auto step_loc = get_step_loc(term.obs, qs.singular_steps[step_idx]);
        auto shifted_near_pt = FarNearLogic<dim>{qs.far_threshold, 1.0}
            .decide(step_loc, term.src_face);
        UnitFacetAdaptiveIntegrator<dim> integrator;
        near_steps[step_idx] = integrator(K, qs.near_tol, term, shifted_near_pt, step_loc);
    }

    return richardson_limit(2, near_steps);
}

template <size_t dim, size_t R, size_t C>
Vec<Vec<Vec<double,C>,R>,dim> 
AdaptiveIntegrationMethod<dim,R,C>::compute_nearfield(const IntegralTerm<dim,R,C>& term,
    const NearestPoint<dim>& nearest_pt) const 
{
    UnitFacetAdaptiveIntegrator<dim> integrator;
    return integrator(K, qs.near_tol, term, nearest_pt, term.obs.loc);
}

template <size_t dim, size_t R, size_t C>
Vec<Vec<Vec<double,C>,R>,dim> 
AdaptiveIntegrationMethod<dim,R,C>::compute_farfield(const IntegralTerm<dim,R,C>& term,
    const NearestPoint<dim>& nearest_pt) const 
{
    auto integrals = zeros<Vec<Vec<Vec<double,C>,R>,dim>>::make();
    for (size_t i = 0; i < qs.src_far_quad.size(); i++) {
        integrals += term.eval_point_influence(K, qs.src_far_quad[i].x_hat) *
                     qs.src_far_quad[i].w;
    }
    return integrals;
}
template struct AdaptiveIntegrationMethod<2,1,1>;
template struct AdaptiveIntegrationMethod<2,2,2>;
template struct AdaptiveIntegrationMethod<3,1,1>;
template struct AdaptiveIntegrationMethod<3,3,3>;

template <size_t dim>
QuadRule<dim-1> choose_sinh_quad(double S, double l, Vec<double,dim-1> singular_pt);

static const auto G2d = gauss(10);
template <>
inline QuadRule<1> choose_sinh_quad<2>(double S, double l, Vec<double,1> singular_pt) {
    assert(l > 0);
    assert(S > 0);
    if ((l / S) > 1) {
        return G2d;
    }
    else {
        size_t n = static_cast<size_t>(10.0 * (1 + std::log(S / l)));
        return sinh_transform(gauss(n), singular_pt[0], l, false);
    }
}

static const auto G3d = tri_gauss(10);
template <>
inline QuadRule<2> choose_sinh_quad<3>(double S, double l, Vec<double,2> singular_pt) {
    assert(l > 0);
    assert(S > 0);
    if ((l / S) > 0.5) {
        return G3d; 
    }
    else {
        size_t n = static_cast<size_t>(8.0 * (1 + std::log(S / l)));
        auto q = sinh_sigmoidal_transform(gauss(2.0 * n), gauss(n), 
            singular_pt[0], singular_pt[1], l, false);
        return q;
    }
}

template <size_t dim, size_t R, size_t C>
Vec<Vec<Vec<double,C>,R>,dim> sinh_integrate(const Kernel<dim,R,C>& k, double tolerance, 
        const IntegralTerm<dim,R,C>& term, const NearestPoint<dim>& near_pt,
        const Vec<double,dim>& nf_obs_pt) 
{
    assert(near_pt.distance > 0);
    auto S = term.src_face.length_scale;
    auto l = near_pt.distance;
    auto q = choose_sinh_quad<dim>(S, l, near_pt.ref_pt);
    auto integrals = zeros<Vec<Vec<Vec<double,C>,R>,dim>>::make();
    for (size_t i = 0; i < q.size(); i++) {
        integrals += term.eval_point_influence(k, q[i].x_hat, nf_obs_pt) * q[i].w;
    }
    return integrals;
}

template <size_t dim, size_t R, size_t C>
Vec<Vec<Vec<double,C>,R>,dim> 
SinhIntegrationMethod<dim,R,C>::compute_singular(const IntegralTerm<dim,R,C>& term,
    const NearestPoint<dim>& nearest_pt) const 
{
    std::vector<Vec<Vec<Vec<double,C>,R>,dim>> steps(qs.n_singular_steps);

    for (int step_idx = 0; step_idx < qs.n_singular_steps; step_idx++) {
        auto step_loc = get_step_loc(term.obs, qs.singular_steps[step_idx]);
        auto shifted_near_pt = FarNearLogic<dim>{qs.far_threshold, 1.0}
            .decide(step_loc, term.src_face);
        steps[step_idx] = sinh_integrate(K, qs.near_tol, term, shifted_near_pt, step_loc);
    }

    return richardson_limit(2, steps);
}

template <size_t dim, size_t R, size_t C>
Vec<Vec<Vec<double,C>,R>,dim> 
SinhIntegrationMethod<dim,R,C>::compute_nearfield(const IntegralTerm<dim,R,C>& term,
    const NearestPoint<dim>& nearest_pt) const 
{
    return sinh_integrate(K, qs.near_tol, term, nearest_pt, term.obs.loc);
}

template <size_t dim, size_t R, size_t C>
Vec<Vec<Vec<double,C>,R>,dim> 
SinhIntegrationMethod<dim,R,C>::compute_farfield(const IntegralTerm<dim,R,C>& term,
    const NearestPoint<dim>& nearest_pt) const 
{
    return adaptive.compute_farfield(term, nearest_pt);
}

template struct SinhIntegrationMethod<2,1,1>;
template struct SinhIntegrationMethod<2,2,2>;
template struct SinhIntegrationMethod<3,1,1>;
template struct SinhIntegrationMethod<3,3,3>;

} //end namespace tbem

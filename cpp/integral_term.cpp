#include <assert.h>
#include "integral_term.h"
#include "richardson.h"
#include "adaptive_quad.h"
#include "sinh_quad.h"
#include "numerics.h"
#include "gte_wrapper.h"

namespace tbem {

template <size_t dim, size_t R, size_t C>
Vec<Vec<Vec<double,C>,R>,dim> IntegralTerm<dim,R,C>::eval_point_influence(
    const Kernel<dim,R,C>& k, const Vec<double,dim-1>& x_hat,
    const Vec<double,dim>& moved_obs_loc) const 
{
    const auto src_pt = ref_to_real(x_hat, src_face.facet);
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
Vec<Vec<Vec<double,C>,R>,dim> 
IntegrationStrategy<dim,R,C>::compute_term(const IntegralTerm<dim,R,C>& term) const
{
    auto nearest_pt = closest_pt_facet(term.obs.loc, term.src_face.facet);
    auto length_scale = term.src_face.length_scale;
    bool nearfield = nearest_pt.distance < far_threshold * length_scale;
    if (nearfield) {
        bool very_close = nearest_pt.distance < singular_threshold * length_scale;
        bool nonzero_richardson_direction =
            hypot2(term.obs.richardson_dir) > 1e-10 * length_scale;
        if (very_close && nonzero_richardson_direction) { 
            return compute_singular(term);
        } else {
            return nearfield_integrator->compute_nearfield(*K, term, nearest_pt);
        }
    }
    return compute_farfield(term);
}

template <size_t dim>
Vec<double,dim> get_step_loc(const ObsPt<dim>& obs, double step_size) {
    const double safe_dist_ratio = 1.0;
    double step_distance = safe_dist_ratio * step_size;
    return obs.loc + step_distance * obs.richardson_dir;
}

template <size_t dim, size_t R, size_t C>
Vec<Vec<Vec<double,C>,R>,dim> 
IntegrationStrategy<dim,R,C>::compute_singular(const IntegralTerm<dim,R,C>& term) const
{
    std::vector<Vec<Vec<Vec<double,C>,R>,dim>> steps(singular_steps.size());

    for (size_t step_idx = 0; step_idx < singular_steps.size(); step_idx++) {
        auto step_loc = get_step_loc(term.obs, singular_steps[step_idx]);
        ObsPt<dim> shifted_obs_pt{step_loc, term.obs.normal, term.obs.richardson_dir};
        auto shifted_nearest_pt = closest_pt_facet(step_loc, term.src_face.facet);
        IntegralTerm<dim,R,C> shifted_term{shifted_obs_pt, term.src_face};
        steps[step_idx] = nearfield_integrator->compute_nearfield(
            *K, shifted_term, shifted_nearest_pt
        );
    }

    return richardson_limit(2.0, steps);
}

template <size_t dim, size_t R, size_t C>
Vec<Vec<Vec<double,C>,R>,dim> 
IntegrationStrategy<dim,R,C>::compute_farfield(const IntegralTerm<dim,R,C>& term) const
{
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


template <size_t dim, size_t R, size_t C>
Vec<Vec<Vec<double,C>,R>,dim> 
AdaptiveIntegrator<dim,R,C>::compute_nearfield(const Kernel<dim,R,C>& K, 
    const IntegralTerm<dim,R,C>& term,
    const NearestPoint<dim>& nearest_pt) 
{
    (void)nearest_pt;
    return UnitFacetAdaptiveIntegrator<dim>()(K, near_tol, term, term.obs.loc);
}

template struct AdaptiveIntegrator<2,1,1>;
template struct AdaptiveIntegrator<2,2,2>;
template struct AdaptiveIntegrator<3,1,1>;
template struct AdaptiveIntegrator<3,3,3>;

template <size_t dim>
QuadRule<dim-1> make_sinh_quad(size_t order, Vec<double,dim-1> singular_pt,
    double scaled_distance);

template <>
QuadRule<1> make_sinh_quad<2>(size_t order, Vec<double,1> singular_pt,
    double scaled_distance)
{
    return sinh_transform(gauss(order), singular_pt[0], scaled_distance, false);
}

template <>
QuadRule<2> make_sinh_quad<3>(size_t order, Vec<double,2> singular_pt,
    double scaled_distance)
{
    return sinh_sigmoidal_transform(
        gauss(2 * order), gauss(order), singular_pt[0], singular_pt[1],
        scaled_distance, false
    );
}

template <size_t dim>
QuadRule<dim-1> choose_sinh_quad(size_t far_order, size_t order_growth_rate,
    double S, double l, Vec<double,dim-1> singular_pt)
{
    double scaled_distance = l / S;
    bool far_enough_that_sinh_quadrature_is_unnecessary = scaled_distance > 0.5;
    if (far_enough_that_sinh_quadrature_is_unnecessary) {
        return gauss_facet<dim>(10); 
    }
    size_t n = far_order + (order_growth_rate * (-std::log(scaled_distance)));
    return make_sinh_quad<dim>(n, singular_pt, scaled_distance);
}

template <size_t dim, size_t R, size_t C>
Vec<Vec<Vec<double,C>,R>,dim> 
SinhIntegrator<dim,R,C>::compute_nearfield(const Kernel<dim,R,C>& K, 
    const IntegralTerm<dim,R,C>& term, const NearestPoint<dim>& nearest_pt)
{
    auto l = nearest_pt.distance;
    assert(l > 0);
    auto S = term.src_face.length_scale;
    auto q = choose_sinh_quad<dim>(sinh_order, sinh_order, S, l, nearest_pt.ref_pt);
    auto integrals = zeros<Vec<Vec<Vec<double,C>,R>,dim>>::make();
    for (size_t i = 0; i < q.size(); i++) {
        integrals += term.eval_point_influence(K, q[i].x_hat, term.obs.loc) * q[i].w;
    }
    return integrals;
}

template struct SinhIntegrator<2,1,1>;
template struct SinhIntegrator<2,2,2>;
template struct SinhIntegrator<3,1,1>;
template struct SinhIntegrator<3,3,3>;

} //end namespace tbem

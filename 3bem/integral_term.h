#ifndef __PPPPPPPPKPKPKPKP_INTEGRAL_TERM_H
#define __PPPPPPPPKPKPKPKP_INTEGRAL_TERM_H

#include <cassert>
#include <cstdlib>
#include <iomanip> 
#include "vec.h"
#include "kernel.h"
#include "quadrature.h"
#include "adaptive_quad.h"
#include "numerics.h"
#include "obs_pt.h"
#include "facet_info.h"
#include "closest_pt.h"

namespace tbem {

/* Data transfer object for computing integral terms. 
 * Values are stored by reference for efficiency's sake.
 * This means that the responsibility for maintaining their lifetime is on the
 * user.
 */
template <size_t dim, size_t R, size_t C>
struct IntegralTerm {
    const QuadStrategy<dim>& qs;
    const Kernel<dim,R,C>& k;
    const ObsPt<dim>& obs;
    const FacetInfo<dim>& src_face;
};

template <size_t dim, size_t R, size_t C>
IntegralTerm<dim,R,C> make_integral_term(const QuadStrategy<dim>& qs,
        const Kernel<dim,R,C>& k, const ObsPt<dim>& obs, const FacetInfo<dim>& src_face)
{
    return IntegralTerm<dim,R,C>{qs, k, obs, src_face};
}


/* Given observation point information and source face information, the
 * fucntion evaluates the influence of a single source quadrature point
 * on the observation point.
 */
template <size_t dim, size_t R, size_t C>
Vec<Vec<Vec<double,C>,R>,dim> eval_point_influence(const Vec<double,dim-1>& x_hat,
                                  const Kernel<dim,R,C>& kernel,
                                  const FacetInfo<dim>& face,
                                  const Vec<double,dim>& obs_loc,
                                  const Vec<double,dim>& obs_n) {
    const auto src_pt = ref_to_real(x_hat, face.face);
    const auto d = src_pt - obs_loc;
    const auto r2 = dot_product(d, d);
    const auto kernel_val = kernel(r2, d, face.normal, obs_n);
    return outer_product(linear_basis(x_hat), kernel_val * face.jacobian);
}

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

template <size_t dim, size_t R, size_t C>
struct IntegrationMethodI {
    typedef Vec<Vec<Vec<double,C>,R>,dim> OutputType;

    virtual OutputType compute(const IntegralTerm<dim,R,C>& term,
        const NearestPoint<dim>& near_pt);
};

template <size_t dim, size_t R, size_t C>
struct FixedIntegrate: public IntegrationMethodI<dim,R,C> {
    typedef Vec<Vec<Vec<double,C>,R>,dim> OutputType;

    const QuadRule<dim-1> q;

    FixedIntegrate(const QuadRule<dim-1>& q):
        q(q)
    {}

    virtual OutputType compute(const IntegralTerm<dim,R,C>& term,
        const NearestPoint<dim>& near_pt) 
    {
        auto integrals = zeros<OutputType>::make();
        for (size_t i = 0; i < term.qs.src_far_quad.size(); i++) {
            integrals += eval_point_influence<dim>(
                            term.qs.src_far_quad[i].x_hat,
                            term.k, term.src_face,
                            term.obs.loc, term.obs.normal
                        ) * term.qs.src_far_quad[i].w;
        }
        return integrals;
    }
};

template <size_t dim, size_t R, size_t C>
struct AdaptiveIntegrate: public IntegrationMethodI<dim,R,C> {

};

// struct AdaptiveIntegrationMethod: public IntegrationMethodI {
//     const QuadRule<dim-1> obs_quad;
//     const QuadRule<dim-1> src_far_quad;
//     
//     const double singular_threshold;
//     const double far_threshold;
//     const int n_singular_steps;
//     const std::vector<double> singular_steps;
//     const double near_tol;
//     QuadStrategy(int obs_order, int src_far_order, int n_singular_steps,
//                  double far_threshold, double near_tol);
// };
// 

template <size_t dim>
struct UnitFacetAdaptiveIntegrator {
    template <size_t R, size_t C>
    Vec<Vec<Vec<double,C>,R>,dim> 
    operator()(const IntegralTerm<dim,R,C>& term, const NearestPoint<dim>& near_pt,
               const Vec<double,dim>& nf_obs_pt);
};


static const auto G = gauss(10);
inline QuadRule<1> choose_2d_quad(double S, double l, double x0) {
    assert(l > 0);
    assert(S > 0);
    static size_t max_n = 0;
    if ((l / S) > 1) {
        return G;
    }
    else {
        size_t n = static_cast<size_t>(10.0 * (1 + std::log(S / l)));
        max_n = std::max(max_n, n);
        // std::cout << max_n << std::endl;
        return sinh_transform(n, x0, l);
    }
}

template <>
struct UnitFacetAdaptiveIntegrator<2> {
    template <size_t R, size_t C>
    Vec<Vec<Vec<double,C>,R>,2> operator()(const IntegralTerm<2,R,C>& term, 
                const NearestPoint<2>& near_pt, const Vec<double,2>& nf_obs_pt) {
        assert(near_pt.distance > 0);
        auto S = term.src_face.length_scale;
        auto l = near_pt.distance;
        auto q = choose_2d_quad(S, l, near_pt.ref_pt[0]);
        auto integrals = zeros<Vec<Vec<Vec<double,C>,R>,2>>::make();
        for (size_t i = 0; i < q.size(); i++) {
            integrals += eval_point_influence<2>(
                            q[i].x_hat,
                            term.k, term.src_face,
                            nf_obs_pt, term.obs.normal
                        ) * q[i].w;
        }
        return integrals;
    }
};

template <>
struct UnitFacetAdaptiveIntegrator<3> {
    template <size_t R, size_t C>
    Vec<Vec<Vec<double,C>,R>,3> operator()(const IntegralTerm<3,R,C>& term, 
            const NearestPoint<3>& near_pt, const Vec<double,3>& nf_obs_pt) {
        // auto q = sinh_sigmoidal_transform(60, 30, near_pt.ref_pt[0],
        //     near_pt.ref_pt[1], near_pt.distance);
        // auto integrals = zeros<Vec<Vec<Vec<double,C>,R>,3>>::make();
        // for (size_t i = 0; i < q.size(); i++) {
        //     integrals += eval_point_influence<3>(
        //                     q[i].x_hat,
        //                     term.k, term.src_face,
        //                     nf_obs_pt, term.obs.normal
        //                 ) * q[i].w;
        // }
        auto correct = adaptive_integrate<Vec<Vec<Vec<double,C>,R>,3>>(
            [&] (double x_hat) {
                if (x_hat == 1.0) {
                    return zeros<Vec<Vec<Vec<double,C>,R>,3>>::make();
                }
                return adaptive_integrate<Vec<Vec<Vec<double,C>,R>,3>>(
                    [&] (double y_hat) {
                        Vec<double,2> q_pt = {x_hat, y_hat};
                        return eval_point_influence<3>(q_pt, term.k, term.src_face,
                                                 nf_obs_pt, term.obs.normal);
                    }, 0.0, 1 - x_hat, term.qs.near_tol);
            }, 0.0, 1.0, term.qs.near_tol);
        // auto error = fabs(correct - integrals) / correct;
        // // std::cout << error << std::endl;
        // assert(all(error < 1e-4 * ones<Vec<Vec<Vec<double,C>,R>,3>>::make()));
        return correct;
    }
};

template <size_t dim, size_t R, size_t C> 
Vec<Vec<Vec<double,C>,R>,dim> 
compute_nearfield(const IntegralTerm<dim,R,C>& term, const NearestPoint<dim>& near_pt,
                   const Vec<double,dim>& nf_obs_pt) {
    UnitFacetAdaptiveIntegrator<dim> integrator;
    return integrator(term, near_pt, nf_obs_pt);
}

template <size_t dim, size_t R, size_t C>
Vec<double,dim> get_step_loc(const IntegralTerm<dim,R,C>& term, int step_idx) {
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

template <size_t dim, size_t R, size_t C> 
Vec<Vec<Vec<double,C>,R>,dim> compute_as_limit(const IntegralTerm<dim,R,C>& term,
    const NearestPoint<dim>& near_pt) 
{
    std::vector<Vec<Vec<Vec<double,C>,R>,dim>> 
        near_steps(term.qs.n_singular_steps);

    for (int step_idx = 0; step_idx < term.qs.n_singular_steps; step_idx++) {
        auto step_loc = get_step_loc(term, step_idx);
        auto shifted_near_pt = FarNearLogic<dim>{term.qs.far_threshold, 3.0}
            .decide(step_loc, term.src_face);
        near_steps[step_idx] = compute_nearfield<dim>(term, shifted_near_pt, step_loc);
    }

    return richardson_limit(near_steps);
}

template <size_t dim, size_t R, size_t C>
Vec<Vec<Vec<double,C>,R>,dim> compute_far_term(const IntegralTerm<dim,R,C>& term) {
    auto integrals = zeros<Vec<Vec<Vec<double,C>,R>,dim>>::make();
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
                       const Vec<Vec<double,dim>,dim>& vs) {
    double res = dist2(pt, vs[0]);
    for (int d = 1; d < dim; d++) {
        res = std::min(res, dist2(pt, vs[d])); 
    }
    return res;
}
/* Compute the full influence of a source facet on an observation point, given
 * a kernel function/Green's function
 */
template <size_t dim, size_t R, size_t C>
Vec<Vec<Vec<double,C>,R>,dim> compute_term(const IntegralTerm<dim,R,C>& term) {
    FarNearLogic<dim> far_near_logic{term.qs.far_threshold, 3.0};
    auto nearest_pt = far_near_logic.decide(term.obs.loc, term.src_face);
    switch (nearest_pt.type) {
        case FarNearType::Singular:
            return compute_as_limit(term, nearest_pt);
            break;
        case FarNearType::Nearfield:
            return compute_nearfield<dim>(term, nearest_pt, term.obs.loc);
            break;
        case FarNearType::Farfield:
            return compute_far_term(term);
            break;
    }
    throw std::exception();
}



} //end namespace tbem

#endif

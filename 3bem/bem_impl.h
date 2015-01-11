#ifndef __QWEQWEJKLJSHDFKJNV_BEM_IMPL_H
#define __QWEQWEJKLJSHDFKJNV_BEM_IMPL_H

/* Included at the bottom of bem.h to separate implementation from
 * definition despite the requirement that templates be fully defined
 * in a header file
 */

template <size_t dim>
ObsPt<dim> ObsPt<dim>::from_face(const Vec<double,dim-1>& ref_loc,
                                 const FacetInfo<dim>& obs_face) {
    const int basis_order = 1;
    return {
        obs_face.length_scale / basis_order,
        ref_to_real(ref_loc, obs_face.face.vertices),
        obs_face.normal,
        obs_face.normal 
    };
}

template <>
FacetInfo<3> FacetInfo<3>::build(const Facet<3>& facet){
    auto unscaled_n = unscaled_normal(facet.vertices);
    auto area = tri_area(unscaled_n);
    auto length_scale = std::sqrt(area);
    auto jacobian = area * inv_ref_facet_area;
    auto normal = unscaled_n / jacobian;
    return FacetInfo<3>{facet, area, length_scale, jacobian, normal};
}

template <>
FacetInfo<2> FacetInfo<2>::build(const Facet<2>& facet){
    auto unscaled_n = unscaled_normal(facet.vertices);
    auto area_scale = hypot2(unscaled_n);
    auto length = std::sqrt(area_scale);
    auto jacobian = length * inv_ref_facet_area;
    auto normal = unscaled_n / length;
    return FacetInfo<2>{facet, area_scale, length, jacobian, normal};
}

template <size_t dim, typename KT>
Vec<typename KT::OperatorType,dim> eval_point_influence(const Vec<double,dim-1>& x_hat,
                                  const KT& kernel,
                                  const FacetInfo<dim>& face,
                                  const Vec<double,dim>& obs_loc,
                                  const Vec<double,dim>& obs_n) {
    const auto src_pt = ref_to_real(x_hat, face.face.vertices);
    const auto d = src_pt - obs_loc;
    const auto r2 = dot_product(d, d);
    const auto kernel_val = kernel(r2, d, face.normal, obs_n);
    return outer_product(linear_basis(x_hat), kernel_val * face.jacobian);
}

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
Vec<typename KT::OperatorType,dim> compute_near_term(const IntegralTerm<dim, KT>& term) {
    const double singular_threshold = 3.0;
    if (term.appx_pt_face_dist_squared < 
            singular_threshold * term.src_face.area_scale) { 
        return compute_as_limit(term);
    } else {
        return compute_adaptively<dim>(term, term.obs.loc);
    }
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

template <size_t dim, typename KT>
Vec<typename KT::OperatorType,dim> compute_term(const IntegralTerm<dim,KT>& term) {
    if (term.appx_pt_face_dist_squared < 
            pow(term.qs.far_threshold, 2) * term.src_face.area_scale) {
        return compute_near_term(term);
    } else {
        return compute_far_term(term);
    }
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

template <size_t dim, typename KT>
std::vector<typename KT::OperatorType> 
integral_equation_vector(const Problem<dim,KT>& p, const QuadStrategy<dim>& qs,
                         const ObsPt<dim>& obs) {
    int n_out_dofs = dim * p.src_mesh.facets.size();
    std::vector<typename KT::OperatorType> result(n_out_dofs);
    for (size_t i = 0; i < p.src_mesh.facets.size(); i++) {
        auto src_face = FacetInfo<dim>::build(p.src_mesh.facets[i]);
        const double dist2 = appx_face_dist2<dim>(obs.loc, src_face.face.vertices);
        auto term = make_integral_term(qs, p.K, obs, src_face, dist2);
        auto integrals = compute_term<dim>(term);
        for (int b = 0; b < dim; b++) {
            result[dim * i + b] = integrals[b];
        }
    }
    return result;
}

template <size_t dim, typename KT>
double eval_integral_equation(const Problem<dim,KT>& p, const QuadStrategy<dim>& qs,
                              const ObsPt<dim>& obs) {
    double result = 0.0;
    auto row = integral_equation_vector(p, qs, obs);
    for (size_t i = 0; i < row.size(); i++) {
        result += row[i] * p.src_strength[i];
    }
    return result;
}

template <size_t dim, typename KT>
std::vector<typename KT::OperatorType> interact_matrix(const Problem<dim,KT>& p,
                                    const QuadStrategy<dim>& qs) 
{
    size_t n_obs_dofs = dim * p.obs_mesh.facets.size();
    size_t n_src_dofs = dim * p.src_mesh.facets.size();
    std::vector<typename KT::OperatorType> matrix(n_obs_dofs * n_src_dofs, 
            zeros<typename KT::OperatorType>::make());
#pragma omp parallel for
    for (size_t obs_idx = 0; obs_idx < p.obs_mesh.facets.size(); obs_idx++) {
        auto obs_face = FacetInfo<dim>::build(p.obs_mesh.facets[obs_idx]);
        for (size_t obs_q = 0; obs_q < qs.obs_quad.size(); obs_q++) {
            auto pt = ObsPt<dim>::from_face(qs.obs_quad[obs_q].x_hat, obs_face);

            const auto row = integral_equation_vector(p, qs, pt);

            const auto basis = linear_basis(qs.obs_quad[obs_q].x_hat);

            for (int v = 0; v < dim; v++) {
                int b = dim * obs_idx + v;
                for (size_t i = 0; i < n_src_dofs; i++) {
                    matrix[b * n_src_dofs + i] +=
                        basis[v] *
                        row[i] *
                        qs.obs_quad[obs_q].w *
                        obs_face.jacobian;
                }
            }
        }
    }
    return matrix;
}

template <size_t dim, typename KT>
std::vector<typename KT::OutType> mass_term(const Problem<dim,KT>& p,
    const QuadStrategy<dim>& qs) 
{
    int n_obs_dofs = dim * p.obs_mesh.facets.size();
    std::vector<typename KT::OutType> integrals(
        n_obs_dofs, 
        zeros<typename KT::OutType>::make()
    );

    for (size_t obs_idx = 0; obs_idx < p.obs_mesh.facets.size(); obs_idx++) {
        auto obs_face = FacetInfo<dim>::build(p.obs_mesh.facets[obs_idx]);
        for (size_t obs_q = 0; obs_q < qs.obs_quad.size(); obs_q++) {
            auto qpt = qs.obs_quad[obs_q];

            int dof = dim * obs_idx;
            Vec<typename KT::InType,dim> face_vals;
            for (int d = 0; d < dim; d++) {
                face_vals[d] = p.src_strength[dof + d];
            }

            auto interp_val = dot_product(linear_basis(qpt.x_hat), face_vals);
            auto kernel_val = p.K.call_with_no_params();
            auto out_val = dot_product(interp_val, kernel_val);

            auto basis = linear_basis(qpt.x_hat);
            for(int v = 0; v < dim; v++) {
                integrals[dof + v] += obs_face.jacobian * basis[v] * 
                                      out_val * qpt.w;
            }
        }
    }
    return integrals;
}

template <typename KT>
std::vector<typename KT::OutType>
bem_mat_mult(const std::vector<typename KT::OperatorType>& A, const KT& k,
             int n_obs_dofs, const std::vector<typename KT::InType>& x) {

    assert(n_obs_dofs * x.size() == A.size());
    std::vector<typename KT::OutType> res(n_obs_dofs, 
                                          zeros<typename KT::OutType>::make());
#pragma omp parallel for
    for (int i = 0; i < n_obs_dofs; i++) {
        for (size_t j = 0; j < x.size(); j++) {
            res[i] += dot_product(x[j], A[i * x.size() + j]);
        }
    }
    return res;
}


template <size_t dim, typename KT>
std::vector<typename KT::OutType> direct_interact(const Problem<dim,KT>& p,
                                                  const QuadStrategy<dim>& qs) {
    auto matrix = interact_matrix(p, qs);
    assert(p.obs_mesh.facets.size() * dim * p.src_strength.size() == matrix.size());
    return bem_mat_mult(matrix, p.K, p.obs_mesh.facets.size() * dim, p.src_strength);
}

template <>
double get_len_scale<3>(Mesh<3>& mesh, int which_face, int q) {
    return std::sqrt(tri_area(mesh.facets[which_face].vertices)) / q;
}

template <>
double get_len_scale<2>(Mesh<2>& mesh, int which_face, int q) {
    return dist(mesh.facets[which_face].vertices[1],
                mesh.facets[which_face].vertices[0]) / q;
}

#endif

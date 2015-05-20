#ifndef ASLJLJQQJJQJQJQJQJQJQJQJQ_NEARFIELD_OPERATOR_H
#define ASLJLJQQJJQJQJQJQJQJQJQJQ_NEARFIELD_OPERATOR_H

#include "sparse_operator.h"
#include "integral_term.h"
#include "obs_pt.h"

namespace tbem {

template <size_t dim, size_t R, size_t C>
using GalerkinNearfieldFnc = std::function<Vec<Vec<Vec<double,C>,R>,dim>(
    const IntegralTerm<dim,R,C>&, const NearestPoint<dim>&
    )>;

template <size_t dim>
std::vector<ObsPt<dim>> nearfield_obs_pts(const Mesh<dim>& obs_mesh,
    const QuadRule<dim-1>& obs_quad, const Mesh<dim>& all_mesh)
{
    std::vector<ObsPt<dim>> out;
    for (size_t obs_idx = 0; obs_idx < obs_mesh.facets.size(); obs_idx++) {
        auto obs_face = FacetInfo<dim>::build(obs_mesh.facets[obs_idx]);
        for (size_t obs_q = 0; obs_q < obs_quad.size(); obs_q++) {
            auto pt = ObsPt<dim>::away_from_nearest_facets(
                obs_quad[obs_q].x_hat, obs_face, all_mesh
            );
            out.push_back(pt);
        }
    }
    return out;
}

template <size_t dim, size_t R, size_t C>
SparseOperator nearfield_helper(const Mesh<dim>& obs_mesh,
    const Mesh<dim>& src_mesh, const IntegrationMethodI<dim,R,C>& mthd,
    const GalerkinNearfieldFnc<dim,R,C>& f, const Mesh<dim>& all_mesh)
{
    const auto& obs_quad = mthd.get_obs_quad();

    size_t n_obs_pts = obs_mesh.n_facets() * obs_quad.size();
    size_t n_src_dofs = src_mesh.n_dofs();

    auto src_facet_info = get_facet_info(src_mesh);
    std::vector<MatrixEntry> entries;
#pragma omp parallel for
    for (size_t obs_idx = 0; obs_idx < obs_mesh.facets.size(); obs_idx++) {
        auto obs_face = FacetInfo<dim>::build(obs_mesh.facets[obs_idx]);
        for (size_t obs_q = 0; obs_q < obs_quad.size(); obs_q++) {
            auto pt = ObsPt<dim>::away_from_nearest_facets(
                obs_quad[obs_q].x_hat, obs_face, all_mesh
            );
            auto pt_idx = obs_idx * obs_quad.size() + obs_q;

            for (size_t i = 0; i < src_mesh.facets.size(); i++) {
                FarNearLogic<dim> far_near_logic{mthd.far_threshold(), 1.0};
                auto nearest_pt = far_near_logic.decide(pt.loc, src_facet_info[i]);
                if (nearest_pt.type == FarNearType::Farfield) {
                    continue; 
                }
                IntegralTerm<dim,R,C> term{pt, src_facet_info[i]};
                auto eval = f(term, nearest_pt);
                for (size_t basis_idx = 0; basis_idx < dim; basis_idx++) {
                    auto src_dof_idx = i * dim + basis_idx; 
                    for (size_t d1 = 0; d1 < R; d1++) {
                        auto row_idx = d1 * n_obs_pts + pt_idx;
                        for (size_t d2 = 0; d2 < C; d2++) {
                            auto col_idx = d2 * n_src_dofs + src_dof_idx;
#pragma omp critical
                            entries.push_back({
                                row_idx, col_idx, eval[basis_idx][d1][d2]
                            });
                        }
                    }
                }
            }
        }
    }

    return SparseOperator::csr_from_coo(
        R * obs_mesh.n_facets() * obs_quad.size(),
        C * n_src_dofs,
        entries
    );
}

template <size_t dim, size_t R, size_t C>
SparseOperator make_nearfield_operator(const Mesh<dim>& obs_mesh,
    const Mesh<dim>& src_mesh, const IntegrationMethodI<dim,R,C>& mthd,
    const Mesh<dim>& all_mesh) 
{
    GalerkinNearfieldFnc<dim,R,C> f = 
        [&] (const IntegralTerm<dim,R,C>& term, const NearestPoint<dim>& pt) {
            return mthd.compute_term(term, pt); 
        };
    return nearfield_helper(obs_mesh, src_mesh, mthd, f, all_mesh);
}

template <size_t dim, size_t R, size_t C>
SparseOperator make_farfield_correction_operator(const Mesh<dim>& obs_mesh,
    const Mesh<dim>& src_mesh, const IntegrationMethodI<dim,R,C>& mthd,
    const Mesh<dim>& all_mesh) 
{
    GalerkinNearfieldFnc<dim,R,C> f = 
        [&] (const IntegralTerm<dim,R,C>& term, const NearestPoint<dim>& pt) {
            return -mthd.compute_farfield(term, pt); 
        };
    return nearfield_helper(obs_mesh, src_mesh, mthd, f, all_mesh);
}


} //end namespace tbem

#endif

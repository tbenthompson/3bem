#ifndef ASLJLJQQJJQJQJQJQJQJQJQJQ_NEARFIELD_OPERATOR_H
#define ASLJLJQQJJQJQJQJQJQJQJQJQ_NEARFIELD_OPERATOR_H

#include "sparse_operator.h"
#include "integral_term.h"
#include "nearest_neighbors.h"
#include "limit_direction.h"
#include "util.h"

namespace tbem {

template <size_t dim>
std::vector<ObsPt<dim>> galerkin_obs_pts(const Mesh<dim>& obs_mesh,
    const QuadRule<dim-1>& obs_quad, const Mesh<dim>& all_mesh)
{
    std::vector<ObsPt<dim>> out;
    NearfieldFacetFinder<dim> nearfield_finder(all_mesh.facets, 1.0);
    for (size_t obs_idx = 0; obs_idx < obs_mesh.facets.size(); obs_idx++) {
        auto obs_face = FacetInfo<dim>::build(obs_mesh.facets[obs_idx]);
        for (size_t obs_q = 0; obs_q < obs_quad.size(); obs_q++) {
            auto ref_loc = obs_quad[obs_q].x_hat;
            auto loc = ref_to_real(ref_loc, obs_face.face);
            auto nf = nearfield_finder.find(loc);

            auto rich_dir = decide_limit_dir(loc, nf, all_mesh.facets, 0.4, 1e-2);

            out.push_back({loc, obs_face.normal, rich_dir});
        }
    }
    return out;
}

template <size_t dim>
std::vector<ObsPt<dim>> interior_obs_pts(const std::vector<Vec<double,dim>>& locs,
    const std::vector<Vec<double,dim>>& normals, const Mesh<dim>& all_mesh)
{
    std::vector<ObsPt<dim>> out;
    NearfieldFacetFinder<dim> nearfield_finder(all_mesh.facets, 1.0);
    for (size_t pt_idx = 0; pt_idx < locs.size(); pt_idx++) {
        auto nf = nearfield_finder.find(locs[pt_idx]);
        auto rich_dir = decide_limit_dir(locs[pt_idx], nf, all_mesh.facets, 0.4, 1e-2);
        out.push_back({locs[pt_idx], normals[pt_idx], rich_dir});
    }
    return out;
}

template <size_t dim, size_t R, size_t C>
using NearfieldFnc = std::function<Vec<Vec<Vec<double,C>,R>,dim>(
    const IntegralTerm<dim,R,C>&, const NearestPoint<dim>&
    )>;

template <size_t dim, size_t R, size_t C>
SparseOperator nearfield_inner_integral(const std::vector<ObsPt<dim>>& obs_pts,
    const Mesh<dim>& src_mesh, const IntegrationStrategy<dim,R,C>& mthd,
    const NearfieldFnc<dim,R,C>& f)
{
    size_t n_src_dofs = src_mesh.n_dofs();
    auto src_facet_info = get_facet_info(src_mesh);

    std::vector<MatrixEntry> entries;
    NearfieldFacetFinder<dim> nearfield_finder(src_mesh.facets, mthd.far_threshold);
#pragma omp parallel for
    for (size_t pt_idx = 0; pt_idx < obs_pts.size(); pt_idx++) {
        auto pt = obs_pts[pt_idx];
        auto nearfield = nearfield_finder.find(pt.loc);
        for (auto i: nearfield.facet_indices) {
            //TODO: This is accessing an internal of integrationstrategy, maybe 
            //that's a bit of a code smell, seems like this should be within the
            //integration strategy
            //TODO: This should be merged in with the nearfieldfacetfinder so that
            //nearfield facets are found once and for all and in a O(logn) time
            //fashion
            //At the moment this is done twice. Working out the best way to organize
            //the data structure will be somewhat complex
            FarNearLogic<dim> far_near_logic{mthd.far_threshold, 1.0};
            auto nearest_pt = far_near_logic.decide(pt.loc, src_facet_info[i]);
            if (nearest_pt.type == FarNearType::Farfield) {
                continue; 
            }

            auto eval = f({pt, src_facet_info[i]}, nearest_pt);
            for (size_t basis_idx = 0; basis_idx < dim; basis_idx++) {
                auto src_dof_idx = i * dim + basis_idx; 
                for (size_t d1 = 0; d1 < R; d1++) {
                    auto row_idx = d1 * obs_pts.size() + pt_idx;
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

    return SparseOperator::csr_from_coo(
        R * obs_pts.size(),
        C * n_src_dofs,
        entries
    );
}

template <size_t dim, size_t R, size_t C>
SparseOperator make_nearfield_operator(
    const std::vector<ObsPt<dim>>& obs_pts, const Mesh<dim>& src_mesh,
    const IntegrationStrategy<dim,R,C>& mthd) 
{
    NearfieldFnc<dim,R,C> f = 
        [&] (const IntegralTerm<dim,R,C>& term, const NearestPoint<dim>& pt) {
            return mthd.compute_term(term, pt); 
        };
    return nearfield_inner_integral(obs_pts, src_mesh, mthd, f);
}

template <size_t dim, size_t R, size_t C>
SparseOperator make_farfield_correction_operator(
    const std::vector<ObsPt<dim>>& obs_pts, const Mesh<dim>& src_mesh,
    const IntegrationStrategy<dim,R,C>& mthd) 
{
    NearfieldFnc<dim,R,C> f = 
        [&] (const IntegralTerm<dim,R,C>& term, const NearestPoint<dim>& pt) {
            return -mthd.compute_farfield(term, pt); 
        };
    return nearfield_inner_integral(obs_pts, src_mesh, mthd, f);
}


} //end namespace tbem

#endif

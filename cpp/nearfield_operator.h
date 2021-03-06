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
    std::vector<Vec<double,dim>> locs;
    std::vector<Vec<double,dim>> normals;
    for (size_t obs_idx = 0; obs_idx < obs_mesh.facets.size(); obs_idx++) {
        auto obs_face = FacetInfo<dim>::build(obs_mesh.facets[obs_idx]);
        for (size_t obs_q = 0; obs_q < obs_quad.size(); obs_q++) {
            locs.push_back(ref_to_real(obs_quad[obs_q].x_hat, obs_face.facet));
            normals.push_back(obs_face.normal);
        }
    }
    return interior_obs_pts(locs, normals, all_mesh);
}

template <size_t dim>
std::vector<ObsPt<dim>> interior_obs_pts(const std::vector<Vec<double,dim>>& locs,
    const std::vector<Vec<double,dim>>& normals, const Mesh<dim>& all_mesh)
{
    std::vector<ObsPt<dim>> out;
    NearfieldFacetFinder<dim> nearfield_finder(all_mesh.facets, 1.0);
    for (size_t pt_idx = 0; pt_idx < locs.size(); pt_idx++) {
        auto nf = nearfield_finder.find(locs[pt_idx]);
        auto rich_dir = decide_limit_dir(
            locs[pt_idx], nf, all_mesh.facets, 0.4, 1e-2
        );
        out.push_back({locs[pt_idx], normals[pt_idx], rich_dir});
    }
    return out;
}

template <size_t dim, size_t R, size_t C>
using NearfieldFnc =
    std::function<Vec<Vec<Vec<double,C>,R>,dim>(const IntegralTerm<dim,R,C>&)>;

template <size_t dim, size_t R, size_t C>
SparseOperator nearfield_inner_integral(const std::vector<ObsPt<dim>>& obs_pts,
    const NearfieldFacetFinder<dim>& nearfield_finder,
    const NearfieldFnc<dim,R,C>& integrate)
{
    //TODO: an idea for logging a bit of stuff
    // logger.log_nearfield_inner_integral(obs_pts, src_mesh, mthd)
    // logger.log_method_used(mthd)
    size_t n_src_dofs = nearfield_finder.n_underlying_dofs();

    std::vector<MatrixEntry> entries;
#pragma omp parallel for
    for (size_t pt_idx = 0; pt_idx < obs_pts.size(); pt_idx++) {
        auto pt = obs_pts[pt_idx];
        auto nearfield = nearfield_finder.find(pt.loc);
        for (size_t i = 0; i < nearfield.facet_indices.size(); i++) {

            auto facet_idx = nearfield.facet_indices[i];
            auto facet_info = nearfield_finder.get_facet_info(facet_idx);
            auto matrix_entries = integrate({pt, facet_info});

            for (size_t basis_idx = 0; basis_idx < dim; basis_idx++) {
                auto src_dof_idx = facet_idx * dim + basis_idx; 
                for (size_t d1 = 0; d1 < R; d1++) {
                    auto row_idx = d1 * obs_pts.size() + pt_idx;
                    for (size_t d2 = 0; d2 < C; d2++) {
                        auto col_idx = d2 * n_src_dofs + src_dof_idx;
#pragma omp critical
                        entries.push_back({
                            row_idx, col_idx, matrix_entries[basis_idx][d1][d2]
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
        [&] (const IntegralTerm<dim,R,C>& term) {
            return mthd.compute_term(term); 
        };
    NearfieldFacetFinder<dim> nearfield_finder(src_mesh.facets, mthd.far_threshold);
    return nearfield_inner_integral(obs_pts, nearfield_finder, f);
}

template <size_t dim, size_t R, size_t C>
SparseOperator make_farfield_correction_operator(
    const std::vector<ObsPt<dim>>& obs_pts, const Mesh<dim>& src_mesh,
    const IntegrationStrategy<dim,R,C>& mthd) 
{
    NearfieldFnc<dim,R,C> f = 
        [&] (const IntegralTerm<dim,R,C>& term) {
            return -mthd.compute_farfield(term); 
        };
    NearfieldFacetFinder<dim> nearfield_finder(src_mesh.facets, mthd.far_threshold);
    return nearfield_inner_integral(obs_pts, nearfield_finder, f);
}


} //end namespace tbem

#endif

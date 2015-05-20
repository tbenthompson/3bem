#ifndef __AAAAAAAAA_BEM_H
#define __AAAAAAAAA_BEM_H

#include <cassert>
#include "nearfield_operator.h"
#include "dense_operator.h"
#include "obs_pt.h"
#include "util.h"
#include "integral_term.h"

namespace tbem {

template <size_t dim>
std::vector<ObsPt<dim>> interior_obs_pts(const std::vector<Vec<double,dim>>& locs,
    const std::vector<Vec<double,dim>>& normals, const Mesh<dim>& all_mesh)
{
    std::vector<ObsPt<dim>> out;
    for (size_t pt_idx = 0; pt_idx < locs.size(); pt_idx++) {
        auto nf = nearest_facets(locs[pt_idx], all_mesh.facets);
        auto rich_dir = decide_richardson_dir(locs[pt_idx], nf);
        auto rich_length = hypot(rich_dir);
        ObsPt<dim> pt{
            rich_length / 5.0, locs[pt_idx], normals[pt_idx], rich_dir / rich_length
        };
        out.push_back(pt);
    }
    return out;
}

template <size_t dim, size_t R, size_t C>
DenseOperator mesh_to_points_operator(const std::vector<Vec<double,dim>>& locs,
    const std::vector<Vec<double,dim>>& normals, const Mesh<dim>& src_mesh,
    const IntegrationMethodI<dim,R,C>& mthd, const Mesh<dim>& all_mesh)
{
    auto n_obs_dofs = locs.size();
    auto n_src_dofs = src_mesh.n_dofs();

    DenseOperator op(R * n_obs_dofs, C * n_src_dofs, 0.0);

    auto src_facet_info = get_facet_info(src_mesh);
    auto obs_pts = interior_obs_pts(locs, normals, all_mesh);
#pragma omp parallel for
    for (size_t pt_idx = 0; pt_idx < obs_pts.size(); pt_idx++) {
        auto pt = obs_pts[pt_idx];
        std::vector<Vec<Vec<double,C>,R>> result(src_mesh.n_dofs());
        FarNearLogic<dim> far_near_logic{mthd.far_threshold(), 1.0};
        for (size_t i = 0; i < src_mesh.facets.size(); i++) {
            IntegralTerm<dim,R,C> term{pt, src_facet_info[i]};
            auto nearest_pt = far_near_logic.decide(pt.loc, src_facet_info[i]);
            auto integrals = mthd.compute_term(term, nearest_pt);
            for (int b = 0; b < dim; b++) {
                result[dim * i + b] = integrals[b];
            }
        }
        auto start_idx = pt_idx * n_src_dofs;
        for (size_t i = 0; i < result.size(); i++) {
            for (size_t d1 = 0; d1 < R; d1++) {
                for (size_t d2 = 0; d2 < C; d2++) {
                    auto row_idx = d1 * n_obs_dofs + pt_idx;
                    auto col_idx = d2 * n_src_dofs + i;
                    op[row_idx * C * n_src_dofs + col_idx] += result[i][d1][d2];
                }
            }
        }
    }
    return op;
}

} // END NAMESPACE tbem
#endif

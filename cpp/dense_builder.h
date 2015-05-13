#ifndef __AAAAAAAAA_BEM_H
#define __AAAAAAAAA_BEM_H

#include <cassert>
#include "dense_operator.h"
#include "obs_pt.h"
#include "util.h"
#include "integral_term.h"

namespace tbem {

template <size_t dim>
std::vector<ObsPt<dim>> setup_obs_pts(const std::vector<Vec<double,dim>>& locs,
    const std::vector<Vec<double,dim>>& normals, const std::vector<Mesh<dim>>& meshes)
{
    auto union_mesh = Mesh<dim>::create_union(meshes);
    std::vector<ObsPt<dim>> out;
    for (size_t i = 0; i < locs.size(); i++) {
        auto nf = nearest_facets(locs[i], union_mesh.facets);
        auto rich_dir = decide_richardson_dir(locs[i], nf);
        auto rich_length = hypot(rich_dir);
        out.push_back({
            rich_length / 5.0, locs[i], normals[i], rich_dir / rich_length
        });
    }
    return out;
}

template <size_t dim, size_t R, size_t C>
DenseOperator mesh_to_points_operator(const std::vector<ObsPt<dim>>& obs_pts,
    const Mesh<dim>& src_mesh, const IntegrationMethodI<dim,R,C>& mthd)
{
    auto n_obs_dofs = obs_pts.size();
    auto n_src_dofs = src_mesh.n_dofs();

    DenseOperator op(R * n_obs_dofs, C * n_src_dofs, 0.0);

    auto src_facet_info = get_facet_info(src_mesh);
// #pragma omp parallel for
    for (size_t pt_idx = 0; pt_idx < n_obs_dofs; pt_idx++) {
        const auto& pt = obs_pts[pt_idx];

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

template <size_t dim, size_t R, size_t C>
DenseOperator dense_integral_operator(const Mesh<dim>& obs_mesh,
    const Mesh<dim>& src_mesh, const IntegrationMethodI<dim,R,C>& mthd)
{
    size_t n_obs_dofs = obs_mesh.n_dofs();
    size_t n_src_dofs = src_mesh.n_dofs();

    DenseOperator op(R * n_obs_dofs, C * n_src_dofs, 0.0);

    auto src_facet_info = get_facet_info(src_mesh);
    const auto& obs_quad = mthd.get_obs_quad();

#pragma omp parallel for
    for (size_t obs_idx = 0; obs_idx < obs_mesh.facets.size(); obs_idx++) {
        auto obs_face = FacetInfo<dim>::build(obs_mesh.facets[obs_idx]);

        std::vector<Vec<Vec<Vec<double,C>,R>,dim>> row(n_src_dofs, 
                zeros<Vec<Vec<Vec<double,C>,R>,dim>>::make());
        for (size_t obs_q = 0; obs_q < obs_quad.size(); obs_q++) {
            //TODO: This is ugly and messy
            auto pt = ObsPt<dim>::away_from_nearest_facets(
                obs_quad[obs_q].x_hat, obs_face, src_mesh
            );

            const auto basis = linear_basis(obs_quad[obs_q].x_hat);
            std::vector<Vec<Vec<double,C>,R>> add_to_row(src_mesh.n_dofs());
            FarNearLogic<dim> far_near_logic{mthd.far_threshold(), 1.0};
            for (size_t i = 0; i < src_mesh.facets.size(); i++) {
                IntegralTerm<dim,R,C> term{pt, src_facet_info[i]};
                auto nearest_pt = far_near_logic.decide(pt.loc, src_facet_info[i]);
                auto integrals = mthd.compute_term(term, nearest_pt);
                for (int b = 0; b < dim; b++) {
                    add_to_row[dim * i + b] = integrals[b];
                }
            }
            for (size_t dof = 0; dof < n_src_dofs; dof++) {
                row[dof] += outer_product(basis,
                    add_to_row[dof] * obs_quad[obs_q].w * obs_face.jacobian);
            }
        }


        for (int obs_basis_idx = 0; obs_basis_idx < dim; obs_basis_idx++) {
            int obs_dof = dim * obs_idx + obs_basis_idx;
            for (size_t src_dof = 0; src_dof < n_src_dofs; src_dof++) {
                auto val_to_add = row[src_dof][obs_basis_idx];
                for (size_t d1 = 0; d1 < R; d1++) {
                    for (size_t d2 = 0; d2 < C; d2++) {
                        auto row_idx = d1 * n_obs_dofs + obs_dof;
                        auto col_idx = d2 * n_src_dofs + src_dof;
                        op[row_idx * C * n_src_dofs + col_idx] += val_to_add[d1][d2];
                    }
                }
            }
        }
    }
    return op;
}

} // END NAMESPACE tbem
#endif

#ifndef __78987898789787_MATRIX_FREE_OPERATOR_BUILDER_H
#define __78987898789787_MATRIX_FREE_OPERATOR_BUILDER_H

#include "dense_builder.h"
#include "quadrature.h"
#include "facet_info.h"
#include "sparse_operator.h"

namespace tbem {


// template <size_t n_rows, size_t n_cols> 
// void reshape_to_add(std::vector<std::vector<SparseMatrixEntry>>& entries, 
//     size_t row, size_t col, const Vec<Vec<double,n_cols>,n_rows>& data) 
// {
//     for (size_t d1 = 0; d1 < n_rows; d1++) {
//         for (size_t d2 = 0; d2 < n_cols; d2++) {
//             entries[d1 * n_cols + d2].push_back({row, col, data[d1][d2]});
//         }
//     }
// }
// 
// 
// template <size_t dim, size_t R, size_t C>
// struct MatrixFreeFarfieldOperator {
//     const Mesh<dim> src_mesh;
//     const Mesh<dim> obs_mesh;
//     const size_t n_q_obs;
//     const size_t n_q_src;
//     const BlockSparseOperator nearfield;
// 
//     MatrixFreeFarfieldOperator(const BoundaryIntegral<dim,R,C>& p, const QuadStrategy<dim>& qs):
//         p(p), qs(qs), nearfield(build_nearfield(p, qs))
//     {}
// 
//     BlockVectorX apply(const BlockVectorX& x) {
//         size_t n_obs_dofs = p.obs_mesh.n_dofs();
//         BlockVectorX farfield(R, VectorX(n_obs_dofs, 0.0));
// 
//         AdaptiveIntegrationMethod<dim,R,C> mthd(qs, p.K);
//         auto src_facet_info = get_facet_info(p.src_mesh);
// #pragma omp parallel for
//         for (size_t obs_idx = 0; obs_idx < p.obs_mesh.facets.size(); obs_idx++) {
//             auto obs_face = FacetInfo<dim>::build(p.obs_mesh.facets[obs_idx]);
//             for (size_t obs_q = 0; obs_q < qs.obs_quad.size(); obs_q++) {
//                 auto pt = ObsPt<dim>::from_face(qs.obs_quad[obs_q].x_hat, obs_face);
// 
//                 auto row = zeros<Vec<double,R>>::make();
//                 for (size_t i = 0; i < p.src_mesh.facets.size(); i++) {
//                     FarNearLogic<dim> far_near_logic{qs.far_threshold, 1.0};
//                     auto nearest_pt = far_near_logic.decide(pt.loc, src_facet_info[i]);
//                     if (nearest_pt.type != FarNearType::Farfield) {
//                         continue; 
//                     }
//                     IntegralTerm<dim,R,C> term{pt, src_facet_info[i]};
//                     auto integrals = mthd.compute_term(term, nearest_pt);
//                     for (int b = 0; b < dim; b++) {
//                         for (size_t d1 = 0; d1 < R; d1++) {
//                             for (size_t d2 = 0; d2 < C; d2++) {
//                                 row[d1] += integrals[b][d1][d2] * x[d2][dim * i + b];
//                             }
//                         }
//                     }
//                 }
// 
//                 const auto basis = linear_basis(qs.obs_quad[obs_q].x_hat);
// 
//                 for (int obs_basis_idx = 0; obs_basis_idx < dim; obs_basis_idx++) {
//                     int obs_dof = dim * obs_idx + obs_basis_idx;
//                     auto val_to_add = basis[obs_basis_idx] * row *
//                         qs.obs_quad[obs_q].w * obs_face.jacobian;
//                     for (size_t d = 0; d < R; d++) {
//                         farfield[d][obs_dof] += val_to_add[d]; 
//                     }
//                 }
//             }
//         }
// 
//         return nearfield.apply(x) + farfield;
//     }
// };
// 
// template <size_t dim, size_t R, size_t C>
// MatrixFreeFarfieldOperator<dim,R,C>
// make_matrix_free(const BoundaryIntegral<dim,R,C>& p, const QuadStrategy<dim>& qs) {
//     return MatrixFreeFarfieldOperator<dim,R,C>(p,qs);
// }

}//end namespace tbem

#endif

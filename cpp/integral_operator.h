#ifndef __ALJSDLAJKHAH_INTEGRAL_OPERATOR_H
#define __ALJSDLAJKHAH_INTEGRAL_OPERATOR_H
#include "sparse_operator.h"
#include "galerkin_operator.h"
#include "interpolation_operator.h"
#include "nbody_operator.h"

namespace tbem {

template <size_t dim, size_t R, size_t C>
BlockSparseOperator galerkin_nearfield(const Mesh<dim>& obs_mesh,
    const Mesh<dim>& src_mesh, const IntegrationMethodI<dim,R,C>& mthd) 
{
    size_t n_obs_dofs = obs_mesh.n_dofs();
    size_t n_src_dofs = src_mesh.n_dofs();
    auto n_blocks = R * C;

    auto src_facet_info = get_facet_info(src_mesh);
    std::vector<std::vector<SparseMatrixEntry>> entries(n_blocks);
    const auto& obs_quad = mthd.get_obs_quad();
#pragma omp parallel for
    for (size_t obs_idx = 0; obs_idx < obs_mesh.facets.size(); obs_idx++) {
        auto obs_face = FacetInfo<dim>::build(obs_mesh.facets[obs_idx]);
        for (size_t obs_q = 0; obs_q < obs_quad.size(); obs_q++) {
            auto pt = ObsPt<dim>::from_face(obs_quad[obs_q].x_hat, obs_face);

            std::vector<std::pair<size_t,Vec<Vec<double,C>,R>>> row;
            for (size_t i = 0; i < src_mesh.facets.size(); i++) {
                FarNearLogic<dim> far_near_logic{mthd.far_threshold(), 1.0};
                auto nearest_pt = far_near_logic.decide(pt.loc, src_facet_info[i]);
                if (nearest_pt.type == FarNearType::Farfield) {
                    continue; 
                }
                IntegralTerm<dim,R,C> term{pt, src_facet_info[i]};
                auto integrals = mthd.compute_term(term, nearest_pt);
                auto farfield_correction = -mthd.compute_farfield(term, nearest_pt);
                auto nearfield_term = integrals - farfield_correction;
                for (int b = 0; b < dim; b++) {
                    row.push_back(std::make_pair(dim * i + b, nearfield_term[b]));
                }
            }

            const auto basis = linear_basis(obs_quad[obs_q].x_hat);

            for (size_t obs_basis_idx = 0; obs_basis_idx < dim; obs_basis_idx++) {
                auto obs_dof = dim * obs_idx + obs_basis_idx;
                for (const auto& e: row) {
                    auto val_to_add = basis[obs_basis_idx] * e.second *
                        obs_quad[obs_q].w * obs_face.jacobian;
#pragma omp critical
                    for (size_t d1 = 0; d1 < R; d1++) {
                        for (size_t d2 = 0; d2 < C; d2++) {
                            entries[d1 * C + d2].push_back(
                                {obs_dof, e.first, val_to_add[d1][d2]}
                            );
                        }
                    }
                }
            }
        }
    }

    std::vector<SparseOperator> ops;
    for (size_t i = 0; i < n_blocks; i++) {
        ops.push_back(SparseOperator(n_obs_dofs, n_src_dofs, entries[i]));
    }
    return BlockSparseOperator(R, C, std::move(ops));
}


template <size_t dim, size_t R, size_t C>
struct BlockIntegralOperator: public BlockOperatorI {
    const BlockSparseOperator nearfield;
    const BlockGalerkinOperator<dim> galerkin;
    const BlockDirectNBodyOperator<dim,R,C> farfield;
    const BlockInterpolationOperator<dim> interp;

    BlockIntegralOperator(const BlockSparseOperator& nearfield,
        const BlockGalerkinOperator<dim>& galerkin,
        const BlockDirectNBodyOperator<dim,R,C>& farfield,
        const BlockInterpolationOperator<dim>& interp):
        nearfield(nearfield),
        galerkin(galerkin),
        farfield(farfield),
        interp(interp)
    {}

    virtual size_t n_block_rows() const {return nearfield.n_block_rows();}
    virtual size_t n_block_cols() const {return nearfield.n_block_cols();}
    virtual size_t n_total_rows() const {return nearfield.n_total_rows();} 
    virtual size_t n_total_cols() const {return nearfield.n_total_cols();}
    virtual BlockVectorX apply(const BlockVectorX& x) const {
        auto interpolated = interp.apply(x);
        auto nbodied = farfield.apply(interpolated);
        auto galerkin_far = galerkin.apply(nbodied);
        auto near_eval = nearfield.apply(x);

        return near_eval + galerkin_far;
    }
};


template <size_t dim, size_t R, size_t C>
BlockIntegralOperator<dim,R,C> integral_operator(const Mesh<dim>& obs_mesh,
    const Mesh<dim>& src_mesh, const IntegrationMethodI<dim,R,C>& mthd) {

    auto nearfield = galerkin_nearfield(obs_mesh, src_mesh, mthd);
    auto nbody_data = nbody_data_from_bem(obs_mesh, src_mesh,
        mthd.get_obs_quad(), mthd.get_src_quad());
    auto farfield = BlockDirectNBodyOperator<dim,R,C>{nbody_data, mthd.get_kernel()};

    auto obs_quad = mthd.get_obs_quad();
    auto src_quad = mthd.get_src_quad();
    BlockGalerkinOperator<dim> galerkin({R,C}, obs_mesh, obs_quad);
    BlockInterpolationOperator<dim> interp({R,C}, src_mesh, src_quad);

    return BlockIntegralOperator<dim,R,C>(nearfield, galerkin, farfield, interp);
}

} // end namespace tbem

#endif

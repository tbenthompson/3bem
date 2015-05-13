#ifndef __ALJSDLAJKHAH_INTEGRAL_OPERATOR_H
#define __ALJSDLAJKHAH_INTEGRAL_OPERATOR_H
#include "sparse_operator.h"
#include "galerkin_operator.h"
#include "interpolation_operator.h"
#include "nbody_operator.h"
#include "integral_term.h"

namespace tbem {

template <size_t dim, size_t R, size_t C>
using GalerkinNearfieldFnc = std::function<Vec<Vec<Vec<double,C>,R>,dim>(
    const IntegralTerm<dim,R,C>&, const NearestPoint<dim>&
    )>;

template <size_t dim, size_t R, size_t C>
SparseOperator galerkin_nearfield_helper(const Mesh<dim>& obs_mesh,
    const Mesh<dim>& src_mesh, const IntegrationMethodI<dim,R,C>& mthd,
    const GalerkinNearfieldFnc<dim,R,C>& f, const Mesh<dim>& all_mesh)
{
    size_t n_obs_dofs = obs_mesh.n_dofs();
    size_t n_src_dofs = src_mesh.n_dofs();
    auto n_blocks = R * C;

    auto src_facet_info = get_facet_info(src_mesh);
    std::vector<MatrixEntry> entries;
    const auto& obs_quad = mthd.get_obs_quad();
#pragma omp parallel for
    for (size_t obs_idx = 0; obs_idx < obs_mesh.facets.size(); obs_idx++) {
        auto obs_face = FacetInfo<dim>::build(obs_mesh.facets[obs_idx]);
        for (size_t obs_q = 0; obs_q < obs_quad.size(); obs_q++) {
            auto pt = ObsPt<dim>::away_from_nearest_facets(
                obs_quad[obs_q].x_hat, obs_face, all_mesh
            );

            std::vector<std::pair<size_t,Vec<Vec<double,C>,R>>> row;
            for (size_t i = 0; i < src_mesh.facets.size(); i++) {
                FarNearLogic<dim> far_near_logic{mthd.far_threshold(), 1.0};
                auto nearest_pt = far_near_logic.decide(pt.loc, src_facet_info[i]);
                if (nearest_pt.type == FarNearType::Farfield) {
                    continue; 
                }
                IntegralTerm<dim,R,C> term{pt, src_facet_info[i]};
                auto eval = f(term, nearest_pt);
                for (int b = 0; b < dim; b++) {
                    row.push_back(std::make_pair(dim * i + b, eval[b]));
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
                            entries.push_back({
                                d1 * n_obs_dofs + obs_dof,
                                d2 * n_src_dofs + e.first,
                                val_to_add[d1][d2]
                            });
                        }
                    }
                }
            }
        }
    }

    return SparseOperator::csr_from_coo(R * n_obs_dofs, C * n_src_dofs, entries);
}

template <size_t dim, size_t R, size_t C>
SparseOperator galerkin_nearfield(const Mesh<dim>& obs_mesh,
    const Mesh<dim>& src_mesh, const IntegrationMethodI<dim,R,C>& mthd,
    const Mesh<dim>& all_mesh) 
{
    GalerkinNearfieldFnc<dim,R,C> f = 
        [&] (const IntegralTerm<dim,R,C>& term, const NearestPoint<dim>& pt) {
            return mthd.compute_term(term, pt); 
        };
    return galerkin_nearfield_helper(obs_mesh, src_mesh, mthd, f, all_mesh);
}

template <size_t dim, size_t R, size_t C>
SparseOperator galerkin_farfield_correction(const Mesh<dim>& obs_mesh,
    const Mesh<dim>& src_mesh, const IntegrationMethodI<dim,R,C>& mthd,
    const Mesh<dim>& all_mesh) 
{
    GalerkinNearfieldFnc<dim,R,C> f = 
        [&] (const IntegralTerm<dim,R,C>& term, const NearestPoint<dim>& pt) {
            return -mthd.compute_farfield(term, pt); 
        };
    return galerkin_nearfield_helper(obs_mesh, src_mesh, mthd, f, all_mesh);
}

template <size_t dim, size_t R, size_t C>
struct IntegralOperator: public OperatorI {
    const SparseOperator nearfield;
    const SparseOperator farfield_correction;
    const GalerkinOperator<dim> galerkin;
    const BlockDirectNBodyOperator<dim,R,C> farfield;
    const InterpolationOperator<dim> interp;

    IntegralOperator(const SparseOperator& nearfield,
        const SparseOperator& farfield_correction,
        const GalerkinOperator<dim>& galerkin,
        const BlockDirectNBodyOperator<dim,R,C>& farfield,
        const InterpolationOperator<dim>& interp):
        nearfield(nearfield),
        farfield_correction(farfield_correction),
        galerkin(galerkin),
        farfield(std::move(farfield)),
        interp(interp)
    {}

    virtual size_t n_rows() const {return nearfield.n_rows();} 
    virtual size_t n_cols() const {return nearfield.n_cols();}
    virtual std::vector<double> apply(const std::vector<double>& x) const {
        auto interpolated = interp.apply(x);
        auto nbodied = farfield.apply(interpolated);
        auto galerkin_far = galerkin.apply(nbodied);
        auto eval = nearfield.apply(x);
        auto correction = farfield_correction.apply(x);
        for (size_t i = 0; i < eval.size(); i++) {
            eval[i] += galerkin_far[i] + correction[i];
        }
        return eval;
    }

    const SparseOperator& get_nearfield_matrix() {
        return nearfield;
    }
};


template <size_t dim, size_t R, size_t C>
IntegralOperator<dim,R,C> integral_operator(const Mesh<dim>& obs_mesh,
    const Mesh<dim>& src_mesh, const IntegrationMethodI<dim,R,C>& mthd,
    const Mesh<dim>& all_mesh) 
{
    auto nearfield = galerkin_nearfield(obs_mesh, src_mesh, mthd, all_mesh);
    auto far_correction = galerkin_farfield_correction(
        obs_mesh, src_mesh, mthd, all_mesh
    );
    auto nbody_data = nbody_data_from_bem(obs_mesh, src_mesh,
        mthd.get_obs_quad(), mthd.get_src_quad());
    auto farfield = BlockDirectNBodyOperator<dim,R,C>{nbody_data, mthd.get_kernel()};

    auto obs_quad = mthd.get_obs_quad();
    auto src_quad = mthd.get_src_quad();
    GalerkinOperator<dim> galerkin({R,C}, obs_mesh, obs_quad);
    InterpolationOperator<dim> interp({R,C}, src_mesh, src_quad);

    return IntegralOperator<dim,R,C>(
        nearfield, far_correction, galerkin, farfield, interp
    );
}

} // end namespace tbem

#endif

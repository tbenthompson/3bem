#ifndef TBEMALJSDLAJKHAH_INTEGRAL_OPERATOR_H
#define TBEMALJSDLAJKHAH_INTEGRAL_OPERATOR_H
#include "sparse_operator.h"
#include "galerkin_operator.h"
#include "interpolation_operator.h"
#include "nbody_operator.h"
#include "integral_term.h"
#include "nearfield_operator.h"
#include "fmm.h"

namespace tbem {

/* This functions in this file creates integral operators for galerkin boundary
 * integral equations
 *
 * The general form is
 * A = G(FI + N + C)
 * where 
 * G is the galerkin op,
 * F is the farfield nbody op,
 * I is the farfield interpolation op,
 * N is the nearfield op,
 * C is the nearfield correction op (to remove overlap with farfield)
 */

template <size_t dim, size_t R, size_t C>
struct IntegralOperator: public OperatorI {
    const SparseOperator nearfield;
    const SparseOperator farfield_correction;
    const SparseOperator galerkin;
    const std::shared_ptr<OperatorI> farfield;
    const SparseOperator interp;

    IntegralOperator(const SparseOperator& nearfield,
        const SparseOperator& farfield_correction,
        const SparseOperator& galerkin,
        const std::shared_ptr<OperatorI>& farfield,
        const SparseOperator& interp):
        nearfield(nearfield),
        farfield_correction(farfield_correction),
        galerkin(galerkin),
        farfield(farfield),
        interp(interp)
    {}

    virtual size_t n_rows() const {return galerkin.n_rows();} 
    virtual size_t n_cols() const {return nearfield.n_cols();}
    virtual std::vector<double> apply(const std::vector<double>& x) const {
        auto nbody_far = galerkin.apply(farfield->apply(interp.apply(x)));
        auto eval = nearfield.apply(x);
        auto correction = farfield_correction.apply(x);
        for (size_t i = 0; i < eval.size(); i++) {
            eval[i] += nbody_far[i] + correction[i];
        }
        return eval;
    }
};

//TODO Lots of ugly duplication in this file
template <size_t dim, size_t R, size_t C>
IntegralOperator<dim,R,C> boundary_operator(const Mesh<dim>& obs_mesh,
    const Mesh<dim>& src_mesh, const IntegrationStrategy<dim,R,C>& mthd,
    const FMMConfig& fmm_config, const Mesh<dim>& all_mesh) 
{
    (void)fmm_config;
    auto near_obs_pts = galerkin_obs_pts(obs_mesh, mthd.obs_near_quad, all_mesh);
    auto nearfield = make_nearfield_operator(near_obs_pts, src_mesh, mthd);
    auto far_obs_pts = galerkin_obs_pts(obs_mesh, mthd.obs_far_quad, all_mesh);
    auto far_correction = make_farfield_correction_operator(far_obs_pts, src_mesh, mthd);
    auto near_galerkin = make_galerkin_operator(R, obs_mesh, mthd.obs_near_quad);

    auto nbody_data = nbody_data_from_bem(
        obs_mesh, src_mesh, mthd.obs_far_quad, mthd.src_far_quad
    );
    // auto farfield = std::make_shared<FMMOperator<dim,R,C>>(
    //     FMMOperator<dim,R,C>(*mthd.K, nbody_data, fmm_config)
    // );
    auto farfield = std::make_shared<DenseOperator>(
        make_direct_nbody_operator(nbody_data, *mthd.K)
    );
    std::shared_ptr<OperatorI> farfield_ptr = farfield;

    auto far_galerkin = make_galerkin_operator(R, obs_mesh, mthd.obs_far_quad);
    auto interp = make_interpolation_operator(C, src_mesh, mthd.src_far_quad);

    return IntegralOperator<dim,R,C>(
        near_galerkin.right_multiply(nearfield),
        far_galerkin.right_multiply(far_correction),
        far_galerkin, farfield_ptr, interp
    );
}

template <size_t dim, size_t R, size_t C>
DenseOperator dense_boundary_operator(const Mesh<dim>& obs_mesh,
    const Mesh<dim>& src_mesh, const IntegrationStrategy<dim,R,C>& mthd,
    const Mesh<dim>& all_mesh)
{

    auto interp = make_interpolation_operator(C, src_mesh, mthd.src_far_quad);
    auto nbody_data = nbody_data_from_bem(
        obs_mesh, src_mesh, mthd.obs_far_quad, mthd.src_far_quad
    );
    auto farfield = make_direct_nbody_operator(nbody_data, *mthd.K);
    auto far_obs_pts = galerkin_obs_pts(obs_mesh, mthd.obs_far_quad, all_mesh);
    auto far_correction = make_farfield_correction_operator(far_obs_pts, src_mesh, mthd);
    auto far_galerkin = make_galerkin_operator(R, obs_mesh, mthd.obs_far_quad);

    auto near_obs_pts = galerkin_obs_pts(obs_mesh, mthd.obs_near_quad, all_mesh);
    auto nearfield = make_nearfield_operator(near_obs_pts, src_mesh, mthd);
    auto near_galerkin = make_galerkin_operator(R, obs_mesh, mthd.obs_near_quad);

    auto out = far_galerkin.right_multiply_with_dense(
            far_correction.add_with_dense(
                interp.left_multiply_with_dense(
                    farfield
                )
            )
        ).add(
            near_galerkin.right_multiply_with_dense(
                nearfield.to_dense()
            )
        );
    return out;
}

template <size_t dim, size_t R, size_t C>
DenseOperator dense_interior_operator(const std::vector<Vec<double,dim>>& locs,
    const std::vector<Vec<double,dim>>& normals, const Mesh<dim>& src_mesh,
    const IntegrationStrategy<dim,R,C>& mthd, const Mesh<dim>& all_mesh)
{
    auto obs_pts = interior_obs_pts(locs, normals, all_mesh);
    auto nearfield = make_nearfield_operator(obs_pts, src_mesh, mthd);
    auto far_correction = make_farfield_correction_operator(obs_pts, src_mesh, mthd);

    auto nbody_src = nbody_src_from_bem(src_mesh, mthd.src_far_quad);
    NBodyData<dim> nbody_data{
        locs, normals, nbody_src.locs, nbody_src.normals, nbody_src.weights
    };
    auto farfield = make_direct_nbody_operator(nbody_data, *mthd.K);
    auto interp = make_interpolation_operator(C, src_mesh, mthd.src_far_quad);

    auto out = far_correction.add_with_dense(
        nearfield.add_with_dense(
            interp.left_multiply_with_dense(
                farfield
            )
        )
    );
    return out;
}

} // end namespace tbem

#endif

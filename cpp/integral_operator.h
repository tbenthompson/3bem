#ifndef __ALJSDLAJKHAH_INTEGRAL_OPERATOR_H
#define __ALJSDLAJKHAH_INTEGRAL_OPERATOR_H
#include "sparse_operator.h"
#include "galerkin_operator.h"
#include "interpolation_operator.h"
#include "nbody_operator.h"
#include "integral_term.h"
#include "nearfield_operator.h"

namespace tbem {

struct IntegralOperator: public OperatorI {
    const SparseOperator nearfield;
    const SparseOperator farfield_correction;
    const SparseOperator galerkin;
    const DenseOperator farfield;
    const SparseOperator interp;

    IntegralOperator(const SparseOperator& nearfield,
        const SparseOperator& farfield_correction,
        const SparseOperator& galerkin,
        const DenseOperator& farfield,
        const SparseOperator& interp):
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
IntegralOperator integral_operator(const Mesh<dim>& obs_mesh,
    const Mesh<dim>& src_mesh, const IntegrationMethodI<dim,R,C>& mthd,
    const Mesh<dim>& all_mesh) 
{
    auto nearfield = galerkin_nearfield(obs_mesh, src_mesh, mthd, all_mesh);
    auto far_correction = galerkin_farfield_correction(
        obs_mesh, src_mesh, mthd, all_mesh
    );
    auto nbody_data = nbody_data_from_bem(obs_mesh, src_mesh,
        mthd.get_obs_quad(), mthd.get_src_quad());
    auto farfield = make_direct_nbody_operator(nbody_data, mthd.get_kernel());

    auto obs_quad = mthd.get_obs_quad();
    auto src_quad = mthd.get_src_quad();
    auto galerkin = make_galerkin_operator(R, obs_mesh, obs_quad);
    auto interp = make_interpolation_operator(C, src_mesh, src_quad);

    return IntegralOperator(
        nearfield, far_correction, galerkin, farfield, interp
    );
}

} // end namespace tbem

#endif

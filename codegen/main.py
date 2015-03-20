from mako.template import Template
import string
import random
import os
header_file = \
'''
#ifndef ${unique_fileid}
#define ${unique_fileid}

#include "dense_operator.h"
#include "block_operator.h"

namespace tbem {

template <size_t dim> struct Mesh;
BlockDenseOperator mesh_to_mesh(const Mesh<3>& obs_mesh, const Mesh<3>& src_mesh);

} //end namespace tbem

#endif
'''

cpp_file = \
'''
#include "${header_filename}"
#include "integral_term.h"
#include "dense_operator.h"
#include "mesh.h"
#include "adaptive_quad.h"
#include "richardson.h"
#include "facet_info.h"

namespace tbem {

struct ${name} {
    Vec<double,3>
    compute_singular(const IntegralTerm<3,1,1>&, const NearestPoint<3>&) const;

    Vec<double,3>
    compute_nearfield(const IntegralTerm<3,1,1>&, const NearestPoint<3>&) const;

    Vec<double,3>
    compute_farfield(const IntegralTerm<3,1,1>&, const NearestPoint<3>&) const;

    Vec<double,3>
    compute_term(const IntegralTerm<3,1,1>& term,
        const NearestPoint<3>& nearest_pt) const
    {
        switch (nearest_pt.type) {
            case FarNearType::Singular:
                return compute_singular(term, nearest_pt);
                break;
            case FarNearType::Nearfield:
                return compute_nearfield(term, nearest_pt);
                break;
            case FarNearType::Farfield:
                return compute_farfield(term, nearest_pt);
                break;
        }
        throw std::exception();
    }
};

Vec<double,3> eval_pt(
    const ObsPt<3>& obs,
    const FacetInfo<3>& src_face,
    const Vec<double,2>& x_hat,
    const Vec<double,3>& moved_obs_loc)
{
    const auto src_pt = ref_to_real(x_hat, src_face.face);
    const auto d = src_pt - moved_obs_loc;
    const auto r2 = dot_product(d, d);
    const double kernel_val = 1.0 / (4.0 * M_PI * std::sqrt(r2));
    return outer_product(linear_basis(x_hat), kernel_val * src_face.jacobian);
}

Vec<double,3>
adaptively(const IntegralTerm<3,1,1>& term, const NearestPoint<3>& near_pt,
    const Vec<double,3>& nf_obs_pt)
{
    return adaptive_integrate<Vec<double,3>>(
        [&] (double x_hat) {
            if (x_hat == 1.0) {
                return zeros<Vec<double,3>>::make();
            }
            return adaptive_integrate<Vec<double,3>>(
                [&] (double y_hat) {
                    Vec<double,2> q_pt = {x_hat, y_hat};
                    return eval_pt(term.obs, term.src_face, q_pt, nf_obs_pt);
                }, 0.0, 1 - x_hat, 1e-4);
        }, 0.0, 1.0, 1e-4);
}

Vec<double,3>
${name}::compute_singular(const IntegralTerm<3,1,1>& term,
    const NearestPoint<3>&) const
{
    std::vector<Vec<double,3>> near_steps(8);

    const double safe_dist_ratio = 5.0;
    for (int step_idx = 0; step_idx < 8; step_idx++) {
        double step_size = std::pow(2.0, -step_idx);
        double step_distance = safe_dist_ratio * term.obs.len_scale * step_size;
        auto step_loc = term.obs.loc + step_distance * term.obs.richardson_dir;
        auto shifted_near_pt = FarNearLogic<3>{3.0, 1.0}
            .decide(step_loc, term.src_face);
        near_steps[step_idx] = adaptively(term, shifted_near_pt, step_loc);
    }

    return richardson_limit(2, near_steps);
}

Vec<double,3>
${name}::compute_nearfield(const IntegralTerm<3,1,1>& term,
    const NearestPoint<3>& nearest_pt) const
{
    return adaptively(term, nearest_pt, term.obs.loc);
}

static auto G = tri_gauss(2);
Vec<double,3>
${name}::compute_farfield(const IntegralTerm<3,1,1>& term,
    const NearestPoint<3>& nearest_pt) const
{
    auto integrals = zeros<Vec<double,3>>::make();
    for (size_t i = 0; i < G.size(); i++) {
        integrals +=
            eval_pt(term.obs, term.src_face, G[i].x_hat, term.obs.loc) * G[i].w;
    }
    return integrals;
}

static auto G_obs = tri_gauss(3);
BlockDenseOperator mesh_to_mesh(const Mesh<3>& obs_mesh, const Mesh<3>& src_mesh) {
    size_t n_obs_dofs = obs_mesh.n_dofs();
    size_t n_src_dofs = src_mesh.n_dofs();
    FarNearLogic<3> far_near_logic{3.0, 1.0};

    BlockDenseOperator block_op{1, 1, {DenseOperator(n_obs_dofs, n_src_dofs, 0.0)}};

    auto src_facet_info = get_facet_info(src_mesh);
    std::cout << "HI" << std::endl;
#pragma omp parallel for
    for (size_t obs_idx = 0; obs_idx < obs_mesh.facets.size(); obs_idx++) {
        auto obs_face = FacetInfo<3>::build(obs_mesh.facets[obs_idx]);

        std::vector<Vec<double,3>> row(n_src_dofs,
                zeros<Vec<double,3>>::make());
        for (size_t obs_q = 0; obs_q < G_obs.size(); obs_q++) {
            auto obs = ObsPt<3>::from_face(G_obs[obs_q].x_hat, obs_face);

            const auto basis = linear_basis(G_obs[obs_q].x_hat);

            ${name} mthd;
            for (size_t i = 0; i < src_mesh.facets.size(); i++) {
                IntegralTerm<3,1,1> term{obs, src_facet_info[i]};
                auto nearest_pt = far_near_logic.decide(obs.loc, src_facet_info[i]);
                auto integrals = mthd.compute_term(term, nearest_pt);
                for (int b = 0; b < 3; b++) {
                    row[3 * i + b] += outer_product(basis,
                        integrals[b] * G_obs[obs_q].w * obs_face.jacobian);
                }
            }
        }


        for (int obs_basis_idx = 0; obs_basis_idx < 3; obs_basis_idx++) {
            int obs_dof = 3 * obs_idx + obs_basis_idx;
            for (size_t src_dof = 0; src_dof < n_src_dofs; src_dof++) {
                block_op.ops[0][obs_dof * n_src_dofs + src_dof] +=
                    row[src_dof][obs_basis_idx];
            }
        }
    }
    return block_op;
}

} // end namespace tbem
'''

def random_fileid(N):
    options = string.ascii_uppercase + string.digits
    return '__' + ''.join(random.choice(options) for _ in range(N))

def header_filename(fileroot):
    return fileroot + '.h'

def cpp_filename(fileroot):
    return fileroot + '.cpp'

def save_file(filepath, text):
    with open(filepath, 'w') as f:
        f.write(text)

def build_header(code_dir, fileroot):
    params = dict()
    params['unique_fileid'] = random_fileid(20) + '_' + fileroot.upper() + '_H'
    params['contents'] = ''
    params['name'] = 'CGIntegrationMethod'
    header_text = Template(header_file).render(**params)
    header_filepath = os.path.join(code_dir, header_filename(fileroot))
    save_file(header_filepath, header_text)

def build_cpp(code_dir, fileroot):
    params = dict()
    params['header_filename'] = header_filename(fileroot)
    params['name'] = 'CGIntegrationMethod'
    text = Template(cpp_file).render(**params)
    cpp_filepath = os.path.join(code_dir, cpp_filename(fileroot))
    save_file(cpp_filepath, text)

def generate():
    code_dir = '3bem'
    fileroot = 'cg_integration'
    build_header(code_dir, fileroot)
    build_cpp(code_dir, fileroot)

if __name__ == '__main__':
    main()

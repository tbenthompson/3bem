#include "fmm.h"

namespace tbem {

template <>
TranslationSurface<2> make_surrounding_surface<2>(size_t expansion_order) 
{
    std::vector<Vec<double,2>> pts;
    std::vector<Vec<double,2>> normals;

    auto theta = linspace(0, 2 * M_PI, expansion_order + 1);
    for (size_t i = 0; i < expansion_order; i++) {
        Vec<double,2> pt{std::cos(theta[i]), std::sin(theta[i])}; 
        Vec<double,2> normal = -pt;
        pts.push_back(pt);
    }

    return {pts, pts};
}

TranslationSurface<3> surrounding_surface_sphere(size_t expansion_order)
{
    std::vector<Vec<double,3>> pts;
    double a = 4 * M_PI / expansion_order;
    double d = std::sqrt(a);
    auto M_theta = static_cast<size_t>(std::round(M_PI / d));
    double d_theta = M_PI / M_theta;
    double d_phi = a / d_theta;
    for (size_t m = 0; m < M_theta; m++) {
        double theta = M_PI * (m + 0.5) / M_theta;
        auto M_phi = static_cast<size_t>(
            std::round(2 * M_PI * std::sin(theta) / d_phi)
        );
        for (size_t n = 0; n < M_phi; n++) {
            double phi = 2 * M_PI * n / M_phi;
            double x = std::sin(theta) * std::cos(phi);
            double y = std::sin(theta) * std::sin(phi);
            double z = std::cos(theta);
            pts.push_back({x, y, z});
        }
    }

    return {pts, pts};
}

template <>
TranslationSurface<3> make_surrounding_surface<3>(size_t expansion_order) 
{
    return surrounding_surface_sphere(expansion_order);
}

template <size_t dim>
std::vector<Vec<double,dim>> 
TranslationSurface<dim>::move(const Box<dim>& box, double r_ref) const
{
    auto cell_r = hypot(box.half_width);
    auto out_r = r_ref * cell_r;
    auto out_center = box.center;

    std::vector<Vec<double,dim>> out_pts(pts.size());
    for (size_t i = 0; i < pts.size(); i++) {
        out_pts[i] = out_r * pts[i] + out_center;
    }

    return out_pts;
}

template struct TranslationSurface<2>;
template struct TranslationSurface<3>;

template <size_t dim, size_t R, size_t C>
std::vector<double>
nbody_matrix(const Kernel<dim,R,C>& K, const NBodyData<dim>& data) 
{
    auto n_pairs = data.obs_locs.size() * data.src_locs.size();
    auto n_blocks = R * C;
    std::vector<double> op(n_pairs * n_blocks);

    for (size_t i = 0; i < data.obs_locs.size(); i++) {
        for (size_t j = 0; j < data.src_locs.size(); j++) {
            auto kernel_val = K(
                data.obs_locs[i], data.src_locs[j],
                data.obs_normals[i], data.src_normals[j]
            );

            for (size_t d1 = 0; d1 < R; d1++) {
                auto row = d1 * data.obs_locs.size() + i;
                for (size_t d2 = 0; d2 < C; d2++) {
                    auto col = d2 * data.src_locs.size() + j;
                    auto matrix_idx = row * C * data.src_locs.size() + col;
                    op[matrix_idx] = kernel_val[d1][d2];
                }
            }
        }
    }

    return op;
}

template <size_t dim, size_t R, size_t C>
std::vector<double>
nbody_eval(const Kernel<dim,R,C>& K, const NBodyData<dim>& data,
           const std::vector<double>& x) 
{
    std::vector<double> out(R * data.obs_locs.size(), 0.0);

    for (size_t i = 0; i < data.obs_locs.size(); i++) {
        for (size_t j = 0; j < data.src_locs.size(); j++) {
            auto kernel_val = data.src_weights[j] * K(
                data.obs_locs[i], data.src_locs[j],
                data.obs_normals[i], data.src_normals[j]
            );

            for (size_t d1 = 0; d1 < R; d1++) {
                auto row = d1 * data.obs_locs.size() + i;
                for (size_t d2 = 0; d2 < C; d2++) {
                    auto col = d2 * data.src_locs.size() + j;
                    out[row] += kernel_val[d1][d2] * x[col];
                }
            }
        }
    }

    return out;
}

template <size_t dim, size_t R, size_t C>
FMMOperator<dim,R,C>::FMMOperator(const Kernel<dim,R,C>& K,
    const NBodyData<dim>& data, const FMMConfig& config):
    K(K),
    data(data),
    surface(make_surrounding_surface<dim>(config.expansion_order)),
    src_oct(build_octree(data.src_locs, config.min_pts_per_cell)),
    obs_oct(build_octree(data.obs_locs, config.min_pts_per_cell)),
    config(config)
{}

template <size_t dim, size_t R, size_t C>
void FMMOperator<dim,R,C>::build_check_to_equiv(const Octree<dim>& cell, 
    CheckToEquiv& created_ops) const
{
    if (created_ops.size() == cell.data.level) {
        auto check_pts = surface.upward_check_points(cell.data.bounds, config.d);
        auto equiv_pts = surface.upward_equiv_points(cell.data.bounds, config.d);

        auto op = nbody_matrix(K, {
            check_pts, surface.normals,
            equiv_pts, surface.normals, {}
        });

        auto svd = svd_decompose(op);
        //TODO: SVD THRESHOLDING!
        set_threshold(svd, 1e-10);
        auto cond = condition_number(svd);
        std::cout << cond << std::endl;
        created_ops.push_back(std::move(svd));
    }

    for (const auto& c: cell.children) {
        if (c == nullptr) {
            continue;
        }
        build_check_to_equiv(*c, created_ops);
    }
}

template <size_t dim, size_t R, size_t C>
std::vector<double> 
FMMOperator<dim,R,C>::apply_src_to_check(const Octree<dim>& cell,
    const BlockVectorX& x) const
{
    assert(cell.is_leaf());
    auto n_src = cell.data.indices.size();

    NBodyData<dim> s2c;
    s2c.src_locs.resize(n_src);
    s2c.src_normals.resize(n_src);
    s2c.src_weights.resize(n_src);
    s2c.obs_locs = surface.upward_check_points(cell.data.bounds, config.d);
    s2c.obs_normals = surface.normals;

    std::vector<double> src_str(n_src * C);
    for (size_t i = 0; i < n_src; i++) {
        s2c.src_locs[i] = data.src_locs[cell.data.indices[i]];
        s2c.src_normals[i] = data.src_normals[cell.data.indices[i]];
        s2c.src_weights[i] = data.src_weights[cell.data.indices[i]];
        for (size_t d = 0; d < C; d++) {
            src_str[d * n_src + i] = x[d][cell.data.indices[i]];
        }
    }

    return nbody_eval(K, s2c, src_str);
}

template <size_t dim, size_t R, size_t C>
std::vector<double>
FMMOperator<dim,R,C>::apply_children_to_check(const Octree<dim>& cell,
    const typename P2MData<dim>::ChildrenType& child_p2m) const
{
    assert(!cell.is_leaf());

    auto check_pts = surface.upward_check_points(cell.data.bounds, config.d);
    auto n_check = check_pts.size();
    auto n_equiv = n_check;

    auto n_children = cell.count_children();
    auto n_src = n_children * n_check;

    NBodyData<dim> c2c;
    c2c.obs_locs = check_pts;
    c2c.obs_normals = surface.normals;
    c2c.src_locs.resize(n_src);
    c2c.src_normals.resize(n_src);
    c2c.src_weights.resize(n_src);

    std::vector<double> src_str(n_src * C);
    size_t child_idx = 0;
    for (size_t c = 0; c < Octree<dim>::split; c++) {
        auto& child = cell.children[c];
        if (child == nullptr) {
            continue;
        }
        auto equiv_pts = surface.upward_equiv_points(child->data.bounds, config.d);
        for (size_t i = 0; i < n_equiv; i++) {
            auto idx = child_idx * n_equiv + i;
            c2c.src_locs[idx] = equiv_pts[i];
            c2c.src_normals[idx] = surface.normals[i];
            c2c.src_weights[idx] = 1.0;
            for (size_t d = 0; d < C; d++) {
                src_str[d * n_src + idx] = child_p2m[c]->data[d * n_equiv + i];
            }
        }
        child_idx++;
    }

    return nbody_eval(K, c2c, src_str);
}

template <size_t dim, size_t R, size_t C>
std::unique_ptr<P2MData<dim>> 
FMMOperator<dim,R,C>::P2M(const Octree<dim>& cell, const BlockVectorX& x,
    CheckToEquiv& check_to_equiv_ops) const
{
    std::vector<double> check_eval;
    typename P2MData<dim>::ChildrenType child_P2M;

    if (cell.is_leaf()) {
        //leaf P2M
        check_eval = apply_src_to_check(cell, x);
    } else {
        //recurse and then M2M
#pragma omp parallel for if(cell.data.level <= 0)
        for (size_t i = 0; i < Octree<dim>::split; i++) {
            if (cell.children[i] == nullptr) {
                continue;
            }
            child_P2M[i] = P2M(*cell.children[i], x, check_to_equiv_ops);
        }
        check_eval = apply_children_to_check(cell, child_P2M);
    }

    auto& svd = check_to_equiv_ops[cell.data.level];
    auto equiv_srcs = svd_solve(svd, check_eval);

    return std::unique_ptr<P2MData<dim>>(new P2MData<dim>{
        equiv_srcs, std::move(child_P2M)
    });
}

template <size_t dim, size_t R, size_t C>
void FMMOperator<dim,R,C>::M2P(BlockVectorX& out, const Octree<dim>& obs_cell,
    const Octree<dim>& src_cell, const P2MData<dim>& p2m) const 
{
    auto equiv_pts = surface.upward_equiv_points(src_cell.data.bounds, config.d);
    auto n_equiv = equiv_pts.size();
    NBodyData<dim> m2p;
    m2p.src_locs = equiv_pts;
    m2p.src_normals = surface.normals;
    m2p.src_weights = std::vector<double>(equiv_pts.size(), 1.0);
    auto n_obs = obs_cell.data.indices.size();
    m2p.obs_locs.resize(n_obs);
    m2p.obs_normals.resize(n_obs);
    for (size_t i = 0; i < n_obs; i++) {
        m2p.obs_locs[i] = data.obs_locs[obs_cell.data.indices[i]];
        m2p.obs_normals[i] = data.obs_normals[obs_cell.data.indices[i]];
    }

    auto res = nbody_eval(K, m2p, p2m.data);

    for (size_t i = 0; i < n_obs; i++) {
        for (size_t d = 0; d < R; d++) {
            out[d][obs_cell.data.indices[i]] += res[d * n_obs + i];
        }
    }
}

template <size_t dim, size_t R, size_t C>
void FMMOperator<dim,R,C>::P2P(BlockVectorX& out, const Octree<dim>& obs_cell,
    const Octree<dim>& src_cell, const BlockVectorX& x) const 
{
    auto n_src = src_cell.data.indices.size();
    auto n_obs = obs_cell.data.indices.size();

    NBodyData<dim> p2p;
    p2p.src_locs.resize(n_src);
    p2p.src_normals.resize(n_src);
    p2p.src_weights.resize(n_src);
    std::vector<double> src_str(C * n_src);
    for (size_t j = 0; j < n_src; j++) {
        p2p.src_locs[j] = data.src_locs[src_cell.data.indices[j]];
        p2p.src_normals[j] = data.src_normals[src_cell.data.indices[j]];
        p2p.src_weights[j] = data.src_weights[src_cell.data.indices[j]];
        for (size_t d = 0; d < C; d++) {
            src_str[d * n_src + j] = x[d][src_cell.data.indices[j]];
        }
    }

    p2p.obs_locs.resize(n_obs);
    p2p.obs_normals.resize(n_obs);
    for (size_t i = 0; i < n_obs; i++) {
        p2p.obs_locs[i] = data.obs_locs[obs_cell.data.indices[i]];
        p2p.obs_normals[i] = data.obs_normals[obs_cell.data.indices[i]];
    }

    auto res = nbody_eval(K, p2p, src_str);

    for (size_t i = 0; i < n_obs; i++) {
        for (size_t d = 0; d < C; d++) {
            out[d][obs_cell.data.indices[i]] += res[d * n_obs + i];
        }
    }
}

template struct FMMOperator<2,1,1>;
template struct FMMOperator<2,2,2>;
template struct FMMOperator<3,1,1>;
template struct FMMOperator<3,3,3>;

} // END namespace tbem

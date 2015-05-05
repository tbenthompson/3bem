#include "fmm.h"

//TODO: Remove the timing code.
#undef TIC
#define TIC
#undef TIC2
#define TIC2
#undef TOC
#define TOC

namespace tbem {

template <>
TranslationSurface<2> make_surrounding_surface<2>(size_t order) 
{
    std::vector<Vec<double,2>> pts;
    std::vector<Vec<double,2>> normals;

    auto theta = linspace(0, 2 * M_PI, order + 1);
    for (size_t i = 0; i < order; i++) {
        Vec<double,2> pt{std::cos(theta[i]), std::sin(theta[i])}; 
        Vec<double,2> normal = -pt;
        pts.push_back(pt);
    }

    return {pts, pts};
}

TranslationSurface<3> surrounding_surface_sphere(size_t order)
{
    std::vector<Vec<double,3>> pts;
    double a = 4 * M_PI / order;
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
TranslationSurface<3> make_surrounding_surface<3>(size_t order) 
{
    return surrounding_surface_sphere(order);
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
           double const* x) 
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
    surface(make_surrounding_surface<dim>(config.order)),
    src_oct(build_octree(data.src_locs, config.min_pts_per_cell)),
    obs_oct(build_octree(data.obs_locs, config.min_pts_per_cell)),
    config(config)
{}


template <size_t dim, size_t R, size_t C>
SVDPtr svd_inverse_nbody(const Kernel<dim,R,C>& K, const NBodyData<dim>& data) 
{
    auto op = nbody_matrix(K, data);
    auto svd = svd_decompose(op);
    //TODO: SVD THRESHOLDING!
    // In some cases, the equivalent surface to check surface operator
    // is poorly conditioned. In this case, truncate the singular values 
    // to solve a regularized least squares version of the problem.
    set_threshold(svd, 1e-10);
    // auto cond = condition_number(svd);
    // assert(cond < 1e10);
    return std::move(svd);
}

template <size_t dim, size_t R, size_t C>
void FMMOperator<dim,R,C>::build_check_to_equiv(const Octree<dim>& cell, 
    CheckToEquiv& up_ops, CheckToEquiv& down_ops) const
{
    if (up_ops.size() == cell.data.level) {
        NBodyData<dim> up_data{
            surface.upward_check_points(cell.data.bounds, config.d), surface.normals,
            surface.upward_equiv_points(cell.data.bounds, config.d), surface.normals,
            {}
        };

        up_ops.push_back(svd_inverse_nbody(K, up_data));

        NBodyData<dim> down_data{
            surface.downward_check_points(cell.data.bounds, config.d), surface.normals,
            surface.downward_equiv_points(cell.data.bounds, config.d), surface.normals,
            {}
        };
        down_ops.push_back(svd_inverse_nbody(K, down_data));
    }

    for (const auto& c: cell.children) {
        if (c == nullptr) {
            continue;
        }
        build_check_to_equiv(*c, up_ops, down_ops);
    }
}

template <size_t dim, size_t R, size_t C>
void FMMOperator<dim,R,C>::P2M(const Octree<dim>& cell,
    const SVDPtr& check_to_equiv, const std::vector<double>& x,
    double* parent_multipoles) const
{
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
            src_str[d * n_src + i] = x[d * data.src_locs.size() + cell.data.indices[i]];
        }
    }

    auto check_eval = nbody_eval(K, s2c, src_str.data());
    auto equiv_srcs = svd_solve(check_to_equiv, check_eval);

#pragma omp critical
    for (size_t i = 0; i < R * config.order; i++) {
        parent_multipoles[i] = equiv_srcs[i];
    }
}

template <size_t dim, size_t R, size_t C>
void FMMOperator<dim,R,C>::M2M(const Octree<dim>& cell,
    const SVDPtr& check_to_equiv, std::vector<double*>& child_multipoles,
    double* parent_multipoles) const
{
    assert(!cell.is_leaf());

    auto check_pts = surface.upward_check_points(cell.data.bounds, config.d);

    auto n_children = cell.n_immediate_children();
    auto n_src = n_children * config.order;

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
        for (size_t i = 0; i < config.order; i++) {
            auto idx = child_idx * config.order + i;
            c2c.src_locs[idx] = equiv_pts[i];
            c2c.src_normals[idx] = surface.normals[i];
            c2c.src_weights[idx] = 1.0;
            for (size_t d = 0; d < C; d++) {
                src_str[d * n_src + idx] = child_multipoles[c][d * config.order + i];
            }
        }
        child_idx++;
    }

    auto check_eval = nbody_eval(K, c2c, src_str.data());
    auto equiv_srcs = svd_solve(check_to_equiv, check_eval);

#pragma omp critical
    for (size_t i = 0; i < R * config.order; i++) {
        parent_multipoles[i] = equiv_srcs[i];
    }
}


template <size_t dim, size_t R, size_t C>
void FMMOperator<dim,R,C>::P2P(const Octree<dim>& obs_cell,
    const Octree<dim>& src_cell, const std::vector<double>& x, std::vector<double>& out) const 
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
            src_str[d * n_src + j] = x[d * data.src_locs.size() + src_cell.data.indices[j]];
        }
    }

    p2p.obs_locs.resize(n_obs);
    p2p.obs_normals.resize(n_obs);
    for (size_t i = 0; i < n_obs; i++) {
        p2p.obs_locs[i] = data.obs_locs[obs_cell.data.indices[i]];
        p2p.obs_normals[i] = data.obs_normals[obs_cell.data.indices[i]];
    }

    auto res = nbody_eval(K, p2p, src_str.data());

    for (size_t i = 0; i < n_obs; i++) {
        for (size_t d = 0; d < R; d++) {
            #pragma omp atomic
            out[d * data.obs_locs.size() + obs_cell.data.indices[i]] += res[d * n_obs + i];
        }
    }
}

template <size_t dim, size_t R, size_t C>
void FMMOperator<dim,R,C>::P2L(const Octree<dim>& obs_cell,
    const Octree<dim>& src_cell, const SVDPtr& check_to_equiv,
    const std::vector<double>& x, double* locals) const
{
    //TODO: This is almost identical to the P2M operator.
    auto n_src = src_cell.data.indices.size();

    NBodyData<dim> s2c;
    s2c.src_locs.resize(n_src);
    s2c.src_normals.resize(n_src);
    s2c.src_weights.resize(n_src);

    s2c.obs_locs = surface.downward_check_points(obs_cell.data.bounds, config.d);
    s2c.obs_normals = surface.normals;

    std::vector<double> src_str(n_src * C);
    for (size_t i = 0; i < n_src; i++) {
        s2c.src_locs[i] = data.src_locs[src_cell.data.indices[i]];
        s2c.src_normals[i] = data.src_normals[src_cell.data.indices[i]];
        s2c.src_weights[i] = data.src_weights[src_cell.data.indices[i]];
        for (size_t d = 0; d < C; d++) {
            src_str[d * n_src + i] = x[d * data.src_locs.size() + src_cell.data.indices[i]];
        }
    }

    auto check_eval = nbody_eval(K, s2c, src_str.data());
    auto equiv_srcs = svd_solve(check_to_equiv, check_eval);

    for (size_t i = 0; i < R * config.order; i++) {
        #pragma omp atomic
        locals[i] += equiv_srcs[i];
    }
}

template <size_t dim, size_t R, size_t C>
void FMMOperator<dim,R,C>::M2P(const Octree<dim>& obs_cell,
    const Octree<dim>& src_cell, double* multipoles, std::vector<double>& out) const
{
    auto equiv_pts = surface.upward_equiv_points(src_cell.data.bounds, config.d);

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

    auto res = nbody_eval(K, m2p, multipoles);

    for (size_t i = 0; i < n_obs; i++) {
        for (size_t d = 0; d < R; d++) {
            #pragma omp atomic
            out[d * data.obs_locs.size() + obs_cell.data.indices[i]] += res[d * n_obs + i];
        }
    }
}

template <size_t dim, size_t R, size_t C>
void FMMOperator<dim,R,C>::M2L(const Octree<dim>& obs_cell, const Octree<dim>& src_cell,
    const SVDPtr& down_check_to_equiv, double* multipoles, 
    double* locals) const
{
    auto equiv_pts = surface.upward_equiv_points(src_cell.data.bounds, config.d);
    NBodyData<dim> m2l;
    m2l.src_locs = equiv_pts;
    m2l.src_normals = surface.normals;
    m2l.src_weights = std::vector<double>(equiv_pts.size(), 1.0);
    
    auto check_pts = surface.downward_check_points(obs_cell.data.bounds, config.d);
    m2l.obs_locs = check_pts;
    m2l.obs_normals = surface.normals;

    auto check_eval = nbody_eval(K, m2l, multipoles);
    auto solved = svd_solve(down_check_to_equiv, check_eval);

    for (size_t i = 0; i < solved.size(); i++) {
        #pragma omp atomic
        locals[i] += solved[i];
    }
}

template <size_t dim, size_t R, size_t C>
void FMMOperator<dim,R,C>::L2L(const Octree<dim>& cell,
    const SVDPtr& check_to_equiv, double* parent_locals,
    std::vector<double*>& child_locals) const
{
    assert(!cell.is_leaf());

    auto equiv_pts = surface.downward_equiv_points(cell.data.bounds, config.d);

    NBodyData<dim> l2l;
    l2l.src_locs = equiv_pts;
    l2l.src_normals = surface.normals;
    l2l.src_weights = std::vector<double>(config.order, 1.0);

    for (size_t c = 0; c < Octree<dim>::split; c++) {
        auto& child = cell.children[c];
        if (child == nullptr) {
            continue;
        }

        auto check_pts = surface.downward_check_points(child->data.bounds, config.d);
        l2l.obs_locs = check_pts;
        l2l.obs_normals = surface.normals;

        auto check_eval = nbody_eval(K, l2l, parent_locals);
        auto equiv_srcs = svd_solve(check_to_equiv, check_eval);

        for (size_t i = 0; i < R * config.order; i++) {
            #pragma omp atomic
            child_locals[c][i] += equiv_srcs[i];
        }
    }
}

template <size_t dim, size_t R, size_t C>
void FMMOperator<dim,R,C>::L2P(const Octree<dim>& cell,
    double* locals, std::vector<double>& out) const
{
    //TODO: Code is essentially identical to M2P
    auto equiv_pts = surface.downward_equiv_points(cell.data.bounds, config.d);

    NBodyData<dim> l2p;
    l2p.src_locs = equiv_pts;
    l2p.src_normals = surface.normals;
    l2p.src_weights = std::vector<double>(equiv_pts.size(), 1.0);

    auto n_obs = cell.data.indices.size();
    l2p.obs_locs.resize(n_obs);
    l2p.obs_normals.resize(n_obs);
    for (size_t i = 0; i < n_obs; i++) {
        l2p.obs_locs[i] = data.obs_locs[cell.data.indices[i]];
        l2p.obs_normals[i] = data.obs_normals[cell.data.indices[i]];
    }

    auto res = nbody_eval(K, l2p, locals);

    for (size_t i = 0; i < n_obs; i++) {
        for (size_t d = 0; d < R; d++) {
            #pragma omp atomic
            out[d * data.obs_locs.size() + cell.data.indices[i]] += res[d * n_obs + i];
        }
    }
}


template <size_t dim, size_t R, size_t C>
std::vector<double> FMMOperator<dim,R,C>::execute_tasks(const FMMTasks<dim>& tasks,
    const std::vector<double>& x, 
    const CheckToEquiv& up_check_to_equiv,
    const CheckToEquiv& down_check_to_equiv) const
{
    // std::cout << "P2M tasks: " << tasks.p2ms.size() << std::endl;
    // std::cout << "P2L tasks: " << tasks.p2ls.size() << std::endl;
    // std::cout << "P2P tasks: " << tasks.p2ps.size() << std::endl;
    // std::cout << "M2M tasks: " << tasks.m2ms.size() << std::endl;
    // std::cout << "M2L tasks: " << tasks.m2ls.size() << std::endl;
    // std::cout << "M2P tasks: " << tasks.m2ps.size() << std::endl;
    // std::cout << "L2L tasks: " << tasks.l2ls.size() << std::endl;
    // std::cout << "L2P tasks: " << tasks.l2ps.size() << std::endl;

    auto n_src_cells = 1 + src_oct.n_children();
    std::vector<double> multipoles(n_src_cells * config.order * C);
    std::vector<double> out(R * data.obs_locs.size(), 0.0);

    auto equiv_start = [&](const Octree<dim>& cell) {
        return cell.data.index * config.order * C;
    };

    TIC
#pragma omp parallel for 
    for (size_t i = 0; i < tasks.p2ms.size(); i++) {
        auto& cell = tasks.p2ms[i].cell;
        auto& check_to_equiv_op = up_check_to_equiv[cell.data.level];
        P2M(cell, check_to_equiv_op, x, multipoles.data() + equiv_start(cell));
    }
    TOC("P2M");

    TIC2
// #pragma omp parallel for 
    for (size_t i = 0; i < tasks.m2ms.size(); i++) {
        auto& cell = tasks.m2ms[i].cell;
        auto& check_to_equiv_op = up_check_to_equiv[cell.data.level];
        std::vector<double*> child_data_ptrs(Octree<dim>::split);
        for (size_t c = 0; c < Octree<dim>::split; c++) {
            if (cell.children[c] == nullptr) {
                continue;
            }
            auto ptr = multipoles.data() + equiv_start(*cell.children[c]);
            child_data_ptrs[c] = ptr;
        }
        auto parent_data_ptr = multipoles.data() + equiv_start(cell);
        M2M(cell, check_to_equiv_op, child_data_ptrs, parent_data_ptr);
    }
    TOC("M2M");

    TIC2
#pragma omp parallel for 
    for (size_t i = 0; i < tasks.p2ps.size(); i++) {
        auto t = tasks.p2ps[i];
        P2P(t.obs_cell, t.src_cell, x, out);
    }
    TOC("P2P");

    TIC2
#pragma omp parallel for
    for (size_t i = 0; i < tasks.m2ps.size(); i++) {
        auto t = tasks.m2ps[i];
        auto data_ptr = multipoles.data() + equiv_start(t.src_cell);
        M2P(t.obs_cell, t.src_cell, data_ptr, out);
    }
    TOC("M2P");

    auto n_obs_cells = 1 + obs_oct.n_children();
    std::vector<double> locals(n_obs_cells * config.order * C, 0.0);

    TIC2
#pragma omp parallel for
    for (size_t i = 0; i < tasks.p2ls.size(); i++) {
        auto t = tasks.p2ls[i];
        auto& check_to_equiv_op = down_check_to_equiv[t.obs_cell.data.level];
        auto data_ptr = locals.data() + equiv_start(t.obs_cell);
        P2L(t.obs_cell, t.src_cell, check_to_equiv_op, x, data_ptr);
    }
    TOC("P2L");

    TIC2
#pragma omp parallel for
    for (size_t i = 0; i < tasks.m2ls.size(); i++) {
        auto t = tasks.m2ls[i];
        auto& check_to_equiv_op = down_check_to_equiv[t.obs_cell.data.level];
        auto multipole_data = multipoles.data() + equiv_start(t.src_cell);
        auto local_data = locals.data() + equiv_start(t.obs_cell);
        M2L(t.obs_cell, t.src_cell, check_to_equiv_op, multipole_data, local_data);
    }
    TOC("M2L");

    TIC2
// #pragma omp parallel for
    for (size_t i = 0; i < tasks.l2ls.size(); i++) {
        auto& cell = tasks.l2ls[i].cell;
        auto& check_to_equiv_op = down_check_to_equiv[cell.data.level + 1];
        std::vector<double*> child_data_ptrs(Octree<dim>::split);
        for (size_t c = 0; c < Octree<dim>::split; c++) {
            if (cell.children[c] == nullptr) {
                continue;
            }
            auto ptr = locals.data() + equiv_start(*cell.children[c]);
            child_data_ptrs[c] = ptr;
        }
        auto parent_data_ptr = locals.data() + equiv_start(cell);
        L2L(cell, check_to_equiv_op, parent_data_ptr, child_data_ptrs);
    }
    TOC("L2L");

    TIC2
#pragma omp parallel for
    for (size_t i = 0; i < tasks.l2ps.size(); i++) {
        auto& cell = tasks.l2ps[i].cell;
        auto data_ptr = locals.data() + equiv_start(cell);
        L2P(cell, data_ptr, out);
    }
    TOC("L2P");

    return out;
}

template <size_t dim, size_t R, size_t C>
std::vector<double> FMMOperator<dim,R,C>::apply(const std::vector<double>& x) const 
{
    assert(x.size() == C * data.src_locs.size());

    TIC
    //TODO: CheckToEquivs can be class members since they only depend on the
    //tree structure.
    CheckToEquiv up_check_to_equiv;
    CheckToEquiv down_check_to_equiv;
    build_check_to_equiv(src_oct, up_check_to_equiv, down_check_to_equiv);
    TOC("Computing check to equiv pseudoinverses");

    TIC2
    FMMTasks<dim> tasks;
    upward_traversal(src_oct, tasks);
    dual_tree(obs_oct, src_oct, tasks);
    downward_traversal(obs_oct, tasks);
    TOC("Traversal/planning");

    TIC2
    auto out = execute_tasks(tasks, x, up_check_to_equiv, down_check_to_equiv);
    TOC("Execution");

    return out;
}

template <size_t dim, size_t R, size_t C>
void FMMOperator<dim,R,C>::upward_traversal(const Octree<dim>& cell,
    FMMTasks<dim>& tasks) const
{
    if (cell.is_leaf()) {
        tasks.p2ms.push_back({cell});
    } else {
        //recurse and then M2M
        for (size_t i = 0; i < Octree<dim>::split; i++) {
            if (cell.children[i] == nullptr) {
                continue;
            }
            upward_traversal(*cell.children[i], tasks);
        }
        tasks.m2ms.push_back({cell});
    }
}

template <size_t dim, size_t R, size_t C>
void FMMOperator<dim,R,C>::dual_tree(const Octree<dim>& obs_cell,
    const Octree<dim>& src_cell, FMMTasks<dim>& tasks) const
{
    auto r_src = hypot(src_cell.data.bounds.half_width);
    auto r_obs = hypot(obs_cell.data.bounds.half_width);
    double r_max = std::max(r_src, r_obs);
    double r_min = std::min(r_src, r_obs);
    auto sep = hypot(obs_cell.data.bounds.center - src_cell.data.bounds.center);
    if (r_max + config.mac * r_min <= config.mac * sep) {
        bool small_src = src_cell.data.indices.size() <= config.order;
        bool small_obs = obs_cell.data.indices.size() <= config.order;
        if (small_obs && small_src) {
            tasks.p2ps.push_back({obs_cell, src_cell});
        } else if (small_obs) {
            tasks.m2ps.push_back({obs_cell, src_cell});
        } else if (small_src) {
            tasks.p2ls.push_back({obs_cell, src_cell});
        } else {
            tasks.m2ls.push_back({obs_cell, src_cell});
        }
        return;
    }

    if (obs_cell.is_leaf() && src_cell.is_leaf()) {
        tasks.p2ps.push_back({obs_cell, src_cell});
        return;
    }

    bool src_is_shallower = obs_cell.data.level > src_cell.data.level;
    bool split_src = (src_is_shallower && !src_cell.is_leaf()) || obs_cell.is_leaf();
    if (split_src) {
        //split src because it is shallower
        for (size_t c = 0; c < Octree<dim>::split; c++) {
            if (src_cell.children[c] == nullptr) {
                continue;
            }
            dual_tree(obs_cell, *src_cell.children[c], tasks);
        }
    } else {
        //split obs
        for (size_t c = 0; c < Octree<dim>::split; c++) {
            if (obs_cell.children[c] == nullptr) {
                continue;
            }
            dual_tree(*obs_cell.children[c], src_cell, tasks);
        }
    }
}

template <size_t dim, size_t R, size_t C>
void FMMOperator<dim,R,C>::downward_traversal(const Octree<dim>& obs_cell,
    FMMTasks<dim>& tasks) const
{
    if (obs_cell.is_leaf()) {
        tasks.l2ps.push_back({obs_cell});
        return;
    } else {
        tasks.l2ls.push_back({obs_cell});
        for (size_t c = 0; c < Octree<dim>::split; c++) {
            if (obs_cell.children[c] == nullptr) {
                continue;
            }
            downward_traversal(*obs_cell.children[c], tasks);
        }
    }
}


template struct FMMOperator<2,1,1>;
template struct FMMOperator<2,2,2>;
template struct FMMOperator<3,1,1>;
template struct FMMOperator<3,3,3>;

} // END namespace tbem

#ifndef __KASDJLKJKJ_FMM_H
#define __KASDJLKJKJ_FMM_H

#include <omp.h>
#include <cassert>
#include <memory>
#include "kernel.h"
#include "nbody_data.h"
#include "octree.h"
#include "numbers.h"
#include "vectorx.h"
#include "util.h"
#include "eigen_wrapper.h"

namespace tbem {

template <size_t dim>
struct TranslationSurface {
    const std::vector<Vec<double,dim>> pts;
    const std::vector<Vec<double,dim>> normals;

    std::vector<Vec<double,dim>> move(const Box<dim>& box, double r_ref) const;

    std::vector<Vec<double,dim>> 
    upward_check_points(const Box<dim>& box, double d) const
    {
        auto r_ref = 4.0 - std::sqrt(2) - 2 * d;
        return move(box, r_ref);
    }

    std::vector<Vec<double,dim>> 
    upward_equiv_points(const Box<dim>& box, double d) const
    {
        auto r_ref = std::sqrt(2) + d;
        return move(box, r_ref);
    }
};
template <size_t dim>
TranslationSurface<dim> make_surrounding_surface(size_t expansion_order);

typedef std::vector<SVDPtr> CheckToEquiv;
template <size_t dim>
using P2MData = OctreeData<dim,std::vector<double>>;

struct FMMConfig {
    const double mac;
    // const double inner_radius;
    // const double outer_radius;
    const size_t expansion_order;
    const size_t min_pts_per_cell;
    const double d;
    // const double alpha;
};

template <size_t dim>
struct M2PTask
{
    const Octree<dim>& obs_cell;
    const Octree<dim>& src_cell;
    const P2MData<dim>& p2m;
};

template <size_t dim>
struct P2PTask
{
    const Octree<dim>& obs_cell;
    const Octree<dim>& src_cell;
};

template <size_t dim>
struct FMMTasks 
{
    std::vector<M2PTask<dim>> m2ps;
    std::vector<P2PTask<dim>> p2ps;
};

template <size_t dim, size_t R, size_t C>
struct FMMOperator {
    void dual_tree(const Octree<dim>& obs_cell,
        const Octree<dim>& src_cell, const P2MData<dim>& p2m,
        FMMTasks<dim>& tasks) const
    {
        auto r_src = hypot(src_cell.data.bounds.half_width);
        auto r_obs = hypot(obs_cell.data.bounds.half_width);
        auto sep = hypot(obs_cell.data.bounds.center - src_cell.data.bounds.center);
        if (sep - r_obs > config.mac * r_src) {
            tasks.m2ps.push_back({obs_cell, src_cell, p2m});
            return;
        }

        if (src_cell.is_leaf() && obs_cell.is_leaf()) {
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
                dual_tree(obs_cell, *src_cell.children[c], *p2m.children[c], tasks);
            }
        } else {
            //split obs
            for (size_t c = 0; c < Octree<dim>::split; c++) {
                if (obs_cell.children[c] == nullptr) {
                    continue;
                }
                dual_tree(*obs_cell.children[c], src_cell, p2m, tasks);
            }
        }
    }

    BlockVectorX execute_tasks(const FMMTasks<dim>& tasks, const BlockVectorX& x) const
    {
        BlockVectorX out(R, VectorX(data.obs_locs.size(), 0.0));
        std::cout << tasks.m2ps.size() << std::endl;

        // label cells
        std::map<size_t,size_t> p2p_cell_map;
        size_t next_cell = 0;
        for (const auto& p2p: tasks.p2ps) {
            auto cell_key = p2p.obs_cell.data.indices[0];
            if (p2p_cell_map.count(cell_key) == 0) {
                p2p_cell_map[cell_key] = next_cell;
                next_cell++;
            }
        }

        std::map<size_t,size_t> m2p_cell_map;
        next_cell = 0;
        for (const auto& m2p: tasks.m2ps) {
            auto cell_key = m2p.obs_cell.data.indices[0];
            if (m2p_cell_map.count(cell_key) == 0) {
                m2p_cell_map[cell_key] = next_cell;
                next_cell++;
            }
        }

        //separate tasks by observer cell
        std::vector<std::vector<P2PTask<dim>>> p2p_sorted(p2p_cell_map.size());
        std::vector<std::vector<M2PTask<dim>>> m2p_sorted(m2p_cell_map.size());
        for (const auto& p2p: tasks.p2ps) {
            auto cell_key = p2p.obs_cell.data.indices[0];
            p2p_sorted[p2p_cell_map[cell_key]].push_back(p2p);
        }
        for (const auto& m2p: tasks.m2ps) {
            auto cell_key = m2p.obs_cell.data.indices[0];
            m2p_sorted[m2p_cell_map[cell_key]].push_back(m2p);
        }

#pragma omp parallel for 
        for (size_t i = 0; i < p2p_sorted.size(); i++) {
            for (const auto& t: p2p_sorted[i]) {
                P2P(out, t.obs_cell, t.src_cell, x);
            }
        }

#pragma omp parallel for
        for (size_t i = 0; i < m2p_sorted.size(); i++) {
            for (const auto& t: m2p_sorted[i]) {
                M2P(out, t.obs_cell, t.src_cell, t.p2m);
            }
        }

        return out;
    }

    BlockVectorX apply(const BlockVectorX& x) const 
    {
        assert(x.size() == C);
        int original_nested = omp_get_nested();
        omp_set_nested(1);

        CheckToEquiv check_to_equiv;
        build_check_to_equiv(src_oct, check_to_equiv);

        const auto p2m = P2M(src_oct, x, check_to_equiv);


        FMMTasks<dim> tasks;
        dual_tree(obs_oct, src_oct, *p2m, tasks);
        auto out = execute_tasks(tasks, x);

        omp_set_nested(original_nested);
        return out;
    }

    const Kernel<dim,R,C>& K;
    const NBodyData<dim> data;
    const TranslationSurface<dim> surface;
    Octree<dim> src_oct;
    Octree<dim> obs_oct;
    const FMMConfig config;

    FMMOperator(const Kernel<dim,R,C>& K, const NBodyData<dim>& data,
        const FMMConfig& config);

    void build_check_to_equiv(const Octree<dim>& cell, 
        CheckToEquiv& created_ops) const;
    
    std::vector<double> apply_src_to_check(const Octree<dim>& cell,
        const BlockVectorX& x) const;

    std::vector<double> apply_children_to_check(const Octree<dim>& cell,
        const typename P2MData<dim>::ChildrenType& child_p2m) const;

    std::unique_ptr<P2MData<dim>> P2M(const Octree<dim>& cell, const BlockVectorX& x,
        CheckToEquiv& check_to_equiv_ops) const;

    void P2P(BlockVectorX& out, const Octree<dim>& obs_cell,
        const Octree<dim>& src_cell, const BlockVectorX& x) const;

    void M2P(BlockVectorX& out, const Octree<dim>& obs_cell,
        const Octree<dim>& src_cell, const P2MData<dim>& p2m) const;
};

} // END namespace tbem
#endif

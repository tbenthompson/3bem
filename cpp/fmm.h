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
        auto r_ref = 0.5; //std::sqrt(2) + d;
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

    void dual_tree(const Octree<dim>& obs_cell,
        const Octree<dim>& src_cell, const P2MData<dim>& p2m,
        FMMTasks<dim>& tasks) const;

    BlockVectorX execute_tasks(const FMMTasks<dim>& tasks, const BlockVectorX& x) const;

    BlockVectorX apply(const BlockVectorX& x) const;
};

} // END namespace tbem
#endif

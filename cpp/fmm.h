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
        auto r_ref = 0.3; //std::sqrt(2) + d;
        return move(box, r_ref);
    }

    std::vector<Vec<double,dim>> 
    downward_check_points(const Box<dim>& box, double d) const
    {
        auto r_ref = std::sqrt(2) + d;
        return move(box, r_ref);
    }

    std::vector<Vec<double,dim>> 
    downward_equiv_points(const Box<dim>& box, double d) const
    {
        auto r_ref = 4.0 - std::sqrt(2) - 2 * d;
        return move(box, r_ref);
    }
};
template <size_t dim>
TranslationSurface<dim> make_surrounding_surface(size_t order);

typedef std::vector<SVDPtr> CheckToEquiv;

struct FMMConfig {
    const double mac;
    // const double inner_radius;
    // const double outer_radius;
    const size_t order;
    const size_t min_pts_per_cell;
    const double d;
    // const double alpha;
};

template <size_t dim>
struct FMMTasks 
{
    struct CellPairTask
    {
        const Octree<dim>& obs_cell;
        const Octree<dim>& src_cell;
    };

    struct CellTask
    {
        const Octree<dim>& cell;
    };

    std::vector<CellTask> p2ms;
    std::vector<CellTask> m2ms;
    std::vector<CellPairTask> m2ps;
    std::vector<CellPairTask> p2ps;
    std::vector<CellPairTask> p2ls;
    std::vector<CellPairTask> m2ls;
    std::vector<CellTask> l2ls;
    std::vector<CellTask> l2ps;
};

/* An implementation of the kernel independent fast multipole method
 * as described in:
 *
 * A Kernel-independent Adaptive Fast Multipole Method in Two and 
 * Three Dimensions. L. Ying, G. Biros, D. Zorin. Journal of Computational
 * Physics, 196(2), 591-626, 2004.
 */
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

    /* Construct the operator relating the influence of a set of
     * a sources on the check surface to the equivalent set of sources on
     * the equivalent surface.
     */
    void build_check_to_equiv(const Octree<dim>& cell, 
        CheckToEquiv& up_ops, CheckToEquiv& down_ops) const;
    
    /* The P2M operator converts sources to equivalent sources in the
     * upward tree traversal.
     */
    void P2M(const Octree<dim>& cell, const SVDPtr& check_to_equiv,
        const BlockVectorX& x, double* parent_multipoles) const;

    /* The M2M operator converts child cell equivalent sources to parent cell
     * equivalent sources in the upward tree traversal.
     */
    void M2M(const Octree<dim>& cell, const SVDPtr& check_to_equiv,
        std::vector<double*>& child_multipoles, double* parent_multipoles) const;

    /* The P2L operator converts from sources straight to the equivalent sources
     * at a farfield observation cell (local coefficients) bypassing the P2M, M2M 
     * and M2L operators in cases where there are insufficient sources to justify
     * using those operators.
     */
    void P2L(const Octree<dim>& obs_cell, const Octree<dim>& src_cell,
        const SVDPtr& check_to_equiv, const BlockVectorX& x,
        double* locals) const;

    /* The M2P operator evaluates the influence of source cell equivalent
     * sources (multipole coefficients) on farfield observation points bypassing
     * the M2L, L2L, L2P operators in cases where there are insufficient observation
     * points to justify using those operators.
     */
    void M2P(const Octree<dim>& obs_cell, const Octree<dim>& src_cell,
        double* multipoles, BlockVectorX& out) const;

    /* The M2L operator converts equivalent sources from an source tree 
     * cell (multipole coefficients) into equivalent sources for an
     * observation tree cell (local coefficients).
     */
    void M2L(const Octree<dim>& obs_cell, const Octree<dim>& src_cell,
        const SVDPtr& down_check_to_equiv, double* multipoles, 
        double* locals) const;

    /* Convert from equivalent sources at a farfield parent observation cell
     * to its child cells.
     */
    void L2L(const Octree<dim>& cell, const SVDPtr& check_to_equiv,
        double* parent_locals, std::vector<double*>& child_locals) const;

    /* The L2P operator evaluates the influence of the equivalent sources from
     * a observation cell on the observation points in that cell
     */
    void L2P(const Octree<dim>& cell, double* locals, BlockVectorX& out) const;

    /* Perform a direct n body calculation between a source and observation
     * cell.
     */
    void P2P(const Octree<dim>& obs_cell, const Octree<dim>& src_cell,
        const BlockVectorX& x, BlockVectorX& out) const;

    /* Determine the sequence of operations necessary to compute the upward
     * equivalent sources
     */
    void upward_traversal(const Octree<dim>& cell, FMMTasks<dim>& tasks) const;

    /* Determine the sequence of FMM operations for the interaction of a set of
     * sources and observers
     */
    void dual_tree(const Octree<dim>& obs_cell, const Octree<dim>& src_cell,
        FMMTasks<dim>& tasks) const;

    /* Determine the sequence of FMM operations for calculating the influence of
     * observation cell equivalent sources on observation points
     */
    void downward_traversal(const Octree<dim>& obs_cell, FMMTasks<dim>& tasks) const;

    BlockVectorX execute_tasks(const FMMTasks<dim>& tasks,
        const BlockVectorX& x, 
        const CheckToEquiv& up_check_to_equiv,
        const CheckToEquiv& down_check_to_equiv) const;

    BlockVectorX apply(const BlockVectorX& x) const;
};

} // END namespace tbem
#endif

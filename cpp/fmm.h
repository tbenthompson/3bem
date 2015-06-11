#ifndef TBEMKASDJLKJKJ_FMM_H
#define TBEMKASDJLKJKJ_FMM_H

#include <cassert>
#include <memory>
#include "kernel.h"
#include "octree.h"
#include "numbers.h"
#include "blas_wrapper.h"
#include "operator.h"

namespace tbem {

template <size_t dim>
struct TranslationSurface {
    const std::vector<Vec<double,dim>> pts;
    const std::vector<Vec<double,dim>> normals;

    std::vector<Vec<double,dim>> move(const Box<dim>& box) const
    {
        auto cell_r = hypot(box.half_width);
        auto out_r = cell_r;
        auto out_center = box.center;

        std::vector<Vec<double,dim>> out_pts(pts.size());
        for (size_t i = 0; i < pts.size(); i++) {
            out_pts[i] = out_r * pts[i] + out_center;
        }

        return out_pts;
    }

    TranslationSurface<dim> scale(double factor) const
    {
        std::vector<Vec<double,dim>> out_pts(pts.size());
        for (size_t i = 0; i < pts.size(); i++) {
            out_pts[i] = factor * pts[i];
        }
        return {out_pts, normals};
    }

    static TranslationSurface<dim> up_check_surface(size_t order, double d) 
    {
        auto r_ref = 4.0 - std::sqrt(2) - 2 * d;
        return make_surrounding_surface(order).scale(r_ref);
    }

    static TranslationSurface<dim> up_equiv_surface(size_t order, double d)
    {
        (void)d;
        auto r_ref = 0.3;//std::sqrt(2) + d;
        return make_surrounding_surface(order).scale(r_ref);
    }

    static TranslationSurface<dim> down_check_surface(size_t order, double d)
    {
        (void)d;
        auto r_ref = 0.3;//std::sqrt(2) + d;
        return make_surrounding_surface(order).scale(r_ref);
    }

    static TranslationSurface<dim> down_equiv_surface(size_t order, double d)
    {
        auto r_ref = 4.0 - std::sqrt(2) - 2 * d;
        return make_surrounding_surface(order).scale(r_ref);
    }

    static TranslationSurface<dim> make_surrounding_surface(size_t order);
};

typedef std::vector<std::vector<double>> CheckToEquiv;

struct FMMConfig {
    const double mac;

    // Note that the order of approximation will only be approximately 
    // satisfied in 3D. For even distribution of point on a sphere in 3D,
    // the translation surface may have a few more or less points. As a 
    // result, config.order is not perfectly representative of the number of
    // multipole or local coefficients either and should not be used as a
    // loop limit.
    const size_t order;

    const size_t min_pts_per_cell;
    const double d;
    //TODO: Add the SVD threshold as a parameter
    //TODO: Add the translation surface inner and outer radii as parameters.
    //
    const bool account_for_small_cells;

    FMMConfig(double mac, size_t order, size_t min_pts_per_cell,
        double d, bool account_for_small_cells):
        mac(mac), order(order), min_pts_per_cell(min_pts_per_cell),
        d(d), account_for_small_cells(account_for_small_cells)
    {}
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

template <size_t dim> struct NBodyData;

/* An implementation of the kernel independent fast multipole method
 * as described in:
 *
 * A Kernel-independent Adaptive Fast Multipole Method in Two and 
 * Three Dimensions. L. Ying, G. Biros, D. Zorin. Journal of Computational
 * Physics, 196(2), 591-626, 2004.
 */
template <size_t dim, size_t R, size_t C>
struct FMMOperator: public OperatorI {
    const std::shared_ptr<Kernel<dim,R,C>> K;
    const NBodyData<dim> data;
    const TranslationSurface<dim> up_equiv_surface;
    const TranslationSurface<dim> up_check_surface;
    const TranslationSurface<dim> down_equiv_surface;
    const TranslationSurface<dim> down_check_surface;
    Octree<dim> src_oct;
    Octree<dim> obs_oct;
    const FMMConfig config;
    const CheckToEquiv up_check_to_equiv;
    const CheckToEquiv down_check_to_equiv;
    FMMTasks<dim> tasks;

    FMMOperator(const Kernel<dim,R,C>& K, const NBodyData<dim>& data,
        const FMMConfig& config);

    /* Construct the operators relating the influence of a set of
     * a sources on the check surface to the equivalent set of sources on
     * the equivalent surface.
     */
    CheckToEquiv build_check_to_equiv(const Octree<dim>& cell,
        size_t start_level, const TranslationSurface<dim>& equiv_surf,
        const TranslationSurface<dim>& check_surf) const;
    
    /* The P2M operator converts sources to equivalent sources in the
     * upward tree traversal.
     */
    void P2M(const Octree<dim>& cell, const std::vector<double>& check_to_equiv,
        const std::vector<double>& x, double* parent_multipoles) const;

    /* The M2M operator converts child cell equivalent sources to parent cell
     * equivalent sources in the upward tree traversal.
     */
    void M2M(const Octree<dim>& cell, const std::vector<double>& check_to_equiv,
        std::vector<double*>& child_multipoles, double* parent_multipoles) const;

    /* The P2L operator converts from sources straight to the equivalent sources
     * at a farfield observation cell (local coefficients) bypassing the P2M, M2M 
     * and M2L operators in cases where there are insufficient sources to justify
     * using those operators.
     */
    void P2L(const Octree<dim>& obs_cell, const Octree<dim>& src_cell,
        const std::vector<double>& check_to_equiv, const std::vector<double>& x,
        double* locals) const;

    /* The M2P operator evaluates the influence of source cell equivalent
     * sources (multipole coefficients) on farfield observation points bypassing
     * the M2L, L2L, L2P operators in cases where there are insufficient observation
     * points to justify using those operators.
     */
    void M2P(const Octree<dim>& obs_cell, const Octree<dim>& src_cell,
        double* multipoles, std::vector<double>& out) const;

    /* The M2L operator converts equivalent sources from an source tree 
     * cell (multipole coefficients) into equivalent sources for an
     * observation tree cell (local coefficients).
     */
    void M2L(const Octree<dim>& obs_cell, const Octree<dim>& src_cell,
        const std::vector<double>& down_check_to_equiv, double* multipoles, 
        double* locals) const;

    /* Convert from equivalent sources at a farfield parent observation cell
     * to its child cells.
     */
    void L2L(const Octree<dim>& cell, const std::vector<double>& check_to_equiv,
        double* parent_locals, std::vector<double*>& child_locals) const;

    /* The L2P operator evaluates the influence of the equivalent sources from
     * a observation cell on the observation points in that cell
     */
    void L2P(const Octree<dim>& cell, double* locals, std::vector<double>& out) const;

    /* Perform a direct n body calculation between a source and observation
     * cell.
     */
    void P2P(const Octree<dim>& obs_cell, const Octree<dim>& src_cell,
        const std::vector<double>& x, std::vector<double>& out) const;

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

    std::vector<double> execute_tasks(const FMMTasks<dim>& tasks,
        const std::vector<double>& x, 
        const CheckToEquiv& up_check_to_equiv,
        const CheckToEquiv& down_check_to_equiv) const;

    std::vector<double> apply(const std::vector<double>& x) const;

    virtual size_t n_rows() const {return data.obs_locs.size() * R;}
    virtual size_t n_cols() const {return data.src_locs.size() * C;}
};

} // END namespace tbem
#endif

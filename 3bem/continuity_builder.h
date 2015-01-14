#ifndef __ASKDJAWERWFJS_CONTINUITY_BUILDER_H

#define __ASKDJAWERWFJS_CONTINUITY_BUILDER_H
#include <unordered_map>

#include "constraint.h"
#include "vertex_iterator.h"
#include "vec_ops.h"
    
namespace tbem {

template <size_t dim>
bool pts_close(Vec<double,dim> a, Vec<double,dim> b, double eps) {
    return all(fabs(a - b) < eps);
}

template <size_t dim>
bool pts_not_close(Vec<double,dim> a, Vec<double,dim> b, double eps) {
    return !pts_close(a, b, eps);
}

template <size_t dim>
using OverlapPair = std::pair<VertexIterator<dim>,VertexIterator<dim>>;

template <size_t dim>
using OverlapMap = std::unordered_multimap<
    VertexIterator<dim>,
    VertexIterator<dim>,
    HashVertexIterator<dim>
>;

template <size_t dim>
OverlapMap<dim>
find_overlapping_vertices(const VertexIterator<dim>& A_begin,
    const VertexIterator<dim>& B_begin, double eps = 1e-15) 
{
    std::vector<OverlapPair<dim>> overlaps;
    for (auto A_it = A_begin; !A_it.is_end(); ++A_it) {
        auto A_pt = *A_it;
        for (auto B_it = B_begin; !B_it.is_end(); ++B_it) {
            auto B_pt = *B_it;
            if (pts_not_close(A_pt, B_pt, eps)) {
                continue;
            }
            overlaps.push_back({A_it, B_it});
        }
    }
    return OverlapMap<dim>(overlaps.begin(), overlaps.end());
}

template <size_t dim>
OverlapMap<dim>
find_overlapping_vertices_same_mesh(const VertexIterator<dim>& A_begin, 
    double eps = 1e-15) 
{
    auto overlap_map = find_overlapping_vertices(A_begin, A_begin);
    std::vector<OverlapPair<dim>> filtered_overlaps;
    for (const auto& o: overlap_map) {
        if (o.first.absolute_index() < o.second.absolute_index()) {
            filtered_overlaps.push_back({o.first, o.second});
        }
    }
    return OverlapMap<dim>(filtered_overlaps.begin(), filtered_overlaps.end());
}

template <size_t dim>
bool continuity_pair_crosses_cut(const VertexIterator<dim>& facet0_it,
    const VertexIterator<dim>& facet1_it, const VertexIterator<dim>& discontinuity_it) 
{
    const auto& facet0 = facet0_it.get_facet();
    const auto& facet1 = facet1_it.get_facet();
    const auto& disc_facet = discontinuity_it.get_facet();

    auto facet0_side = which_side_facet(disc_facet.vertices, facet0.vertices);
    auto facet1_side = which_side_facet(disc_facet.vertices, facet1.vertices);

    if (facet0_side == facet1_side) {
        return false;
    }
    return true;
}

template <size_t dim>
OverlapMap<dim> mesh_continuity(const VertexIterator<dim>& begin_iter) 
{
    return find_overlapping_vertices_same_mesh(begin_iter);
}

template <size_t dim>
OverlapMap<dim> cut_at_intersection(const OverlapMap<dim>& continuity,
    const VertexIterator<dim>& mesh_begin_iter,
    const VertexIterator<dim>& cut_begin_iter) 
{
    auto out_continuity = continuity;

    auto cut_pairs = find_overlapping_vertices(mesh_begin_iter, cut_begin_iter);
    for (const auto& p: cut_pairs) {
        const auto& continuity_it = p.first; 
        const auto& discontinuity_it = p.second;

        auto range = out_continuity.equal_range(continuity_it);
        for (auto potential_cut = range.first;
            potential_cut != range.second;) 
        {
            if (continuity_pair_crosses_cut(potential_cut->first,
                potential_cut->second, discontinuity_it)) 
            {
                out_continuity.erase(potential_cut++);
            } else {
                ++potential_cut;
            }
        }
    }

    return out_continuity;
}

template <size_t dim>
std::vector<ConstraintEQ> convert_to_constraints(const OverlapMap<dim>& continuity) {
    std::vector<ConstraintEQ> constraints;
    for (const auto& p: continuity) {
        auto c = continuity_constraint(p.first.absolute_index(),
            p.second.absolute_index());
        constraints.push_back(c);
    }
    return constraints;
}

template <size_t dim> 
std::vector<ConstraintEQ> form_neighbor_bcs(
    const VertexIterator<dim>& continuous_mesh,
    const VertexIterator<dim>& neighbor_mesh,
    const std::vector<double>& values
)
{
    auto overlaps = find_overlapping_vertices(continuous_mesh, neighbor_mesh);
    std::vector<ConstraintEQ> constraints;
    for (const auto& o: overlaps) {
        auto constrained_dof = o.first.absolute_index();
        auto bc_val = values[o.second.absolute_index()];
        auto c = boundary_condition(constrained_dof, bc_val);
        constraints.push_back(c);    
    }
    return constraints;
}

} //END namespace tbem

#endif

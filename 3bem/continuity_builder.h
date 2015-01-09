#ifndef __ASKDJAWERWFJS_CONTINUITY_BUILDER
#define __ASKDJAWERWFJS_CONTINUITY_BUILDER
#include <unordered_map>

#include "constraint.h"
#include "mesh.h"
#include "vertex_iterator.h"
    
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
        if (o.first < o.second) {
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
struct ContinuityBuilder {
    const VertexIterator<dim> mesh_iter;
    OverlapMap<dim> continuity_map;

    ContinuityBuilder(const VertexIterator<dim>& begin_iter):
        mesh_iter(begin_iter)
    {
        continuity_map = find_overlapping_vertices_same_mesh(begin_iter);
    }

    void apply_discontinuities(const VertexIterator<dim>& cut_begin_iter) {
        auto cut_pairs = find_overlapping_vertices(mesh_iter, cut_begin_iter);
        for (const auto& p: cut_pairs) {
            const auto& continuity_it = p.first; 
            const auto& discontinuity_it = p.second;

            auto range = continuity_map.equal_range(continuity_it);
            for (auto potential_cut = range.first; potential_cut != range.second;) 
            {
                if (continuity_pair_crosses_cut(potential_cut->first,
                    potential_cut->second, discontinuity_it)) 
                {
                    continuity_map.erase(potential_cut++);
                } else {
                    ++potential_cut;
                }
            }
        }
    }

    std::vector<ConstraintEQ> build() {
        std::vector<ConstraintEQ> constraints;
        for (const auto& p: continuity_map) {
            auto c = continuity_constraint(p.first.absolute_index(),
                p.second.absolute_index());
            constraints.push_back(c);
        }
        return constraints;
    }
};

} //END namespace tbem

#endif

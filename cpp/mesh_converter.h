#ifndef TBEM_MESH_FORMAT_CONVERTER_H
#define TBEM_MESH_FORMAT_CONVERTER_H

#include "vec.h"
#include "vertex_iterator.h"
#include "continuity_builder.h"
#include <vector>
#include <map>

namespace tbem {

template <size_t dim>
struct PtIndexMesh {
    std::vector<Vec<double,dim>> points;
    std::vector<Vec<size_t,dim>> facets;
};

template <size_t dim>
PtIndexMesh<dim> convert_facet_to_pt_index(
    const std::vector<Vec<Vec<double,dim>,dim>>& facets)
{
    Mesh<dim> m{facets};
    auto overlaps = find_overlapping_vertices_same_mesh(m.begin());
    auto n_input_vertices = m.n_facets() * dim;

    std::map<size_t,size_t> in_out_vertex_mapping;
    std::vector<Vec<double,dim>> out_points;
    std::vector<Vec<size_t,dim>> out_pt_idxs(m.n_facets());
    size_t next_out_idx = 0;
    for (auto v = m.begin(); v != m.end(); ++v)
    {
        if (in_out_vertex_mapping.count(v.absolute_index()) > 0) 
        {
            size_t out_idx = in_out_vertex_mapping.find(v.absolute_index())->second;
            out_pt_idxs[v.facet_idx][v.vertex_idx] = out_idx;
            continue;
        }
        in_out_vertex_mapping.insert({v.absolute_index(), next_out_idx});
        out_points.push_back(v.get_vertex());
        out_pt_idxs[v.facet_idx][v.vertex_idx] = next_out_idx;
        auto this_v_overlaps = overlaps.equal_range(v);
        for (auto o = this_v_overlaps.first; o != this_v_overlaps.second; ++o)
        {
            in_out_vertex_mapping.insert({o->second.absolute_index(), next_out_idx});
        }
        next_out_idx++;
    }

    return {out_points, out_pt_idxs};
}

template <size_t NP, size_t dim>
std::vector<Vec<Vec<double,dim>,NP>> convert_pt_index_to_facet(
    const std::vector<Vec<size_t,NP>>& pt_idxs,
    const std::vector<Vec<double,dim>>& pts)
{
    std::vector<Vec<Vec<double,dim>,NP>> facets(pt_idxs.size());
    for (size_t i = 0; i < pt_idxs.size(); i++) {
        for (size_t v = 0; v < NP; v++) {
            facets[i][v] = pts[pt_idxs[i][v]];
        }
    }
    return facets;
}

} //end namespace tbem

#endif

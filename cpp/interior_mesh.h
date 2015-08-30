#ifndef TBEM_INTERIOR_MESH_H
#define TBEM_INTERIOR_MESH_H

#include <unordered_map>
#include <unordered_set>

namespace tbem {

typedef std::unordered_map<size_t,std::unordered_set<size_t>> PtTriMap;

//I think the data structures here are suboptimal.

template <size_t dim>
PtTriMap build_pt_to_tri_map(const std::vector<Vec<size_t,dim>>& tris)
{
    PtTriMap pt_to_tri_map;
    for (size_t i = 0; i < tris.size(); i++) {
        for (size_t d = 0; d < dim; d++) {
            pt_to_tri_map[tris[i][d]].insert(i);
        }
    }
    return pt_to_tri_map;
}

bool are_adjacent_tris(const Vec<size_t,3> tri_A, const Vec<size_t,3> tri_B)
{
    size_t touching_verts = 0;
    for (size_t d1 = 0; d1 < 3; d1++) {
        for (size_t d2 = 0; d2 < 3; d2++) {
            if (tri_A[d1] == tri_B[d2]) {
                touching_verts += 1;
                continue;
            }
        }
    }
    if (touching_verts == 2) {
        return true;
    }
    return false;
}

std::unordered_set<size_t> find_adjacent_tris(const std::vector<Vec<size_t,3>>& tris,
    const Vec<size_t,3>& query_tri, const PtTriMap& pt_to_tri_map)
{
    std::unordered_set<size_t> adjacent_tris;
    for (size_t d = 0; d < 3; d++) {
        for (size_t touching_idx: pt_to_tri_map.at(query_tri[d])) {
            if (are_adjacent_tris(tris[touching_idx], query_tri)) {
                adjacent_tris.insert(touching_idx);
            }
        }
    }
    return adjacent_tris;
}

bool across_bdry_edge(const Vec<size_t,3>& tri_A, const Vec<size_t,3>& tri_B,
    const PtTriMap& pt_bdry_facet_map)
{
    std::unordered_map<size_t,size_t> facet_touch_counts;
    for (size_t d = 0; d < 3; d++) {
        if (pt_bdry_facet_map.count(tri_A[d]) > 0) {
            for (auto bdry_facet: pt_bdry_facet_map.at(tri_A[d])) {
                facet_touch_counts[bdry_facet] += 1;
            }
        }
        if (pt_bdry_facet_map.count(tri_B[d]) > 0) {
            for (auto bdry_facet: pt_bdry_facet_map.at(tri_B[d])) {
                facet_touch_counts[bdry_facet] += 1;
            }
        }
    }
    for (auto it = facet_touch_counts.begin(); it != facet_touch_counts.end(); ++it) {
        if (it->second == 4) {
            return true;
        }
    }
    return false;
}

//TODO: This is written with only 2D in mind.
std::vector<size_t> identify_regions(const std::vector<Vec<size_t,3>>& tris,
    const std::vector<Vec<size_t,2>>& bdry_facets)
{
    auto pt_to_tri_map = build_pt_to_tri_map(tris);
    auto pt_to_facet_map = build_pt_to_tri_map(bdry_facets);

    auto all_tri_idxs = range<size_t>(tris.size());
    std::unordered_set<size_t> tris_remaining(
        all_tri_idxs.begin(), all_tri_idxs.end()
    );

    std::vector<size_t> regions(tris.size());

    size_t next_label = 0;
    while (!tris_remaining.empty()) {
        std::set<size_t> processed_tris;
        std::vector<size_t> remaining_connected_tris;
        remaining_connected_tris.push_back(*tris_remaining.begin());

        while (!remaining_connected_tris.empty()) {
            auto cur_tri = remaining_connected_tris.back();
            remaining_connected_tris.pop_back();
            processed_tris.insert(cur_tri);
            regions[cur_tri] = next_label;
            tris_remaining.erase(cur_tri);

            auto adj_tris = find_adjacent_tris(tris, tris[cur_tri], pt_to_tri_map); 
            for (auto adj_idx: adj_tris) {
                if (processed_tris.count(adj_idx) > 0) {
                    continue;
                }
                if (!across_bdry_edge(tris[adj_idx], tris[cur_tri], pt_to_facet_map)) {
                    remaining_connected_tris.push_back(adj_idx);
                }
            }
        }


        next_label++;
    }

    return regions;
}

} //end namespace tbem
#endif

#include "containing_tri.h"
#include "mesh_converter.h"
#include "geometry.h"
#include "gte_wrapper.h"

namespace tbem {

std::pair<bool,size_t> find_containing_tri_idx(const Vec<double,2>& query_pt,
    const std::vector<Vec<size_t,3>>& tris,
    const std::vector<Vec<double,2>>& pts)
{
    auto tri_verts = convert_pt_index_to_facet(tris, pts);
    auto balls = make_facet_balls(tri_verts);
    auto intersections = intersect_balls_all_pairs({{query_pt, 0.0}}, balls);
    for (size_t i = 0; i < intersections.size(); i++) {
        auto tri_idx = intersections[i].second;
        if(is_point_in_triangle(query_pt, tri_verts[tri_idx])) {
            return {true, tri_idx};
        }
    }
    return {false, 0};
}

} //end namespace tbem

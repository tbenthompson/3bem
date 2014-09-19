#include "bem.h"
#include "numerics.h"
#include <iostream>
#include <algorithm>
#include <unordered_map>

void refine_edge(Mesh& new_mesh, std::array<int, 2> seg) {
    // Find the new vertex and segments.
    const auto v0 = new_mesh.vertices[seg[0]];
    const auto v1 = new_mesh.vertices[seg[1]];
    const std::array<double,2> midpt = {(v0[0] + v1[0]) / 2, (v0[1] + v1[1]) / 2};

    int next_vert = new_mesh.vertices.size();
    const std::array<int,2> new_seg0 = {seg[0], next_vert};
    const std::array<int,2> new_seg1 = {next_vert, seg[1]};

    new_mesh.vertices.push_back(midpt);
    new_mesh.segments.push_back(new_seg0);
    new_mesh.segments.push_back(new_seg1);
}

Mesh refine_mesh(const Mesh& m, std::vector<int> refine_these) {
    Mesh new_mesh;
    //Keep all the vertices (refinement should never remove a vertex).
    new_mesh.vertices = m.vertices;
    
    // Sort the refined edges so that we only have to check the
    // next one at any point in the loop.
    std::sort(refine_these.begin(), refine_these.end());

    // The next refined edge.
    int current = 0;

    for (int i = 0; i < (int)m.segments.size(); i++) {
        if (i == refine_these[current]) {
            current += 1;
            refine_edge(new_mesh, m.segments[i]);
        } else {
            new_mesh.segments.push_back(m.segments[i]);
        }
    }

    return new_mesh;
}

Subsegments build_ref_subsegs(QuadratureRule& quad) {
    Subsegments ref;
    ref.ref_left.resize(quad.size());
    ref.ref_center.resize(quad.size());
    ref.ref_right.resize(quad.size());

    ref.ref_left[0] = -1.0;
    ref.ref_center[0] = quad[0].first;
    for (unsigned int i = 1; i < quad.size(); i++) {
        ref.ref_center[i] = quad[i].first;
        ref.ref_left[i] = (quad[i - 1].first + quad[i].first) / 2.0;
        ref.ref_right[i - 1] = ref.ref_left[i];
        /* std::cout << ref_subsegs.left[i] << std::endl; */
    }
    ref.ref_right[quad.size() - 1] = 1.0;

    return ref;
}

Subsegments get_src_obs(Mesh& m, QuadratureRule& quad_rule)
{
    auto subsegs = build_ref_subsegs(quad_rule);
    for (unsigned int i = 0; i < m.segments.size(); i++) {
        const auto v0 = m.vertices[m.segments[i][0]];
        const auto v1 = m.vertices[m.segments[i][1]];
        for (unsigned int k = 0; k < quad_rule.size(); k++) {
            const double x = ref_to_real(quad_rule[k].first, v0[0], v1[0]);
            const double y = ref_to_real(quad_rule[k].first, v0[1], v1[1]);
            subsegs.center.push_back({x, y});
            subsegs.owner.push_back(i);
        }
    }
    return subsegs;
}

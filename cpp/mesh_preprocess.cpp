#include "mesh_preprocess.h"
#include "geometry.h"
#include "intersect_balls.h"
#include "gte_wrapper.h"
#include <algorithm>

namespace tbem {

double intersection_eps = 1e-14;


template <size_t dim>
std::vector<FacetIntersection<dim>> MeshPreprocessor<dim>::find_intersections(
    const std::vector<Vec<Vec<double,dim>,dim>>& facetsA,
    const std::vector<Vec<Vec<double,dim>,dim>>& facetsB)
{
    auto ballsA = make_facet_balls(facetsA);
    auto ballsB = make_facet_balls(facetsB);
    auto ball_intersections = intersect_balls_all_pairs(ballsA, ballsB);
    std::vector<FacetIntersection<dim>> out;
    for (size_t i = 0; i < ball_intersections.size(); i++) {
        // Due to finite precision arithmetic, some endpoint intersections will
        // not be found. To get around this, expand each facet by a very small
        // amount (order machine epsilon) so that the intersection will
        // definitely be within the facet. Even if it catches some intersections
        // that don't actually exist, this is good, because it's okay to fail with
        // ridiculous geometries that vary of scales of order machine epsilon.
        auto A_idx = ball_intersections[i].first;
        auto B_idx = ball_intersections[i].second;
        auto expandedA = expand_facet(facetsA[A_idx], intersection_eps);
        auto expandedB = expand_facet(facetsB[B_idx], intersection_eps);
        auto real_intersection = facet_facet_intersection(expandedA, expandedB);
        if (real_intersection.size() == 0) {
            continue;
        }
        out.push_back({
            ball_intersections[i].first,
            ball_intersections[i].second,
            real_intersection
        });
    }
    return out;
}

template <size_t dim>
bool intersection_is_endpoint(const Vec<Vec<double,dim>,dim>& f,
    const FacetIntersection<dim>& intersection) 
{
    double threshold = intersection_eps * 10 * facet_ball(f).radius;
    for (size_t d = 0; d < dim; d++) {
        if (hypot(f[d] - intersection.pts[0]) < threshold) {
            return true;
        }
    }
    return false;
}

template <size_t dim>
std::vector<Vec<Vec<double,dim>,dim>> split_facet(const Vec<Vec<double,dim>,dim>& f,
    const FacetIntersection<dim>& intersection);

template <>
std::vector<Vec<Vec<double,2>,2>> split_facet<2>(const Vec<Vec<double,2>,2>& f,
    const FacetIntersection<2>& intersection) 
{
    return {
        {f[0], intersection.pts[0]},
        {intersection.pts[0], f[1]}
    };
}

template <>
std::vector<Vec<Vec<double,3>,3>> split_facet<3>(const Vec<Vec<double,3>,3>& f,
    const FacetIntersection<3>& intersection) 
{
    (void)f;
    (void)intersection;
    //TODO: Not implemented yet.
    return {};
}

template <size_t dim>
std::vector<Vec<Vec<double,dim>,dim>> 
MeshPreprocessor<dim>::split_facets_at_intersections(
    const std::vector<Vec<Vec<double,dim>,dim>>& facetsA,
    const std::vector<FacetIntersection<dim>>& intersections)
{
    //TODO: splitting segments is easy, splitting triangles will be harder. i
    //think it should be done in two steps. First, split the triangle into a
    //a new triangle and a quadrilateral. The split the quadrilateral into two
    //triangles. Then, use some other mesh preprocessing step to split/fix ugly
    //triangles.

    std::vector<Vec<Vec<double,dim>,dim>> out_facets;

    auto sorted_intersections = intersections;
    std::sort(sorted_intersections.begin(), sorted_intersections.end(),
        [](const FacetIntersection<dim>& i0, const FacetIntersection<dim>& i1) {
            return i0.facet_idx_A < i1.facet_idx_A;
        });

    size_t current = 0;

    // TODO: Very similar to the Mesh::refine function in mesh.cpp. Can I consolidate
    // the two functions somehow in order to prevent bugs?
    // TODO: BAD
    // TODO: BAD
    // TODO: BAD
    // TODO: BAD
    for (size_t i = 0; i < facetsA.size(); i++) {
        if (current == sorted_intersections.size()) {
            out_facets.push_back(facetsA[i]);
            continue;
        }
        bool facet_is_intersected = i == sorted_intersections[current].facet_idx_A;
        if (facet_is_intersected) {
            if (!intersection_is_endpoint(facetsA[i], sorted_intersections[current])) {
                auto split = split_facet(facetsA[i], sorted_intersections[current]);
                for (auto s: split) {
                    out_facets.push_back(s);
                }
            } else {
                out_facets.push_back(facetsA[i]);
            }
            current++;
        } else {
            out_facets.push_back(facetsA[i]);
        }
    }
    return out_facets;
}

template struct MeshPreprocessor<2>;
template struct MeshPreprocessor<3>;
} // end namespace tbem

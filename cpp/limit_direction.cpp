#include "limit_direction.h"
#include "geometry.h"
#include "nearest_neighbors.h"
#include "gte_wrapper.h"
#include "intersect_balls.h"

namespace tbem {

template <size_t dim>
std::vector<Ball<dim>> make_facet_balls(const std::vector<Vec<Vec<double,dim>,dim>>& f)
{
    std::vector<Ball<dim>> out(f.size());
    for (size_t i = 0; i < f.size(); i++) {
        out[i] = facet_ball(f[i]);
    }
    return out;
}

template <size_t dim>
NearfieldFacetFinder<dim>::NearfieldFacetFinder(
    const std::vector<Vec<Vec<double,dim>,dim>>& facets):
    facets(facets),
    facet_balls(make_facet_balls(facets)),
    oct(make_octree<dim>(facet_balls, n_facets_per_leaf))
{}

template <size_t dim> 
NearfieldFacets<dim> 
NearfieldFacetFinder<dim>::find(const Vec<double,dim>& pt)
{
    //steps:
    //-- find the nearest facet (nearest neighbors search)
    assert(facets.size() > 0);
    auto nearest_neighbor = nearest_facet_brute_force(pt, facets);
    auto closest_facet_idx = nearest_neighbor.idx;
    auto closest_pt = nearest_neighbor.pt;

    //-- determine the search radius for nearfield facets
    auto closest_facet = facets[closest_facet_idx];
    auto search_ball = facet_ball(closest_facet);
    if (search_ball.radius < nearest_neighbor.distance) {
        return {{}, closest_facet, closest_pt, nearest_neighbor.distance};
    }

    //-- find the facet balls intersecting the sphere with that radius
    auto near_ball_indices = intersect_balls(search_ball, facet_balls, oct);

    //-- filter the facet balls to find the facets that actually intersected 
    //by the sphere (rather than just the surrounding ball)
    std::vector<Vec<Vec<double,dim>,dim>> close_facets;
    for (size_t facet_idx: near_ball_indices) {
        if (facet_idx == closest_facet_idx) {
            continue;
        }
        auto ref_pt = closest_pt_facet(pt, facets[facet_idx]);
        auto mesh_pt = ref_to_real(ref_pt, facets[facet_idx]);
        auto dist_to_mesh = dist(pt, mesh_pt);
        if (dist_to_mesh <= search_ball.radius) {
            close_facets.push_back(facets[facet_idx]);
        } 
    }

    return {close_facets, closest_facet, closest_pt, nearest_neighbor.distance};
}

template struct NearfieldFacetFinder<2>;
template struct NearfieldFacetFinder<3>;

template <size_t dim>
Vec<double,dim> decide_limit_dir(const Vec<double,dim>& pt,
    const NearfieldFacets<dim>& nearfield_facets, double safety_factor, 
    double epsilon) 
{
    const auto& closest_facet = nearfield_facets.nearest_facet;
    // TODO: The length scale here needs to be exactly double the length scale
    // used in finding the nearest facets, which means that the functions should
    // be unified. Basically, the function to find nearfield_facets deserves to be
    // in this file.
    double len_scale = facet_ball(closest_facet).radius * 2;

    // If the point isn't lying precisely on any facet, then no limit 
    // needs to be taken, in which case the limit direction is the 0 vector.
    double threshold = len_scale * epsilon;
    if (nearfield_facets.distance > threshold) {
        return zeros<Vec<double,dim>>::make();
    }

    // In the plane of the facet, the limit direction will point towards the 
    // centroid of the facet. Out of the plane of the facet, the limit direction
    // will point in the same direction as the facet's normal, with a magnitude
    // similar to the length scale of the facet.
    auto close_center = centroid(closest_facet);
    auto close_dir = unscaled_normal(closest_facet); 
    auto end_pt = close_center + normalized(close_dir) * len_scale;

    // Check if the vector intersects any facets.
    // If it does, back up to halfway between the intersection point and the
    // singular point
    for (auto f: nearfield_facets.facets) {
        assert(f != closest_facet);
        std::vector<Vec<double,dim>> intersections =
            seg_facet_intersection<dim>(f, {end_pt, pt});
        assert(intersections.size() != 2);
        if (intersections.size() == 1 && intersections[0] != pt) {
            end_pt = pt + (intersections[0] - pt) / 2.0;
        }
    }

    // Scaling from facet length scale down after doing the 
    // intersection tests is safer because a wider range of
    // facet intersections are checked and avoided.
    return (end_pt - pt) * safety_factor;
}

template Vec<double,2> 
decide_limit_dir(const Vec<double,2>&, const NearfieldFacets<2>&, double, double);

template Vec<double,3>
decide_limit_dir(const Vec<double,3>& pt, const NearfieldFacets<3>&, double, double);

} // end namespace tbem

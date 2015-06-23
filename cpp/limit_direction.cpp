#include "limit_direction.h"
#include "geometry.h"
#include "nearest_neighbors.h"
#include "gte_wrapper.h"
#include "intersect_balls.h"

namespace tbem {

template <size_t dim>
NearfieldFacetFinder<dim>::NearfieldFacetFinder(
    const std::vector<Vec<Vec<double,dim>,dim>>& facets,
    double far_threshold):
    nn_data(facets, 20),
    far_threshold(far_threshold)
{}



template <size_t dim> 
NearfieldFacets<dim> 
NearfieldFacetFinder<dim>::find(const Vec<double,dim>& pt)
{
    //steps:
    //-- find the nearest facet (nearest neighbors search)
    auto nearest_neighbor = nearest_facet(pt, nn_data);
    auto closest_facet_idx = nearest_neighbor.idx;

    if (nn_data.facets.size() == 0) {
        return {
            {}, 0, zeros<Vec<double,dim-1>>::make(),
            zeros<Vec<double,dim>>::make(), 0
        };
    }

    //-- determine the search radius for nearfield facets
    auto closest_facet = nn_data.facets[closest_facet_idx];
    auto closest_facet_ball = facet_ball(closest_facet);
    Ball<dim> search_ball{
        closest_facet_ball.center, closest_facet_ball.radius * far_threshold
    };

    //-- find the facet balls intersecting the sphere with that radius
    auto near_ball_indices = intersect_balls(
        search_ball, nn_data.facet_balls, nn_data.oct
    );

    //-- filter the facet balls to find the facets that actually intersected 
    //by the sphere (rather than just the surrounding ball)
    std::vector<size_t> close_facet_indices;
    for (size_t facet_idx: near_ball_indices) {
        auto ref_pt = closest_pt_facet(pt, nn_data.facets[facet_idx]);
        auto mesh_pt = ref_to_real(ref_pt, nn_data.facets[facet_idx]);
        auto dist_to_mesh = dist(pt, mesh_pt);
        if (dist_to_mesh <= search_ball.radius) {
            close_facet_indices.push_back(facet_idx);
        } 
    }

    return {
        close_facet_indices, closest_facet_idx, nearest_neighbor.ref_pt,
        nearest_neighbor.pt, nearest_neighbor.distance
    };
}

template struct NearfieldFacetFinder<2>;
template struct NearfieldFacetFinder<3>;

template <size_t dim>
Vec<double,dim> backup_halfway_from_intersection(Vec<double,dim> end_pt, 
    double len_scale, const Vec<double,dim>& pt,
    const NearfieldFacets<dim>& nearfield_facets, 
    const std::vector<Vec<Vec<double,dim>,dim>>& facets)
{
    for (auto f_idx: nearfield_facets.facet_indices) {
        if (f_idx == nearfield_facets.nearest_facet_idx) {
            continue;
        }
        auto& f = facets[f_idx];
        std::vector<Vec<double,dim>> intersections =
            seg_facet_intersection<dim>(f, {end_pt, pt});
        assert(intersections.size() != 2);
        if (intersections.size() == 1 && intersections[0] != pt) {
            auto to_intersection_dir = intersections[0] - pt;
            if (hypot(to_intersection_dir) < 1e-12 * len_scale) {
                continue;
            }
            end_pt = pt + to_intersection_dir / 2.0;
        }
    }
    return end_pt;
}

template <size_t dim>
Vec<double,dim> decide_limit_dir(const Vec<double,dim>& pt,
    const NearfieldFacets<dim>& nearfield_facets, 
    const std::vector<Vec<Vec<double,dim>,dim>>& facets, 
    double safety_factor, 
    double epsilon) 
{
    if (nearfield_facets.facet_indices.size() == 0) {
        return zeros<Vec<double,dim>>::make();
    }

    const auto& closest_facet = facets[nearfield_facets.nearest_facet_idx];

    // The length scale here needs to be exactly double the length scale
    // used in finding the nearest facets.
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

    // Make sure that the end point is within the interior of the mesh (the
    // line segment from the centroid of the element to the end_pt does not
    // intersect any edges)
    // TODO:
    auto interior_pt = end_pt; 
    // interior_pt = backup_halfway_from_intersection(
    //     end_pt, len_scale, close_center, nearfield_facets, facets
    // );

    // Check if the vector intersects any facets.
    // If it does, back up to halfway between the intersection point and the
    // singular point
    auto backup_pt = backup_halfway_from_intersection(
        interior_pt, len_scale, pt, nearfield_facets, facets
    );

    // Scaling from facet length scale down after doing the 
    // intersection tests is safer because a wider range of
    // facet intersections are checked and avoided.
    auto limit_dir = (backup_pt - pt) * safety_factor;
    return limit_dir;
}

template Vec<double,2> 
decide_limit_dir(const Vec<double,2>&, const NearfieldFacets<2>&,
    const std::vector<Vec<Vec<double,2>,2>>&, double, double);

template Vec<double,3> 
decide_limit_dir(const Vec<double,3>&, const NearfieldFacets<3>&,
    const std::vector<Vec<Vec<double,3>,3>>&, double, double);
} // end namespace tbem

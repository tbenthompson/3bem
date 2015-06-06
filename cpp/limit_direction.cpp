#include "limit_direction.h"
#include "geometry.h"
#include "nearest_neighbors.h"
#include "gte_wrapper.h"

namespace tbem {

template <size_t dim>
Vec<double,dim> decide_limit_dir(const Vec<double,dim>& pt,
    const NearestFacets<dim>& nearest_facets, double epsilon, 
    double safety_factor) 
{
    const auto& closest_facet = nearest_facets.nearest_facet;
    // TODO: The length scale here needs to be exactly double the length scale
    // used in finding the nearest facets, which means that the functions should
    // be unified. Basically, the function to find nearest_facets deserves to be
    // in this file.
    double len_scale = facet_ball(closest_facet).radius * 2;

    // If the point isn't lying precisely on any facet, then no limit 
    // needs to be taken, in which case the limit direction is the 0 vector.
    double threshold = len_scale * epsilon;
    if (nearest_facets.distance > threshold) {
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
    for (auto f: nearest_facets.facets) {
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

template 
Vec<double,2> decide_limit_dir(const Vec<double,2>& pt,
    const NearestFacets<2>& nearest_facets, double epsilon, double safety_factor);
template 
Vec<double,3> decide_limit_dir(const Vec<double,3>& pt,
    const NearestFacets<3>& nearest_facets, double epsilon, double safety_factor);

} // end namespace tbem

#include "limit_direction.h"
#include "geometry.h"
#include "nearest_neighbors.h"
#include "boost_geometry_wrapper.h"

namespace tbem {

template <size_t dim>
Vec<double,dim> decide_limit_dir(const Vec<double,dim>& pt,
    const NearestFacets<dim>& nearest_facets, double epsilon) 
{
    const auto& closest_facet = nearest_facets.nearest_facet;
    const double diam_to_radius = 0.5;
    double len_scale = hypot(closest_facet[0] - closest_facet[1]) * diam_to_radius;
    double threshold = len_scale * epsilon;
    if (nearest_facets.distance > threshold) {
        return zeros<Vec<double,dim>>::make();
    }

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

    return end_pt - pt;
}

template 
Vec<double,2> decide_limit_dir(const Vec<double,2>& pt,
    const NearestFacets<2>& nearest_facets, double epsilon);
template 
Vec<double,3> decide_limit_dir(const Vec<double,3>& pt,
    const NearestFacets<3>& nearest_facets, double epsilon);

} // end namespace tbem

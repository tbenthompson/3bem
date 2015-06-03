#include "limit_direction.h"
#include "geometry.h"
#include "nearest_neighbors.h"
#include "boost_geometry_wrapper.h"

namespace tbem {

template <size_t dim>
Vec<double,dim> away_from_facet_dir(const Vec<double,dim>& pt,
    const Vec<Vec<double,dim>,dim>& f) 
{
    auto facet_normal = unscaled_normal(f); 
    auto which_side = which_side_point(f, pt);
    if (which_side == INTERSECT) {
        return facet_normal;
    } else if (which_side == FRONT) {
        return facet_normal;
    } else {
        return -facet_normal;
    }
}

template <size_t dim>
Vec<double,dim> decide_limit_dir(const Vec<double,dim>& pt,
    const NearestFacets<dim>& nearest_facets) 
{
    if (nearest_facets.distance > 0) {
        return zeros<Vec<double,dim>>::make();
    }

    auto close_center = centroid(nearest_facets.nearest_facet);
    auto close_dir = away_from_facet_dir(pt, nearest_facets.nearest_facet);
    std::cout << "CLOSEDIR: " << close_dir << std::endl;
    std::cout << "NORMAL: " << unscaled_normal(nearest_facets.nearest_facet) << std::endl;
    auto end_pt = close_center + close_dir;

    // Check if the vector intersects any facets.
    // If it does, back up to halfway between the intersection point and the
    // singular point
    for (auto f: nearest_facets.facets) {
        assert(f != nearest_facets.nearest_facet);
        std::vector<Vec<double,dim>> intersections =
            seg_facet_intersection<dim>(f, {end_pt, pt});
        assert(intersections.size() != 2);
        if (intersections.size() == 1 && intersections[0] != pt) {
            end_pt = pt + (intersections[0] - pt) / 2.0;
        }
    }
    std::cout << end_pt << std::endl;
    std::cout << end_pt - pt << std::endl;
    std::cout << " " << std::endl;

    return end_pt - pt;
}

template 
Vec<double,2> decide_limit_dir(const Vec<double,2>& pt,
    const NearestFacets<2>& nearest_facets);
template 
Vec<double,3> decide_limit_dir(const Vec<double,3>& pt,
    const NearestFacets<3>& nearest_facets);

} // end namespace tbem

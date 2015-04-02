#ifndef __ALNANSNNADNN_FACET_INFO_H
#define __ALNANSNNADNN_FACET_INFO_H

#include "mesh.h"
#include "vec.h"
#include "numerics.h"

namespace tbem {

template <size_t dim>
struct FacetInfo {
    const Facet<dim> face;
    const double area_scale;
    const double length_scale;
    const double jacobian;
    const Vec<double,dim> normal;

    static FacetInfo<dim> build(const Facet<dim>& facet);
};

template <>
inline FacetInfo<3> FacetInfo<3>::build(const Facet<3>& facet) {
    auto unscaled_n = unscaled_normal(facet);
    auto area = tri_area(unscaled_n);
    auto length_scale = std::sqrt(area);
    auto jacobian = area * inv_ref_facet_area<3>();
    auto normal = unscaled_n / jacobian;
    return FacetInfo<3>{facet, area, length_scale, jacobian, normal};
}

template <>
inline FacetInfo<2> FacetInfo<2>::build(const Facet<2>& facet) {
    auto unscaled_n = unscaled_normal(facet);
    auto area_scale = hypot2(unscaled_n);
    auto length = std::sqrt(area_scale);
    auto jacobian = length * inv_ref_facet_area<2>();
    auto normal = unscaled_n / length;
    return FacetInfo<2>{facet, area_scale, length, jacobian, normal};
}

template <size_t dim>
std::vector<FacetInfo<dim>> get_facet_info(const Mesh<dim>& m) {
    std::vector<FacetInfo<dim>> out;
    for (const auto& f: m.facets) {
        out.push_back(FacetInfo<dim>::build(f));
    }
    return out;
}

} //end namespace tbem

#endif

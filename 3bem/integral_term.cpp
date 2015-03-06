#include "integral_term.h"

namespace tbem {

template <size_t dim>
NearestPoint<dim> FarNearLogic<dim>::decide(const Vec<double,dim>& pt,
    const FacetInfo<dim>& facet) 
{
    auto near_ref_pt = closest_pt_facet(pt, facet.face);
    auto near_pt = ref_to_real(near_ref_pt, facet.face);
    auto exact_dist2 = dist2(near_pt, pt);
    auto appx_dist2 = exact_dist2;//appx_face_dist2(obs.loc, facet.face);
    bool nearfield = appx_dist2 < pow(far_threshold, 2) * facet.area_scale;
    if (nearfield) {
        bool singular = appx_dist2 < singular_threshold * facet.area_scale;
        if (singular) { 
            return {
                near_ref_pt, near_pt, std::sqrt(exact_dist2), FarNearType::Singular
            };
        } else {
            return {
                near_ref_pt, near_pt, std::sqrt(exact_dist2), FarNearType::Nearfield
            };
        }
    } else {
        return {near_ref_pt, near_pt, std::sqrt(exact_dist2), FarNearType::Farfield};
    }
}

template struct FarNearLogic<2>;
template struct FarNearLogic<3>;

} //end namespace tbem

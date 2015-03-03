#ifndef __ALSDJLA123_OBS_PT_H
#define __ALSDJLA123_OBS_PT_H

#include "vec.h"
#include "facet_info.h"

namespace tbem {

template <size_t dim>
struct ObsPt {
    static ObsPt<dim> from_face(const Vec<double,dim-1>& ref_loc,
        const FacetInfo<dim>& obs_face)
    {
        const int basis_order = 1;
        return {
            obs_face.length_scale / basis_order,
            ref_to_real(ref_loc, obs_face.face),
            obs_face.normal,
            obs_face.normal 
        };
    }

    const double len_scale;
    const Vec<double,dim> loc;
    const Vec<double,dim> normal;
    const Vec<double,dim> richardson_dir;
};


} //end namespace tbem
#endif

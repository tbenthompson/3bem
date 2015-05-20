#ifndef __ALSDJLA123_OBS_PT_H
#define __ALSDJLA123_OBS_PT_H

#include "vec.h"
#include "facet_info.h"
#include "nearest_neighbors.h"

namespace tbem {

template <size_t dim>
struct ObsPt {
    static ObsPt<dim> away_from_nearest_facets(const Vec<double,dim-1>& ref_loc,
        const FacetInfo<dim>& obs_face, const Mesh<dim>& mesh)
    {
        auto loc = ref_to_real(ref_loc, obs_face.face);
        auto nf = nearest_facets(loc, mesh.facets);
        double interior_shift = facet_ball(obs_face.face).radius * 1e-8;
        auto slightly_interior_loc = loc + interior_shift * obs_face.normal;
        auto rich_dir = decide_richardson_dir(slightly_interior_loc, nf);
        auto rich_length = hypot(rich_dir);
        return {rich_length, loc, obs_face.normal, rich_dir / rich_length};
    }

    const double len_scale;
    const Vec<double,dim> loc;
    const Vec<double,dim> normal;
    const Vec<double,dim> richardson_dir;

    ObsPt(double len_scale, const Vec<double,dim>& loc,
          const Vec<double,dim>& normal, const Vec<double,dim>& richardson_dir):
        len_scale(len_scale), loc(loc), normal(normal), richardson_dir(richardson_dir)
    {}
};


} //end namespace tbem
#endif

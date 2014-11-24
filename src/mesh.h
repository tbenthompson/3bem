#ifndef __MESH_3D_H
#define __MESH_3D_H
#include <vector>
#include <array>
#include <unordered_map>
#include <iostream>
#include <algorithm>
#include <GL/gl.h>
#include "vec.h"
#include "constraint.h"

template <int dim>
struct Facet {
    const Vec<Vec<double,dim>,dim> vertices;
};

template <int dim>
struct Mesh {
    typedef std::function<Vec<double,dim>(Vec<double,dim>)> RefineFnc;
    const std::vector<Facet<dim>> facets;
    const bool has_refine_mod;
    const RefineFnc refine_mod;

    Mesh<dim> refine(const std::vector<int>& refine_these) const;
    Mesh<dim> refine() const;
    Mesh<dim> refine_repeatedly(unsigned int times) const;

    static Mesh<dim> from_vertices_faces(const std::vector<Vec<double,dim>>& vertices,
                         const std::vector<std::array<int,dim>>& facets,
                         bool has_refine_mod = false,
                         const typename Mesh<dim>::RefineFnc& refine_mod = nullptr);
};

#endif

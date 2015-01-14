#ifndef __BLAHBLABHLBHALBHALBHA_FUNCTION_H
#define __BLAHBLABHLBHALBHALBHA_FUNCTION_H
#include "vec.h"

namespace tbem {

template <typename T, size_t dim>
struct FacetFunction {
    const Vec<T,dim> vertices;
};

template <typename T, size_t dim>
struct FunctionDOFIterator;

template <typename T, size_t dim>
struct Function {
    const std::vector<FacetFunction<T,dim>> facets;

    size_t n_facets() const;
    size_t n_dofs() const;

    FunctionDOFIterator<T,dim> begin() const;
    FunctionDOFIterator<T,dim> end() const;

    Function<T,dim> refine(const std::vector<int>& refine_these) const;
    Function<T,dim> refine() const;
    Function<T,dim> refine_repeatedly(unsigned int times) const;

    static Function<T,dim> form_union(const std::vector<Function<T,dim>>& others);

    static Function<T,dim> from_vertices_faces(const std::vector<T>& vertices,
        const std::vector<std::array<int,dim>>& facets_by_vert_idx);
};

} // End namespace tbem

#endif

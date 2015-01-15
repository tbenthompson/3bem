#ifndef __BLAHBLABHLBHALBHALBHA_FUNCTION_H
#define __BLAHBLABHLBHALBHALBHA_FUNCTION_H
#include "vec_ops.h"

namespace tbem {

template <typename T, size_t dim>
struct FunctionDOFIterator;

template <typename T, size_t dim>
struct Function {
    const std::vector<Vec<T,dim>> facets;

    size_t n_facets() const;
    size_t n_dofs() const;

    FunctionDOFIterator<T,dim> begin() const;
    FunctionDOFIterator<T,dim> end() const;

    Function<T,dim> refine(const std::vector<int>& refine_these) const;
    Function<T,dim> refine() const;
    Function<T,dim> refine_repeatedly(unsigned int times) const;

    static Function<T,dim> create_union(const std::vector<Function<T,dim>>& others);

    static Function<T,dim> from_vertices_faces(const std::vector<T>& vertices,
        const std::vector<std::array<int,dim>>& facets_by_vert_idx);
};

template <size_t dim, typename T>
Function<T,dim> vec_to_fnc(const std::vector<T>& dofs) {
    return {reinterpret_vector<Vec<T,dim>>(dofs)};
}

template <size_t dim, typename T>
std::vector<T> fnc_to_vec(const Function<T,dim>& fnc) {
    return reinterpret_vector<T>(fnc.facets);
}

} // End namespace tbem

#endif

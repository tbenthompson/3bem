#ifndef __INTTEST_DISLOC_SHARED_H
#define __INTTEST_DISLOC_SHARED_H

#include "3bem.h"
#include "elastic_kernels.h"
#include "constraint.h"

using namespace tbem;

//TODO: Refactor more of rect_dislocation and planestrain_fault into here

template <size_t dim>
ConstraintMatrix surf_fault_constraints(const VertexIterator<dim>& surf_it,
    const VertexIterator<dim>& fault_it) 
{
    auto continuity = mesh_continuity(surf_it);
    auto cut_cont = cut_at_intersection(continuity, surf_it, fault_it);
    auto constraints = convert_to_constraints(cut_cont);
    auto constraint_matrix = from_constraints(constraints);
    return constraint_matrix;
}

#endif

#ifndef __INTTEST_DISLOC_SHARED_H
#define __INTTEST_DISLOC_SHARED_H

#include "3bem.h"
#include "elastic_kernels.h"
#include "constraint_matrix.h"

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

template <size_t dim>
std::vector<double> solve(
    const BlockOperator& lhs,
    const Function& rhs, 
    const BlockDOFMap& dof_map,
    const Mesh<dim>& mesh,
    const ConstraintMatrix& constraint_matrix) {

    size_t count = 0;
    return solve_system(rhs, 1e-5,
        [&] (std::vector<double>& x, std::vector<double>& y) {
            std::cout << "iteration " << count << std::endl;
            count++;

            auto x_fncs = expand(dof_map, x);
            BlockFunction x_vec(dim);
            for (size_t d = 0; d < dim; d++) {
                x_vec[d] = distribute_vector(constraint_matrix, x_fncs[d], mesh.n_dofs());
            }
            auto y_vec = apply_operator(lhs, x_vec);
            BlockFunction condensed(dim);
            for (size_t d = 0; d < dim; d++) {
                condensed[d] = condense_vector(constraint_matrix, y_vec[d]);
            }
            auto out = concatenate(dof_map, condensed);
            for (std::size_t i = 0; i < out.size(); i++) {
                y[i] = out[i];
            }
        }
    );
}

#endif

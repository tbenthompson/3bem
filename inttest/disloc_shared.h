#ifndef __INTTEST_DISLOC_SHARED_H
#define __INTTEST_DISLOC_SHARED_H

#include "3bem.h"
#include "elastic_kernels.h"
#include "constraint_matrix.h"
#include "operator.h"

using namespace tbem;

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
BlockVectorX make_rhs(const Mesh<dim>& surface, const Mesh<dim>& fault,
        const ConstraintMatrix& constraint_matrix, const QuadStrategy<dim>& qs,
        const ElasticHypersingular<dim>& hyp, const BlockVectorX& du) {
    TIC
    auto p_rhs = make_problem<dim>(surface, fault, hyp);
    auto rhs_op = make_matrix_free(p_rhs, qs);
    auto all_dofs_rhs = rhs_op.apply(du);

    BlockVectorX condensed(dim);
    for (size_t d = 0; d < dim; d++) {
        condensed[d] = condense_vector(constraint_matrix, all_dofs_rhs[d]);
    };
    TOC("Building RHS");
    return condensed;
}

template <size_t dim>
MatrixFreeFarfieldOperator<dim,ElasticHypersingular<dim>>
make_lhs(const ElasticHypersingular<dim>& hyp, const Mesh<dim>& mesh,
    const QuadStrategy<dim>& qs )
{
    TIC
    auto p_lhs = make_problem<dim>(mesh, mesh, hyp);
    auto lhs = make_matrix_free(p_lhs, qs);
    TOC("Building LHS matrices");
    return lhs;
}


template <size_t dim>
struct DislocationProblem {
    MatrixFreeFarfieldOperator<dim,ElasticHypersingular<dim>> lhs;
    ConstraintMatrix constraint_matrix;
    ElasticHypersingular<dim> hyp;
    Mesh<dim> surface;
    Mesh<dim> fault;
    BlockVectorX condensed_rhs;
    BlockDOFMap dof_map;
    VectorX rhs;

    DislocationProblem(ElasticHypersingular<dim>& hyp, 
            const Mesh<dim>& surface, const Mesh<dim>& fault,
            const QuadStrategy<dim>& qs, const BlockVectorX& du):
        lhs(make_lhs(hyp, surface, qs)),
        constraint_matrix(surf_fault_constraints(surface.begin(), fault.begin())),
        hyp(hyp),
        surface(surface),
        fault(fault),
        condensed_rhs(make_rhs(surface, fault, constraint_matrix, qs, hyp, du)),
        dof_map(block_dof_map_from_functions(condensed_rhs)),
        rhs(concatenate(dof_map, condensed_rhs))
    {}

    void solve()
    {
        size_t count = 0;
        auto reduced_soln = solve_system(rhs, 1e-5,
            [&] (std::vector<double>& x, std::vector<double>& y) {
                std::cout << "iteration " << count << std::endl;
                count++;

                auto x_fncs = expand(dof_map, x);
                BlockVectorX x_vec(dim);
                for (size_t d = 0; d < dim; d++) {
                    x_vec[d] = distribute_vector(
                        constraint_matrix, x_fncs[d], surface.n_dofs()
                    );
                }
                auto y_vec = lhs.apply(x_vec);
                BlockVectorX condensed(dim);
                for (size_t d = 0; d < dim; d++) {
                    condensed[d] = condense_vector(constraint_matrix, y_vec[d]);
                }
                auto out = concatenate(dof_map, condensed);
                for (std::size_t i = 0; i < out.size(); i++) {
                    y[i] = out[i];
                }
            }
        );

        auto disp_reduced_vec = expand(dof_map, reduced_soln);

        BlockVectorX soln(dim);
        for (size_t d = 0; d < dim; d++) {
            soln[d] = distribute_vector(
                constraint_matrix, disp_reduced_vec[d], surface.n_dofs()
            );
        };

        auto file = HDFOutputter("test_out/rect_dislocation_u.hdf5");
        out_surface(file, surface, soln);
    }
};

#endif

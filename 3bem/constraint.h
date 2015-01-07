#ifndef __UYIWRKSLSHHHS_CONSTRAINT_H
#define __UYIWRKSLSHHHS_CONSTRAINT_H
#include <utility>
#include <vector>
#include <unordered_map>
#include "vec.h"

namespace tbem {

/* A list of matrix constraints is generated from the mesh
 * connectivity and boundary conditions. These constraints
 * are represented by an integer referring to the relevant
 * degree of freedom and a double that multiplies that 
 * dofs value in the linear system. A rhs value is represented
 * by a dof id of -1. 
 * So, for example:
 * 4x + 3y = 2 ====> 4x + 3y - 2 = 0
 * would be a constraint list consisting of the tuples
 * [(x_dof, 4), (y_dof, 3), (RHS, -2)]
 *
 * Constraints are used to ensure continuity between
 * elements and to enforce displacement discontinuities between elements
 * that connect to another element with a displacement discontinuity 
 * boundary condition or solution type.
 *
 * Continuity constraints for a Lagrange interpolating basis that has a
 * node at the boundary will simply be of the form [(e1_dof, 1), (e2_dof, -1)]
 * because the two boundary dofs are set equal.
 */ 
template <typename T>
using DOFWeight = std::pair<int,T>;

/* A constraint is stored as a list of DOFs and weights, where the
* first dof is implicitly the constrained dof.
*/
struct Constraint {
    std::vector<DOFWeight<double>> dof_weights;
    double rhs_value;
    friend std::ostream& operator<<(std::ostream& os, const Constraint& c);
};

/* Store a constraint matrix as a map from the constrained dof to
 * the constraint on that dof. The constrained dof will be the last
 * one in the constraint. This allows looping through a vector and
 * applying constraints without worrying about cycles (IS THIS SENTENCE TRUE?).
 */
struct ConstraintMatrix {
    typedef std::unordered_map<int,Constraint> MapT;
    const MapT c_map;

    /* Accepts a reduced DOF vector and returns the full DOF vector. */
    template <typename T>
    std::vector<T> get_all(const std::vector<T>& in, int total_dofs) const; 

    /* Accepts a full DOF vector and returns the reduced DOF vector.
     */
    template <typename T>
    std::vector<T> get_reduced(const std::vector<T>& all) const;

    /* Adds to an entry in a vector according to the constraint matrix.
     * For example, if DOF #2 is constrained and I add to DOF #2, the value is 
     * applied to those DOFs that #2 is constrained to be equal to.
     */
    template <typename T>
    void add_vec_with_constraints(const DOFWeight<T>& entry,
                                  std::vector<T>& rhs) const;

    ConstraintMatrix add_constraints(const std::vector<Constraint>& constraints);

    bool is_constrained(int dof) const;

    friend std::ostream& operator<<(std::ostream& os, const ConstraintMatrix& cm);

    static 
    ConstraintMatrix from_constraints(const std::vector<Constraint>& constraints);
};

std::vector<DOFWeight<double>>::iterator find_unconstrained_dof(
        const ConstraintMatrix::MapT& c_map, Constraint& constraint);

bool is_constrained(const ConstraintMatrix::MapT& c_map, int dof);

/* Constrain two degrees of freedom to be identical. */
Constraint continuity_constraint(int dof1, int dof2);

/* Set a boundary condition. A very simple constraint that looks like:
 * value[dof] = value
 */
Constraint boundary_condition(int dof, double value);

template <typename T, unsigned long dim>
using Vec = std::array<T,dim>;
template <typename T, int dim>
struct MeshField;
template <int dim>
using Mesh = MeshField<Vec<double,dim>,dim>;

/* Find the overlapping vertices for the given mesh and produce continuity
 * constraints. 
 */
template <int dim>
std::vector<Constraint> mesh_continuity(const Mesh<dim>& m, double eps = 1e-10);

/* Find the vertices overlapping between a discontinuity mesh and a surface
 * on which continuity constraints have been applied. 
 * Note that continuity constraints and ONLY continuity constraints should
 * have been applied prior to using this function.
 */
template <int dim> 
ConstraintMatrix apply_discontinuities(const Mesh<dim>& surface,
                                       const Mesh<dim>& disc,
                                       const ConstraintMatrix& c_mat,
                                       double eps = 1e-10);


template <typename T>
std::vector<T> ConstraintMatrix::get_all(const std::vector<T>& in,
                                              int total_dofs) const {
    std::vector<T> out(total_dofs); 
    int next_in = 0;
    for (int i = 0; i < total_dofs; i++) {
        auto dof_and_constraint = c_map.find(i);
        if (dof_and_constraint == c_map.end()) {
            out[i] = in[next_in];
            next_in++;
            continue;
        }

        auto dof_constraint = dof_and_constraint->second.dof_weights;
        auto out_val = constant<T>::make(dof_and_constraint->second.rhs_value);
        for (std::size_t j = 0; j < dof_constraint.size() - 1; j++) {
            out_val -= dof_constraint[j].second * out[dof_constraint[j].first];
        }
        auto this_dof_weight = (dof_constraint.end() - 1)->second;
        out[i] = out_val / this_dof_weight;
    }
    return out;
}

template <typename T>
void ConstraintMatrix::add_vec_with_constraints(const DOFWeight<T>& entry,
                                                std::vector<T>& vec) const {
    auto dof_and_constraint = c_map.find(entry.first);
    if (dof_and_constraint == c_map.end()) {
        vec[entry.first] += entry.second;
        return;
    }

    auto constraint = dof_and_constraint->second;
    auto dof_weights = constraint.dof_weights;
    int n_dof_weights = dof_weights.size();
    double this_weight = dof_weights[n_dof_weights - 1].second;
    for (int i = 0; i < n_dof_weights - 1; i++) {
        T recurse_weight = -dof_weights[i].second * entry.second / this_weight;
        DOFWeight<T> new_entry{dof_weights[i].first, recurse_weight};
        add_vec_with_constraints(new_entry, vec);
    }
}

template <typename T>
std::vector<T> ConstraintMatrix::get_reduced(const std::vector<T>& all) const {
    std::vector<T> condensed_dofs(all.size(), zeros<T>::make());
    for (std::size_t i = 0; i < all.size(); i++) {
        add_vec_with_constraints(DOFWeight<T>{i, all[i]}, condensed_dofs);
    }

    std::vector<T> out;
    for (std::size_t i = 0; i < all.size(); i++) {
        auto dof_and_constraint = c_map.find(i);
        if (dof_and_constraint == c_map.end()) {
            out.push_back(condensed_dofs[i]);
            continue;
        }
    }

    return out;
}

} //END NAMESPACE tbem
#endif

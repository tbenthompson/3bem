#ifndef __UYIWRKSLSHHHS_CONSTRAINT_H
#define __UYIWRKSLSHHHS_CONSTRAINT_H
#include <utility>
#include <vector>
#include <unordered_map>

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
typedef std::pair<int, double> DOFWeight;

/* A constraint is stored as a list of DOFs and weights, where the
* first dof is implicitly the constrained dof.
*/
struct Constraint {
    std::vector<DOFWeight> dof_constraints;
    double rhs_value;
    friend std::ostream& operator<<(std::ostream& os, const Constraint& c);
};

/* Store a constraint matrix as a map from the constrained dof to
 * the constraint on that dof. The constrained dof will be the last
 * one in the constraint. This allows looping through a vector and
 * applying constraints without worrying about cycles.
 */
struct ConstraintMatrix {
    typedef std::unordered_map<int,Constraint> MapT;
    const MapT c_map;

    ConstraintMatrix add_constraints(const std::vector<Constraint>& constraints);

    /* Accepts a reduced DOF vector and returns a full DOF vector
     */
    std::vector<double> get_all(const std::vector<double>& in, int total_dofs) const; 

    /* Accepts a full DOF vector and returns the reduced DOF vector.
     */
    std::vector<double> get_reduced(const std::vector<double>& all) const;

    void add_vec_with_constraints(const DOFWeight& entry, std::vector<double>& rhs) const;
    std::vector<double> condense(const std::vector<double>& all) const;

    friend std::ostream& operator<<(std::ostream& os, const ConstraintMatrix& cm);

    static ConstraintMatrix from_constraints(const std::vector<Constraint>& constraints);
};

/* Constrain two degrees of freedom to be identical. */
Constraint continuity_constraint(int dof1, int dof2);

/* This creates a constraint representing the offset between the value
 * of two DOFs. Primarily useful for imposing a fault slip as the 
 * displacement jump when a fault intersects the surface.
 * The direction of offset is from dof1 to dof2. So,
 * value[dof1] - value[dof2] = -offset
 */
Constraint offset_constraint(int dof1, int dof2, double offset);

/* Set a boundary condition. A very simple constraint that looks like:
 * value[dof] = value
 */
Constraint boundary_condition(int dof, double value);

template <int dim>
class Mesh;

/* Find the overlapping vertices for the given mesh and produce continuity
 * constraints. 
 * TODO: Handle hanging nodes?
 */
template <int dim>
std::vector<Constraint> mesh_continuity(const Mesh<dim>& m, double eps = 1e-10);

/* Find the vertices overlapping between a discontinuity mesh and a surface
 * on which continuity constraints have been applied. 
 * Note that continuity constraints and ONLY continuity constraints should
 * have been applied prior to using this function.
 * TODO: Handle more general cases where its not just the vertices that 
 * overlap.
 */
template <int dim> 
ConstraintMatrix apply_discontinuities(const Mesh<dim>& surface,
                                       const Mesh<dim>& disc,
                                       const ConstraintMatrix& c_mat,
                                       double eps = 1e-10);
#endif

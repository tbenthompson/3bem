// #include "common.h"
// #include "constraint.h"
// #include <armadillo>
// 
// namespace codim1 {
// 
// Constraint continuity_constraint(int dof1, int dof2) {
//     return {
//         DOFWeight(dof1, 1.0),
//         DOFWeight(dof2, -1.0)
//     };
// }
// 
// Constraint offset_constraint(int dof1, int dof2, double offset) {
//     Constraint c = continuity_constraint(dof1, dof2);
//     c.push_back(DOFWeight(RHS, offset));
//     return c;
// }
// 
// void add_constraint(ConstraintMatrix& cm, Constraint c) {
//     int constrained = c[0].first;
//     cm[constrained] = c;
// }
// 
// void non_constrained_row(MatrixEntry entry,
//                          mat& mat, vec& rhs,
//                          ConstraintMatrix cm) {
// 
//     auto dof_and_constraint = cm.find(entry.col);
//     if (dof_and_constraint == cm.end()) {
//         mat(entry.row, entry.col) += entry.value;
//         return;
//     }
// 
//     Constraint constraint = dof_and_constraint->second;
//     for(unsigned int i = 1; i < constraint.size(); i++) {
//         MatrixEntry new_entry({entry.row,
//                                constraint[i].first, 
//                                -constraint[i].second * entry.value});
//         add_mat_with_constraints(new_entry, mat, rhs, cm);
//     }
// }
// 
// void constrained_row(MatrixEntry entry,
//                      mat& mat, vec& rhs,
//                      ConstraintMatrix cm) {
//     Constraint constraint = cm[entry.row];
//     for(unsigned int i = 1; i < constraint.size(); i++) {
//         MatrixEntry new_entry({constraint[i].first, 
//                                entry.col,
//                                -constraint[i].second * entry.value});
//         add_mat_with_constraints(new_entry, mat, rhs, cm);
//     }
// }
// 
// void add_mat_with_constraints(MatrixEntry entry,
//                               mat& mat, vec& rhs,
//                               ConstraintMatrix cm) {
//     // TODO: Implement for inhomogenous constraints.
//     DBGMSG(std::cout, "Add to matrix with row dof: " << entry.row <<
//                       " and col dof: " << entry.col << 
//                       " and value: " << entry.value);
// 
//     // The only way that entry.row == RHS is if this function is called
//     // after recursively applying 
//     if (entry.row == RHS) {
//         return;
//     }
// 
//     if (entry.col == RHS) {
//         add_rhs_with_constraints(DOFWeight({entry.row, entry.value}),
//                                  rhs, cm);
//         return;
//     }
// 
//     auto constraint = cm.find(entry.row);
//     if (constraint == cm.end()) {
//         non_constrained_row(entry, mat, rhs, cm);
//         return;
//     }
//     constrained_row(entry, mat, rhs, cm);
//     return;
// }
// 
// /* TODO: This could be consolidated with the add_mat... function */
// void add_rhs_with_constraints(DOFWeight entry,
//                               vec& rhs, 
//                               ConstraintMatrix cm) {
//     DBGMSG(std::cout, "Add to RHS with entry dof: " << entry.first <<
//                       " and value: " << entry.second);
// 
//     auto dof_and_constraint = cm.find(entry.first);
//     if (dof_and_constraint == cm.end()) {
//         rhs(entry.first) += entry.second;
//         return;
//     }
// 
//     Constraint constraint = dof_and_constraint->second;
//     for(unsigned int i = 1; i < constraint.size(); i++) {
//         DOFWeight new_entry({constraint[i].first, 
//                              -constraint[i].second * entry.second});
//         add_rhs_with_constraints(new_entry, rhs, cm);
//     }
// }
// 
// 
// } // end namespace codim1

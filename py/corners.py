from tbempy import *
from tbempy.TwoD import *
from laplace2d import theta_u, make_theta_dudn
import scipy.sparse.linalg as sp_la
import numpy as np

center = [5, 0]
theta_dudn = make_theta_dudn(center)
r = 3.0
refine_level = 7;
linear_solve_tol = 1e-9;
switch_dof = 2000
qs = QuadStrategy(4, 4, 10, 3.0, 1e-5);
surface = circle_mesh(center, r, refine_level);
n_dofs = surface.n_dofs();
if switch_dof > n_dofs:
    switch_dof = n_dofs

continuity = mesh_continuity(surface.begin());
pot_constraints = convert_to_constraints(continuity)
pot_bc_constraints = interpolate_bc_constraints(
    surface, range(0, switch_dof), theta_u
)
pot_constraints.extend(pot_bc_constraints)
pot_cm = from_constraints(pot_constraints)

flux_constraints = interpolate_bc_constraints(
    surface, range(switch_dof, n_dofs), theta_dudn
)
flux_cm = from_constraints(flux_constraints)

single_kernel = LaplaceSingle()
double_kernel = LaplaceDouble()
id = IdentityScalar()

single_layer = make_boundary_integral(surface, surface, single_kernel)
single_op = mesh_to_mesh_operator(single_layer, qs)

double_layer = make_boundary_integral(surface, surface, double_kernel)
double_op = mesh_to_mesh_operator(double_layer, qs)

mass = make_boundary_integral(surface, surface, id)
mass_op = mass_operator(mass, qs)

rhs_pot = VectorX([0.0] * n_dofs)
rhs_flux = VectorX([0.0] * n_dofs)
condensed = BlockVectorX([
    condense_vector(pot_cm, rhs_pot),
    condense_vector(flux_cm, rhs_flux)
])

dof_map = block_dof_map_from_functions(condensed)
rhs = concatenate(dof_map, condensed)

def mv(v):
    vec_v = VectorX(v)
    both = expand(dof_map, vec_v)
    pot = distribute_vector(pot_cm, both.storage[0], surface.n_dofs())
    flux = distribute_vector(flux_cm, both.storage[1], surface.n_dofs())
    single_eval = single_op.apply_scalar(flux)
    double_eval = double_op.apply_scalar(pot)
    mass_eval = mass_op.apply_scalar(pot)
    out = mass_eval - single_eval + double_eval
    pot_out = condense_vector(pot_cm, out)
    flux_out = condense_vector(flux_cm, out)
    print("Iteration: " + str(mv.it))
    mv.it += 1
    res = np.array(concatenate(dof_map, BlockVectorX([pot_out, flux_out])).storage)
    return res
mv.it = 0

np_rhs = np.array(rhs.storage)
np_rhs[0] += 2 * linear_solve_tol
A = sp_la.LinearOperator((np_rhs.shape[0], np_rhs.shape[0]),
                         matvec = mv, dtype = np.float64)
def callback(x): print(np.max(np.abs(x)))
res = sp_la.gmres(A, np_rhs, tol = linear_solve_tol, restart = 200, callback = callback)
both_solns = expand(dof_map, VectorX(res[0]))
soln_pot = distribute_vector(pot_cm, both_solns.storage[0], surface.n_dofs())
soln_flux = distribute_vector(flux_cm, both_solns.storage[1], surface.n_dofs())

from antiplane_fault import get_vertices
import matplotlib.pyplot as plt
xs = get_vertices(surface, 0)
plt.plot(xs, np.array(soln_pot.storage))
plt.show()

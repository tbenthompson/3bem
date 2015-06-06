from tbempy.TwoD import *
from dislocation import *
import numpy as np
import scipy.sparse.linalg as sp_la

def full_space():
    fault = line_mesh([0, -1], [0, 0])
    slip = np.ones(fault.n_dofs())
    surface = line_mesh([-10, -0.5], [10, -0.5]).refine_repeatedly(8)

    qs = QuadStrategy(5, 8, 4.0, 1e-13);

    double_kernel = LaplaceDouble()
    mthd = make_adaptive_integration_mthd(qs, double_kernel)
    double_layer = integral_operator(surface, fault, mthd, fault)

    def fnc(x):
        if x[0] == 0.0:
            return 0.0
        op = dense_interior_operator([x], [[0, 1]], fault, mthd, fault)
        val = op.apply(slip)
        exact = np.arctan(0.5 / x[0]) / np.pi
        assert(np.abs(exact - val) < 2e-4)
        return val[0]
    u = interpolate(surface, fnc)

def get_qs():
    return QuadStrategy(5, 10, 3.0, 1e-5)

def solve_half_space(slip, fault, surface):
    constraint_matrix = faulted_surface_constraints(tbempy.TwoD, surface, fault, 1)
    qs = get_qs()

    all_mesh = Mesh.create_union([surface, fault])

    hypersingular_kernel = LaplaceHypersingular()
    hypersingular_mthd = make_adaptive_integration_mthd(qs, hypersingular_kernel)
    rhs_op = integral_operator(surface, fault, hypersingular_mthd, all_mesh)
    full_rhs = (rhs_op.apply(slip))
    rhs = condense_vector(constraint_matrix, full_rhs)

    lhs_op = integral_operator(surface, surface, hypersingular_mthd, all_mesh)

    def mv(v):
        distributed = distribute_vector(constraint_matrix, v, surface.n_dofs())
        applied = lhs_op.apply(distributed)
        condensed = condense_vector(constraint_matrix, applied)
        mv.it += 1
        return condensed
    mv.it = 0

    A = sp_la.LinearOperator((rhs.shape[0], rhs.shape[0]),
                             matvec = mv, dtype = np.float64)
    res = sp_la.gmres(A, rhs, tol = 1e-6)
    soln = distribute_vector(constraint_matrix, res[0], surface.n_dofs())
    return soln

def half_space(refine):
    fault = line_mesh([0, -1], [0, 0])
    slip = np.ones(fault.n_dofs())
    surface = line_mesh([50, 0.0], [-50, 0.0]).refine_repeatedly(refine)
    soln = solve_half_space(slip, fault, surface)
    xs = surface.facets[:, :, 0].reshape((surface.n_facets() * 2))
    indices = [i for i in range(len(xs)) if 0 < np.abs(xs[i]) < 10]
    xs = xs[indices]
    exact = np.arctan(1.0 / xs) / np.pi
    error = np.sqrt(np.mean((exact - soln[indices]) ** 2))
    return error

def half_space_interior(refine):
    fault = line_mesh([0, -1], [0, 0])
    slip = np.ones(fault.n_dofs())
    surface = line_mesh([50, 0.0], [-50, 0.0]).refine_repeatedly(refine)
    qs = get_qs()
    all_mesh = Mesh.create_union([surface, fault])

    soln = solve_half_space(slip, fault, surface)

    xs = np.linspace(-5, 5, 100)
    ys = np.linspace(-5, 0, 100)

    normalsx = np.array([[1, 0]] * xs.shape[0] * ys.shape[0])
    normalsy = np.array([[0, 1]] * xs.shape[0] * ys.shape[0])
    pts = []
    for x in xs:
        for y in ys:
            pts.append([x, y])
    pts = np.array(pts)

    double_kernel = LaplaceDouble()
    hypersingular_kernel = LaplaceHypersingular()
    double_mthd = make_adaptive_integration_mthd(qs, double_kernel)
    hypersingular_mthd = make_adaptive_integration_mthd(qs, hypersingular_kernel)

    mtpo = dense_interior_operator
    disp = mtpo(pts, normalsx, fault, double_mthd, all_mesh).apply(slip) -\
           mtpo(pts, normalsx,  surface, double_mthd, all_mesh).apply(soln)

    tracx = mtpo(pts, normalsx, fault, hypersingular_mthd, all_mesh).apply(slip)\
        - mtpo(pts, normalsx, surface, hypersingular_mthd, all_mesh).apply(soln)
    tracy = mtpo(pts, normalsy, fault, hypersingular_mthd, all_mesh).apply(slip)\
        - mtpo(pts, normalsy, surface, hypersingular_mthd, all_mesh).apply(soln)

    x = np.array([p[0] for p in pts])
    y = np.array([p[1] for p in pts])

    d = 1
    s = 1
    shear_modulus = 1.0
    exact_uz = (-s / (2 * np.pi)) * (
        np.arctan((y - d) / x) -
        np.arctan((y + d) / x))
    exact_tracx = -(s * shear_modulus) / (2 * np.pi) * (
        ((y + d) / (x ** 2 + (y + d) ** 2)) -
        ((y - d) / (x ** 2 + (y - d) ** 2)))
    exact_tracy = (s * shear_modulus) / (2 * np.pi) * (
        (x / (x ** 2 + (y + d) ** 2)) -
        (x / (x ** 2 + (y - d) ** 2)))
    l2_error_disp = np.sqrt(np.mean((disp - exact_uz) ** 2))
    l2_error_tracx = np.sqrt(np.mean((tracx - exact_tracx) ** 2))
    l2_error_tracy = np.sqrt(np.mean((tracy - exact_tracy) ** 2))

    return l2_error_disp, l2_error_tracx, l2_error_tracy

def test_fullspace():
    full_space()

def test_halfspace():
    error = half_space(5)
    assert(error < 0.025)

def test_halfspace_interior():
    disp_error, tracx_error, tracy_error = half_space_interior(7)
    assert(disp_error < 6e-4)
    assert(tracx_error < 2e-3)
    assert(tracy_error < 2e-3)

def test_convergence():
    prev = half_space(4)
    for i in range(5,8):
        error = half_space(i)
        ratio = prev / error
        assert(ratio > 2.0)

if __name__ == "__main__":
    half_space(7)
    half_space_interior(7)

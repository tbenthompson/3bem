from antiplane_fault import solve_half_space
import tbempy.TwoD as tbem
import numpy as np

def test_refinement_exact_error():
    iterations = refiner(exact_errors)
    assert(iterations == 8)

def test_refinement_estimated_error():
    iterations = refiner(estimated_errors)
    assert(iterations == 10)

def refiner(error_fnc):
    fault = tbem.line_mesh([0, -1], [0, 0])
    slip = np.ones(fault.n_dofs())
    L = 100.0
    surface = tbem.line_mesh([L / 2, 0.0], [-L / 2, 0.0]).refine_repeatedly(1)

    refine_fraction = 0.3
    iterations = 0
    all_steps = []
    while surface.n_facets() < 4000:
        n_facets = surface.n_facets()
        soln = solve_half_space(slip, fault, surface)
        all_steps.append((surface, soln))
        # plot_soln(surface, soln)
        facet_error, _ = error_fnc(slip, fault, surface, soln)
        _, l2_error = exact_errors(slip, fault, surface, soln)
        l2_error /= L
        if (l2_error < 0.001):
            break

        worst_facets = sorted(range(n_facets), key = lambda k: facet_error[k])
        n_refine = int(np.ceil(n_facets * refine_fraction))
        refine_me = worst_facets[-n_refine:]

        surface = surface.refine(refine_me)
        iterations += 1
    # animate(all_steps)
    return iterations

def exact_errors(slip, fault, surface, soln):
    xs, exact = get_exact(surface)
    dof_error = np.abs(exact - soln)
    return local_and_global_errors(dof_error, surface)

def get_exact(surface):
    xs = surface.facets[:, :, 0].reshape((surface.n_facets() * 2))
    xs = np.where(xs != 0, xs, 0.00001)
    exact = np.arctan(1.0 / xs) / np.pi
    return xs, exact

def estimated_errors(slip, fault, surface, soln):
    all_mesh = tbem.Mesh.create_union([surface, fault])

    double_mthd = tbem.make_adaptive_integrator(
        1e-5, 6, 6, 5, 10, 3.0, tbem.LaplaceDouble()
    )
    fmm_config = tbem.FMMConfig(0.3, 30, 250, 0.05, True)
    rhs_op = tbem.boundary_operator(surface, fault, double_mthd, fmm_config, all_mesh)
    rhs = rhs_op.apply(slip)

    lhs_op = tbem.boundary_operator(surface, surface, double_mthd, fmm_config, all_mesh)
    lhs = lhs_op.apply(soln)

    return local_and_global_errors(np.abs(lhs - rhs), surface)

def local_and_global_errors(dof_error, surface):
    facet_length = np.abs(surface.facets[:, 0, 0] - surface.facets[:, 1, 0])
    facet_error = dof_error[::2] + dof_error[1::2]
    facet_error *= facet_length
    l2_error = np.sqrt(np.mean(facet_error ** 2))
    return facet_error, l2_error

def plot_soln(surface, soln, show = True):
    xs, exact = get_exact(surface)
    import matplotlib.pyplot as plt
    minx = np.min(xs)
    maxx = np.max(xs)
    x = np.linspace(minx, maxx, 500)
    # plt.plot(x, np.arctan(1.0 / x) / np.pi, 'r*-')
    plt.plot(xs, soln, 'b*-')
    if show:
        plt.show()

def animate(all_steps):
    import matplotlib.pyplot as plt
    def f(step_idx):
        surface = all_steps[step_idx][0]
        soln = all_steps[step_idx][1]
        plt.cla()
        xs, exact = get_exact(surface)
        plt.plot(xs, soln, 'b*-')
        plt.ylim([-0.7, 0.7])
    fig = plt.figure(figsize = (5,4))
    from matplotlib import animation
    anim = animation.FuncAnimation(fig, f, frames = len(all_steps))
    anim.save('antiplane_refine.gif', writer = 'imagemagick', fps = 1)



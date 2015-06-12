from tbempy.ThreeD import *
import numpy as np

def unused_operator(obs_mesh, src_mesh):
    mthd = make_adaptive_integrator(1e-5, 5, 10, 3.0, LaplaceHypersingular())
    fmm_config = FMMConfig(0.3, 30, 250, 0.05, True)
    rhs_op = boundary_operator(obs_mesh, src_mesh, mthd, fmm_config, src_mesh)

# Just check that these don't produce an error
def test_empty_operator():
    unused_operator(Mesh(np.zeros((0, 3, 3))), Mesh(np.zeros((0, 3, 3))))

def test_empty_source_operator():
    unused_operator(
        Mesh(np.array([[[0, 0, 0], [1, 0, 0], [0, 1, 0]]])),
        Mesh(np.zeros((0, 3, 3)))
    )

if __name__ == '__main__':
    test_empty_operator()

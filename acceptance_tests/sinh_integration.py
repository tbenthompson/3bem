from tbempy.TwoD import make_sinh_integrator, ElasticHypersingular
from planestrain_fault import planestrain_with_integration_mthd

def test_sinh_integration():
    mthd = make_sinh_integrator(12, 3, 8, 3.0, ElasticHypersingular(30e9, 0.25))
    planestrain_with_integration_mthd(mthd)

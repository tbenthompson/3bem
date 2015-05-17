from tbempy.TwoD import make_sinh_integration_mthd
from planestrain_fault import planestrain_with_integration_mthd

def test_sinh_integration():
    def make_mthd(qs, K):
        return make_sinh_integration_mthd(12, qs, K)
    planestrain_with_integration_mthd(make_mthd)

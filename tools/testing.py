import numpy as np

def check_error_norm(est, exact, limit):
    error_norm = np.linalg.norm(exact - est, 2)
    exact_norm = np.linalg.norm(exact, 2)
    assert(error_norm / exact_norm < limit)

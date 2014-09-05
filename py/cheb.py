import numpy as np
import matplotlib.pyplot as plt

def cheb_polys(x, n_max, a = -1, b = 1):
    x_hat = 2 * ((x - a) / (b - a)) - 1
    try:
        length = x_hat.shape[0]
    except AttributeError:
        length = 1
    except IndexError:
        length = 1

    res = np.empty((length, n_max + 1))
    res[:, 0] = 1.0
    if n_max == 0:
        return res

    if len(x_hat.shape) > 1:
        res[:, 1:2] = x_hat
    else:
        res[:, 1] = x_hat
    if n_max == 1:
        return res
    for i in range(2, n_max + 1):
        res[:, i] = 2 * x_hat * res[:, i - 1] - res[:, i - 2]
    return res

def s_n(x, y, n, a = -1, b = 1):
    x_cheb = cheb_polys(x, n - 1, a, b)
    y_cheb = cheb_polys(y, n - 1, a, b)
    # chop the first value because it is always one
    x_cheb = x_cheb[:, 1:]
    y_cheb = y_cheb[:, 1:]
    return (1.0 / n) + (2.0 / n) * np.sum(x_cheb * y_cheb, axis = 1)

def cheb_pts(n, a = -1, b = 1):
    m = np.arange(1, n + 1)
    return a + (b - a) * 0.5 * (np.cos((2 * m - 1) * np.pi / (2 * n)) + 1)

def test_cheb_poly():
    theta = np.linspace(0.0, 1.0, 4)
    cos_theta = np.cos(theta)
    cheb_test = cheb_polys(cos_theta, 10)
    np.testing.assert_almost_equal(cheb_polys(0.5, 0), [[1]])
    exact = np.cos(np.outer(np.arange(11), theta).T)
    np.testing.assert_almost_equal(exact, cheb_test)

def test_cheb_interp():
    fnc = lambda x: np.sin(16 * x)
    n_pts_interp = 32

    n_locs = 300
    a = -1.0
    b = 1.0
    eval_locs = np.linspace(-1, 1, n_locs)
    exact = fnc(eval_locs)

    cheb_nodes = cheb_pts(n_pts_interp, a = a, b = b)

    est = np.zeros_like(exact)
    for i in range(n_locs):
        for m in range(n_pts_interp):
            est[i] += fnc(cheb_nodes[m]) * \
                      s_n(cheb_nodes[m], eval_locs[i], n_pts_interp, a, b)
    log_error = np.log(np.abs(exact - est) + 1e-16) / np.log(10)
    assert((log_error < 1e-6).all())
    plt.plot(eval_locs, exact, 'k-')
    plt.plot(eval_locs, est, 'b-')
    plt.plot(cheb_nodes, fnc(cheb_nodes), 'o')
    plt.show()

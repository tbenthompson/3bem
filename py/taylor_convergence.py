import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize

sequence = [0.219601, 0.133246, 0.0633529, 0.0372489, 0.0276542, 0.0203656, 0.0129682, 0.00630824, 0.0010842, -0.00253946, -0.00467146, -0.0055293, -0.00538493, -0.00453787, -0.00328608, -0.00189743, -0.000588182, 0.000488822, 0.00124668, 0.00165927, 0.00175052, 0.00158024, 0.00122883, 0.000782741, 0.000322262, -8.73964e-05, -0.000401329, -0.000596559, -0.000670554, -0.000637454, -0.000522897, -0.000358389, -0.000176002, -4.02464e-06, 0.000136066, 0.000231165, 0.000276707, 0.000275679, 0.000236815, 0.000172365, 9.58087e-05, 1.98161e-05, -4.53085e-05, -9.25897e-05, -0.000118726, -0.000123879, -0.000111042, -8.51304e-05, -5.19813e-05, -1.7401e-05, 1.36172e-05, 3.74303e-05, 5.20064e-05, 5.69245e-05, 5.31556e-05, 4.26929e-05, 2.81068e-05, 1.20969e-05, -2.89729e-06, -1.49855e-05, -2.29921e-05, -2.64967e-05, -2.57657e-05, -2.16017e-05, -1.51446e-05, -7.66056e-06, -3.47656e-07, 5.81634e-06, 1.01716e-05, 1.24143e-05, 1.25801e-05, 1.09839e-05, 8.13097e-06, 4.61746e-06, 1.03354e-06, -2.11622e-06, -4.46818e-06, -5.82676e-06, -6.16424e-06, -5.59756e-06, -4.34918e-06, -2.69967e-06, -9.40002e-07, 6.70059e-07, 1.93257e-06, 2.72878e-06, 3.02281e-06, 2.85353e-06, 2.31732e-06, 1.54582e-06, 6.82306e-07, -1.39856e-07, -8.13919e-07, -1.27025e-06, -1.48006e-06, -1.45297e-06, -1.2297e-06, -8.71582e-07, -4.48941e-07, -3.00326e-08, 3.28015e-07]
seq = np.array(sequence).astype(np.float64)
n = np.arange(seq.shape[0]).astype(np.float64)

def plot_oscillations():
    def f(x):
        A = x[0]
        B = x[1]
        omega = x[2]
        phi = x[3]
        start = 20
        return np.sum(((np.cos(n[start:] * omega - phi) ** 2 + (B * seq[start:] / (A ** n[start:])) ** 2) - 1) ** 2)



    res = minimize(f, [0.91, 1 / 0.012, 2 * np.pi / 16.5, 0.0])
    print res
    x = res.x

    plt.plot(n, np.abs(seq) / (0.91 ** n))
    plt.show()
    plt.plot(np.cos(n * x[2]), x[1] * seq / (x[0] ** n))
    # plt.plot(n, x[1] * seq / (x[0] ** n))
    # plt.plot(n, np.cos(n * x[2] + x[3]))
    plt.ylim([-1, 1])
    plt.show()

def more_plots():
    plt.plot(n, seq / (0.91 ** n))
    plt.show()
    log_mag = np.log(np.abs(seq)) / np.log(10)
    ratio = (log_mag[1:] - log_mag[:-1]) / (n[1:] - n[:-1])
    print np.mean(ratio)
    plt.plot(log_mag)
    plt.show()

def make_it_alternating():
    omega = 2 * np.pi / 16.5
    C = seq / np.cos(n * omega)
    print C
    # an = (-1) ** n * seq * np.exp(1j * n * (np.pi) * (1 - 2 / 16.5))
    # s = 0
    # for a in an:
    #     s += np.real(a)
    #     print s
    # print np.sum(seq)
    # print np.sum(an)

if __name__ == "__main__":
    # plot_oscillations()
    more_plots()
    # make_it_alternating()

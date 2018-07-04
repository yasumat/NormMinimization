import numpy as np
import matplotlib.pyplot as plt
from normmin import solve

'''
Example of robust ling fitting by sparse regression
'''


def true_fun(x):
    return 0.6 * x + 0.1


if __name__=='__main__':
    n_samples = 30    # number of data points
    n_outliers = 10
    # generate data points
    np.random.seed(0)
    x = np.sort(np.random.rand(n_samples))
    y = true_fun(x)
    # perturb y for simulating outliers
    indices = np.random.choice(n_samples, n_outliers)
    y[indices] = -y[indices] * np.random.rand(n_outliers)

    # matrix dimensions (m x n matrix)
    m = n_samples
    n = 2

    # Create a linear regression problem
    # y = a x + b, estimate a and b
    # y = A w, A = [x | 1], w = [a | b]^T
    A = np.ones((m, n))
    b = y
    for r in range(m):
        A[r, 0] = x[r]

    # (Case 1): Conventional least-squares fitting
    w1, residue1, ite1 = solve(A_list=[A], b_list=[b], lambda_list=[1], p_list=[2])
    # (Case 2): L1 sparse regression
    w2, residue2, ite2 = solve(A_list=[A], b_list=[b], lambda_list=[1], p_list=[1])

#    plt.style.use('fivethirtyeight')
    f1 = np.poly1d(w1)
    f2 = np.poly1d(w2)
    xp = np.linspace(0, 1, 100)
    plt.plot(x, true_fun(x), 'k-.', label='True function')
    plt.plot(xp, f1(xp), 'g-', label='L2')
    plt.plot(xp, f2(xp), 'r-', label='L1')

    plt.scatter(x, y, edgecolor='b', label='data')
    plt.legend()
    plt.xlabel('x')
    plt.ylabel('y')
    plt.xlim((0, 1))
    plt.ylim((-1, 1))
    plt.show()

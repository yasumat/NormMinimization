import numpy as np
import matplotlib.pyplot as plt
import normapprox as na

'''
Example of polynomial fitting with a regularizer using Norm Approximation
   y = w_1 * x^d + w_2 * x^(d-1) ... + w_d * x + w_{d+1}
'''


def true_fun(X):
    return np.cos(1.5 * np.pi * X)


if __name__=='__main__':
    degree = 15    # degree of polynomial
    n_samples = 30    # number of points

    # matrix dimensions (m x n matrix)
    m = n_samples
    n = degree + 1

    # generate data points
    np.random.seed(0)
    x = np.sort(np.random.rand(n_samples))
    y = true_fun(x) + np.random.randn(n_samples) * 0.15

    # Create a linear regression problem
    A = np.ones((m, n))
    b = y
    for r in range(m):
        for c in range(0, n-1):
            A[r, c] = np.power(x[r], n-1-c)

    # Case 1: Conventional least-squares fitting
    w1, residue1, ite1 = na.solve(A_list=[A], b_list=[b], lambda_list=[1], p_list=[2])
    # Case 2: Least-squares fitting with L2 regularization (special form of Ridge regression)
    w2, residue2, ite2 = na.solve(A_list=[A, np.identity(n)], b_list=[b, np.zeros(n)],
                                  lambda_list=[1, 1e-4], p_list=[2, 2])
    plt.style.use('fivethirtyeight')
    f1 = np.poly1d(w1)
    f2 = np.poly1d(w2)
    xp = np.linspace(0, 1, 100)
    plt.scatter(x, y, edgecolor='b', label='data')
    plt.plot(x, true_fun(x), 'k-.', label='True function')
    plt.plot(xp, f1(xp), 'g-', label='L2')
    plt.plot(xp, f2(xp), 'r-', label='L2 with regularization')
    plt.legend()
    plt.xlabel('x')
    plt.ylabel('y')
    plt.xlim((0, 1))
    plt.ylim((-2, 2))
    plt.show()


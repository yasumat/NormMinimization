import numpy as np
import matplotlib.pyplot as plt
import normapprox as na

'''
Example of polynomial fitting with different regularizers using Norm Approximation
   y = w_1 * x^d + w_2 * x^(d-1) ... + w_d * x + w_{d+1}
'''


def true_fun(X):
    return np.cos(1.5 * np.pi * X)


if __name__=='__main__':
    degree = 15    # degree of polynomial
    n_samples = 30    # number of points

    np.random.seed(0)
    x = np.sort(np.random.rand(n_samples))
    y = true_fun(x) + np.random.randn(n_samples) * 0.15

    A = np.ones((n_samples, degree+1))
    b = y
    for r in range(n_samples):
        for c in range(0, degree):
            A[r, c] = np.power(x[r], degree-c)

    # least-squares fitting
    w1, residue1, ite1 = na.solve(A_list=[A], b_list=[b], lambda_list=[1e-1], p_list=[2])
    # least-squares fitting with regularization
    w2, residue2, ite2 = na.solve(A_list=[A, np.identity(degree+1)], b_list=[b, np.zeros(degree+1)],
                                  lambda_list=[1000, 1], p_list=[2, 2])
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


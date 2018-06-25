import numpy as np
from normapprox import solve
import matplotlib.pyplot as plt


'''
Example of polynomial fitting with different regularizers using Norm Approximation
   y = w_0 + w_1 * x + w_2 * x^2 + ... + w_d * x^d
'''

if __name__=='__main__':
    d = 10    # degree of polynomial
    n = 5    # number of points

    np.random.seed(1)
    x = np.random.rand(n)
    y = np.random.rand(n)
    plt.plot(x, y, '*')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.show()

    A = np.ones((n, d+1))
    b = np.zeros(())
    for r in range(n):
        for c in range(1, d+1):
            A[r, c] = np.power(x[r], c)




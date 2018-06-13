import numpy as np
import warnings


def solve(A_list=None, b_list=None, lambda_list=None, p_list=None, max_iter=10000, tol=1.0e-8):
    """
    Solves a general norm approximation problem
        minimize \sum_i 1/p_i * \lambda_i * ||A_i x - b_i||_{p_i}^{p_i}
    :param A_list: List of matrices A_k \in \mathbb{R}^{m_k \times n}
    :param b_list: List of vectors b_k \in \mathbb{R}^{m_k}
    :param lambda_list: List of weighting scalar \lambda_k \in \mathbb{R}
    :param p_list: List of norm indicators p_k \ in \mathbb{R}
    :param max_iter: Maximum number of iterations
    :param tol: Tolerance
    :return: x \in \mathbb{R}^n
    """

    alpha = 0.0    # small value for regularizing the weighted least squares
    eps = 1.0e-8    # small value for avoiding zero-division in weight update

    if A_list is None or b_list is None or lambda_list is None or p_list is None:
        raise ValueError("Required argument is empty")
    K = len(A_list)
    if K != len(b_list) or K != len(lambda_list) or K != len(p_list):
        raise ValueError("Inconsistent arguments")
    if any(lambda_list[a] < 0 for a in range(len(lambda_list))):
        raise ValueError("lambda needs to be all positive")
    if any(p_list[a] < 0 for a in range(len(p_list))):
        raise ValueError("p needs to be all positive")
    n = A_list[0].shape[1]    # domain dim
    xold = np.ones((n, 1))
    In = np.identity(n)
    w_list = []
    e_list = []
    for k in range(K):
        w_list.append(np.identity(A_list[k].shape[0]))
        e_list.append(np.zeros((A_list[k].shape[0], 1)))

    print(e_list)
    ite = 0
    while ite < max_iter:
        ite = ite + 1
        A = alpha * In
        b = np.zeros((n, 1))
        # Create a normal equation of the problem
        for k in range(K):
            A = A + p_list[k] * lambda_list[k] * \
                np.dot(np.dot(np.dot(A_list[k].T, w_list[k].T), w_list[k]), A_list[k])
            b = b + p_list[k] * lambda_list[k] * \
                np.dot(np.dot(np.dot(A_list[k].T, w_list[k].T), w_list[k]), b_list[k])

        x = np.linalg.solve(A, b)
        for k in range(K):
            e_list[k] = b_list[k] - A_list[k].dot(x)

#        print('e_list:', e_list)
        if np.linalg.norm(x - xold) < tol:
            return x
        else:
            xold = x
        # update weights
        for k in range(K):
            print("k, E: ", k, 1.0 / np.asarray(np.maximum(np.power(np.fabs(e_list[k]), 1.0 - p_list[k] / 2.0), eps))[:, 0])
            w_list[k] = np.diag(
               np.asarray(1.0 / np.maximum(np.power(np.fabs(e_list[k]), 1.0 - p_list[k] / 2.0), eps))[:, 0])

    warnings.warn("Exceeded the maximum number of iterations")
    return x


def func01():
    np.random.seed(5)
    m = 20
    n = 10
    A = np.random.rand(m, n)
    x_gt = 3.0 * np.random.randn(n)
    inds = np.arange(n)
    np.random.shuffle(inds)
#    x_gt[inds[5:]] = 0
    b = np.dot(A, x_gt)
    b = np.reshape(b, (m, 1))
    A_list = [A, np.identity(n)]
    b_list = [b, np.zeros((n, 1))]
    p_list = [2, 1]
    lambda_list = [10.0, 1.0e-5]

    x1 = solve(A_list, b_list, p_list, lambda_list)

    print(np.reshape(x_gt, (n, 1)))
    print(x1)


def func02():
    np.random.seed(5)
    m = 20
    n = 10
    A = np.random.rand(m, n)
    x_gt = 3.0 * np.random.randn(n)
    b = np.dot(A, x_gt)

    x = np.linalg.lstsq(A, b)[0]

    print(x)
    print(x_gt)

def func03():
    print(np.power(0, 1))
    np.random.seed(5)
    m = 5
    n = 2
    A = np.random.rand(m, n)
    x_gt = 3.0 * np.random.randn(n)
    inds = np.arange(n)
    np.random.shuffle(inds)
    x_gt[inds[5:]] = 0
    b = np.dot(A, x_gt)
    b = np.reshape(b, (m, 1))
    A_list = [A]
    b_list = [b]
    p_list = [2]
    lambda_list = [10.0]

    x1 = solve(A_list, b_list, p_list, lambda_list)

    print(np.reshape(x_gt, (n, 1)))
    print(x1)

if __name__=='__main__':
    func01()

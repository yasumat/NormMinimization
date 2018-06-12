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

    alpha = 1.0e-8    # small value for regularizing the weighted least squares
    eps = 1.0e-8    # small value for avoiding zero-division in weight update

    if A_list is None or b_list is None or lambda_list is None or p_list is None:
        raise ValueError("Required argument is empty")

    K = len(A_list)
    if K != len(b_list) or K != len(lambda_list) or K != len(p_list):
        raise ValueError("Inconsistent arguments")
    if K != len(np.where(lambda_list > 0.0)[0]):
        raise ValueError("lambda needs to be all positive")
    if K != len(np.where(p_list > 0.0)[0]):
        raise ValueError("p needs to be all positive")
    n = A_list[0].shape[1]    # domain dim
    xold = np.zeros((n,1))
    In = np.identity(n)
    w_list = []
    e_list = []
    for k in range(K):
        w_list.append(np.identity(A_list[k].shape[0]))
        e_list.append(np.zeros((A_list[k].shape[0], 1)))

    ite = 0
    while ite < max_iter:
        ite = ite + 1
        A = alpha * In
        b = np.zeros((n, 1))
        for k in range(K):
            A = A + lambda_list[k] * np.dot(np.dot(A_list[k].T, w_list[k]), A_list[k]) * p_list[k]
            b = b + lambda_list[k] * np.dot(np.dot(A_list[k].T, w_list[k]), b_list[k]) * p_list[k]

        x = np.linalg.solve(A, b)
        for k in range(K):
            e_list[k] = b_list[k] - A_list[k].dot(x)

        if np.linalg.norm(x - xold) < tol:
            return x
        else:
            xold = x
        # update weights
        for k in range(K):
            w_list[k] = np.diag(np.asarray(1.0 / np.maximum(np.power(np.fabs(e_list[k]), 2.0 - p_list[k]), eps))[:, 0])
    warnings.warn("Exceeded the maximum number of iterations")


def test01():
    np.random.seed(5)
    m = 10
    n = 20
    A = np.random.rand(m, n)
    x_gt = 3.0 * np.random.randn(n)
    inds = np.arange(n)
    np.random.shuffle(inds)
    x_gt[inds[5:]] = 0
    b = np.dot(A, x_gt)
    b = np.reshape(b, (m, 1))
    print(b)
    print(A)


if __name__=='__main__':
    test01()

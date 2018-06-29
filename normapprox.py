import numpy as np
import warnings


def residue_norm(e_list=None):
    residue = 0.0
    for e in e_list:
        residue = residue + np.linalg.norm(e)
    return residue


def solve(A_list=None, b_list=None, lambda_list=None, p_list=None, max_iter=10000, tol=1.0e-8):
    """
    Solves a general norm approximation problem
        minimize_x \sum_k \lambda_k * ||A_k x - b_k||_{p_k}^{p_k}

    :param A_list: List of matrices A_k \in \mathbb{R}^{m_k \times n}
    :param b_list: List of vectors b_k \in \mathbb{R}^{m_k}
    :param lambda_list: List of weighting scalar \lambda_k \in \mathbb{R}
    :param p_list: List of norm indicators p_k \ in \mathbb{R}
    :param max_iter: Maximum number of iterations
    :param tol: Tolerance
    :return: x \in \mathbb{R}^n
    :return: residue norm
    :return: number of iterations spent for computation
    """

    alpha = 0    #1.0e-8   # small value for regularizing the weighted least squares
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
    for index in range(len(b_list)):
        if b_list[index].ndim == 1:
            b_list[index] = b_list[index][:, np.newaxis]

    n = A_list[0].shape[1]    # domain dim
    x_old = np.ones((n, 1))
    In = np.identity(n)
    w_list = []
    e_list = []
    for k in range(K):
        w_list.append(np.identity(A_list[k].shape[0]))
        e_list.append(np.zeros((A_list[k].shape[0], 1)))

    ite = 0
    while ite < max_iter:
        ite = ite + 1
        C = alpha * In    # n \times n matrix
        d = np.zeros((n, 1))    # n-vector
        # Create a normal equation of the problem
        for k in range(K):
            C = C + p_list[k] * lambda_list[k] * np.dot(np.dot(A_list[k].T, w_list[k]), A_list[k])
            d = d + p_list[k] * lambda_list[k] * np.dot(np.dot(A_list[k].T, w_list[k]), b_list[k])
        x = np.linalg.solve(C, d)
        for k in range(K):
            e_list[k] = b_list[k] - A_list[k].dot(x)
        # stopping criteria
        if np.linalg.norm(x - x_old) < tol:
            return x.ravel(), residue_norm(e_list), ite
        else:
            print(x[0], x[1])
            x_old = x
        # update weights
        for k in range(K):
            w_list[k] = np.diag(
               np.asarray(1.0 / np.maximum(np.power(np.fabs(e_list[k]), 2.0 - p_list[k]), eps))[:, 0])

    warnings.warn("Exceeded the maximum number of iterations")
    return x.ravel(), residue_norm(e_list), ite

import uncertainties.unumpy as unp
import numpy as np
from scipy.linalg import inv

def least_squares(A, measured, errors):
    size = A.shape[0]
    W = np.identity(size)
    row, col = np.diag_indices_from(W)
    W[row, col] = np.array(1/errors**2)
    # print(A, W)

    B = inv((A.T).dot(W).dot(A))
    C = (A.T).dot(W).dot(measured)

    mu = B.dot(C)
    mu_err = np.sqrt(np.diag(inv((A.T).dot(W).dot(A))))
    res = unp.uarray(mu, np.abs(mu_err))
    return res
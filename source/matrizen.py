import numpy as np


def laengs_mat(Z):
    A = np.array([[1, Z], [0, 1]])
    return A


def quer_mat(Z):
    A = np.array([[1, 0], [1/Z, 1]])
    return A
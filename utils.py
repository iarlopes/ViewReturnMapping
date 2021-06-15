"""
This file contains util functions used in the package ViewReturnMapping.

Igor A. Rodrigues Lopes, Jun 2021, Initial coding.
"""

import numpy as np


def additive_decomposition(Cauchy):
    """Perform additive decomposition of stress tensor."""
    # Get the pressure
    p = (Cauchy[0,0] + Cauchy[1,1] + Cauchy[2,2]) / 3.0
    # Get deviatoric stress
    s = Cauchy - p * np.eye(3)
    #
    return s, p


def symMatrix2vec(matrix):
    """Store a symmetric matrix in vector format."""
    return np.array([matrix[0,0],
                     matrix[1,1],
                     matrix[2,2],
                     matrix[0,1],
                     matrix[1,2],
                     matrix[0,2]])


def vec2symMatrix(vec):
    """Recover symmetric matrix from its vector format."""
    return np.array([[vec[0], vec[3], vec[5]],
                     [vec[3], vec[1], vec[4]],
                     [vec[5], vec[4], vec[2]]])

def coordsPiPlane(xyz):
    """Projects x,y,z coordinates to the pi-plane, according to the isometric perspective."""
    # Build transformation matrix to project to the isometric perspective
    Tmatrix = 1/np.sqrt(6) * np.array([[np.sqrt(3), 0         , -np.sqrt(3)],
                                       [-1        , 2         , -1]])
    return np.matmul(Tmatrix, xyz)
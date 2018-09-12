# -*- coding: utf-8 -*-


"""
Created on Tue Jun 12 2018
Last edit: Tue Jun 12 2018
@author: kfinc

"""

import numpy as np


def calculate_lsn_edges (A, labels):
    """Function calculates number of edges between and within predefined large-scale networks (LSNs).
    The function takes binary symetrical adjacency matrix, module assignment of each ROI and calculate number of edges between and within each
    large-scale network.

    Parameters
    ------------
    array: N x N binary ajdacency matrix
    array: N-length vector with module assignment for each node

    Returns
    ------------
    array: M x M matrix with number of edges between each module

    """
    columns = np.unique(labels)
    lsn_matrix = np.zeros((len(labels), len(columns)))
    lsn_edges = np.zeros((len(columns), len(columns)))

    for col in range(len(columns)):
        module = columns[col, ]
        for row in range(len(labels)):
            if (labels[row, ] == module):
                lsn_matrix[row, col] = 1
            else:
                lsn_matrix[row, col] = 0

    lsn_edges = lsn_matrix.T @ A @ lsn_matrix
    return lsn_edges


def allegiance_matrix(M):
    """Calculates M x M allegiance matrix from M x N matrix, where M represents nodes, and N represents time windows. 
    Each value on allegiance matrix represents the probability that node i and node j have been assigned to the same 
    functional community (Bassett et al., 2014).
    
    Parameters
    ------------
    array: M x N module assignment matrix 

    Returns
    ------------
    array: M x M allegiance matrix
        
    """
    
    roi_n = len(M[:,0])
    A = np.zeros((roi_n, roi_n))
    
    for i in range(roi_n):
        for j in range(i):
            if i == j:
                continue
            else:
                vector = M[i, :] == M[j, :]
                p = vector.mean()
                A[i, j] = p
    A = A + A.T
    return A


def allegiance_matrices_4d(M):
    """Calculates 4D array composed of allegiance matrices for each subject and each condition/session."
    
    Parameters
    ------------
    array: 4D array S(subjecys) x C(condition/session) x M(node) x N(window) 

    Returns
    ------------
    array: 4D array S(subjecys) x C(condition/session) x M(node) x M(node), where M x M is allegiance matrix 
        
    """
    
    AM = np.zeros((46, 4, 264, 264))
    sub_n = len(M[:,0,0,0])
    ses_n = len(M[0,:,0,0])
    
    for sub in range(sub_n):
        for ses in range(ses_n):
            AM[sub, ses, :, :] = allegiance_matrix(M[sub, ses, :, :])
    return AM    
    
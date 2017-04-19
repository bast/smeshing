"""
MATLAB compatibility methods

unique_rows    : Similar to unique(..., 'rows')
setdiff_rows   : Similar to setdiff(..., 'rows')
"""

__all__ = [
    'setdiff_rows',
    'unique_rows',
    ]

import numpy as np
import scipy.sparse as spsparse
import scipy.interpolate as spinterp

def setdiff_rows(A, B, return_index=False):
    """
    Similar to MATLAB's setdiff(A, B, 'rows'), this returns C, I
    where C are the row of A that are not in B and I satisfies
    C = A[I,:].

    Returns I if return_index is True.
    """
    A = np.require(A, requirements='C')
    B = np.require(B, requirements='C')

    assert A.ndim == 2, "array must be 2-dim'l"
    assert B.ndim == 2, "array must be 2-dim'l"
    assert A.shape[1] == B.shape[1], \
           "arrays must have the same number of columns"
    assert A.dtype == B.dtype, \
           "arrays must have the same data type"

    # NumPy provides setdiff1d, which operates only on one dimensional
    # arrays. To make the array one-dimensional, we interpret each row
    # as being a string of characters of the appropriate length.
    orig_dtype = A.dtype
    ncolumns = A.shape[1]
    dtype = np.dtype((np.character, orig_dtype.itemsize*ncolumns))
    C = np.setdiff1d(A.view(dtype), B.view(dtype)) \
        .view(A.dtype) \
        .reshape((-1, ncolumns), order='C')
    if return_index:
        raise NotImplementedError
    else:
        return C

def unique_rows(A, return_index=False, return_inverse=False):
    """
    Similar to MATLAB's unique(A, 'rows'), this returns B, I, J
    where B is the unique rows of A and I and J satisfy
    A = B[J,:] and B = A[I,:]

    Returns I if return_index is True
    Returns J if return_inverse is True
    """
    A = np.require(A, requirements='C')
    assert A.ndim == 2, "array must be 2-dim'l"

    orig_dtype = A.dtype
    ncolumns = A.shape[1]
    dtype = np.dtype((np.character, orig_dtype.itemsize*ncolumns))
    B, I, J = np.unique(A.view(dtype),
                        return_index=True,
                        return_inverse=True)

    B = B.view(orig_dtype).reshape((-1, ncolumns), order='C')

    # There must be a better way to do this:
    if (return_index):
        if (return_inverse):
            return B, I, J
        else:
            return B, I
    else:
        if (return_inverse):
            return B, J
        else:
            return B

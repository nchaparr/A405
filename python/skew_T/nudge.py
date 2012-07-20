import numpy as np

def nudge(Vec):
    """
    
    Returns an array with the same entires as Vec (1D array), except
    all entries that are equal (within a tolerance) are "nudged" (one of the
    entires is increased by a percentage).
    
    Tests
    - - - -
    >>> p = np.array([1.,1.,2.,3.,3.,4.])
    >>> pnew = nudge(p)
    >>> ptest = np.array([1, 1.001, 2, 3, 3.001, 4])
    >>> np.alltrue(abs(ptest - pnew)) < 1.e-8
    True
    
    """
    
    newVec = Vec
    hit, = np.where(np.abs(np.diff(Vec)) < 1.e-8)
    newVec[hit+1] = Vec[hit] + 1.e-3*Vec[hit]
    return newVec

def _test():
    import doctest
    doctest.testmod()

if __name__ == "__main__":
     _test()	
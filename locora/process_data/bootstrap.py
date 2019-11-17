import numpy as np

# http://nbviewer.jupyter.org/gist/aflaxman/6871948
def bootstrap(X, n=None):
    """ Bootstrap resample an array_like
    Parameters
    ----------
    X : array_like
      data to resample
    n : int, optional
      length of resampled array, equal to len(X) if n==None
    Results
    -------
    returns X_resamples
    """
    if n == None:
        n = len(X)
    elif n > len(X):
        n = len(X)
        
    resample_i = np.floor(np.random.rand(n)*len(X)).astype(int)
    X_resample = X[resample_i]
    return X_resample


def moving_block_bootstrap(X, n=None, b=None):
    """ Bootstrap resample an array_like
    Parameters
    ----------
    X : array_like
      data to resample
    n : int, optional
      total number of resampled arrays, equal to 10 if n==None
    b : int, optional
      length of resampled blocks, equal to len(X)/10 if b==None

    Results
    -------
    returns X_resamples with shape (n*b)
    """
    l = len(X)
    
    if type(n) == type(None):
        n = 10

    if type(b) == type(None):
        b = l/10

    resample_i      = np.random.randint(low=0, high=l-b+1, size=n)
    block_resample  = np.zeros((n,b), dtype=int)
    block_resample += resample_i.reshape((n,1))
    block_resample += np.arange(b, dtype=int)
    X_resample      = X[block_resample.reshape(b*n)]

    return X_resample

def block_bootstrap(X, n=None, b=None):
    
    l = len(X)
    
    if type(n) == type(None):
        n = 10

    if type(b) == type(None):
        b = l/10
        
    block_resample  = np.arange(l/b*b, dtype=int).reshape((l/b,b))
    resample_i      = np.random.randint(low=0, high=l/b, size=n)
    X_resample      = X[block_resample[resample_i].reshape(n*b)]

    return X_resample
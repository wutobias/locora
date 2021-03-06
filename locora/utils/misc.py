import numpy as np

def are_you_numpy(a):

    """
    Returns True if a is an instance of numpy.
    False otherwise.
    """

    return type(a).__module__ == np.__name__
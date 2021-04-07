import numpy as np

class A(object):
    __getitem__ = np.array
A = A()

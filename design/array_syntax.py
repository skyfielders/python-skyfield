from numpy import array
from time import time

class A(object):
    __getitem__ = array
A = A()

t0 = time()
x = array((1, 2, 3))
tt = time() - t0
print('Old way:', tt)

t0 = time()
x = A[1, 2, 3]
tt = time() - t0
print('A[] way:', tt)



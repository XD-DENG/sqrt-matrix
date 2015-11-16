import numpy
import scipy.linalg

a=numpy.array([[1, 1, 0, 0],
       [1, 0, 2, 0],
       [1, 0, 0, 3],
       [2, 2, 2, 2]])

x = scipy.linalg.sqrtm(a)

print type(x)
print x
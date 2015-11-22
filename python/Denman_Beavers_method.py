import numpy
from scipy.linalg import inv


def Denman_Beavers_sqrtm(A):
    Y = A
    Z = numpy.eye(len(A))

    error = 1
    error_tolerance = 1.5e-8

    flag = 1
    while(error > error_tolerance):
        Y_old = Y
        Y = (Y_old + inv(Z))/2
        Z = (Z + inv(Y_old))/2
        error_matrix = abs(Y - Y_old)
        error = 0
        # detect the maximum value in the error matrix
        for i in range(len(A)):
            temp_error = max(error_matrix[i])
            if(temp_error > error):
                error = temp_error


        flag = flag + 1


    print("Iteration Times: ", flag)
    return Y

a=numpy.array([[1, 1, 0, 0],
       [1, 0, 2, 0],
       [1, 0, 0, 3],
       [2, 2, 2, 2]])
print(Denman_Beavers_sqrtm(a))
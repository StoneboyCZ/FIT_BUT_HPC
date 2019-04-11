import numpy as np
import scipy as sc

def ode_linear(t,y,A,b):
    dy = np.zeros((2,1))
    #print(dy)
    print(A[0,0]*y[0]+A[0,1]*y[1])
    dy[0] = (A[0,0]*y[0]+A[0,1]*y[1])+b[0]
    dy[1] = (A[1,0]*y[0]+A[1,1]*y[1])+b[1]
    return dy

import numpy as np
import scipy.sparse as sparse

def specifyProblem():
    A = np.array(
        [
        [0,1],
        [-1,0]
        ]
    )


    b = np.array([[0,0]])

    y = np.array([[0,1]])
 
    return A,b,y
    


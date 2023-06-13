import numpy as np
from mpmath import mp

def integrate(A, b, y0, dt, tspan, eps, dps):
    '''
        Function perfoms numerical integration using the MTSM.
        
        Parameters:
        :param A: Jacobian matrix describing the system.
        :param b: RHS vector
        :param y0: initial conditions of the system
        :param dt float: Size of the integration step.
        :param tspan: Simulation time range. 
        :param dps: Number of valid decimal placesThe number of bits used for calculation. 
    '''    

    mp.dps = dps

    # number of equations
    ne = A.rows
    steps = round((tspan[1] - tspan[0])/dt)
   
    # number of steps used for the stopping rule
    nStopping = 3
    # result array
    y = mp.matrix(ne, steps+1)
    ORD = mp.matrix(1, steps+1)

    # copy the initial condition into the y
    y[:,0] = y0
    
    stopping = mp.matrix(1, nStopping)
    stopping[0,0:nStopping-1] = 1e10

    ORD[0,0] = 0

    for i in range(1,steps+1,1):
        dy = mp.matrix(ne, 1)
        y[:, i] = y[:, i-1]
        
        k = 0
        dy[:, k] = y[:, i]

        stopping[0,k] = mp.fabs(max(dy[:,k]))
        
        while mp.norm(stopping) > eps:
            Ay = A*dy[:,k]
            k = k+1  
            if k == 1:
                dy[:,k] = dt*(Ay + b)     
            else:
                dy[:,k] = mp.fdiv(dt,k)*(Ay)
            y[:,i] = y[:,i] + dy[:,k]
            stopping[0,k%nStopping] = mp.fabs(max(dy[:,k]))
            #print(stopping)
        #print(y[:,i])
        ORD[0,i] = k
    return y,ORD

    
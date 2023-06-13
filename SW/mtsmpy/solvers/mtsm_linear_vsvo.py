import numpy as np
import math

def integrate(A, b, y0, dt, tspan, eps):
    """
        Function performs numerical integration using the MTSM.

        Parameters:
        :param A: Jacobian matrix describing the system.
        :param b: RHS vector
        :param y0: initial conditions of the system
        :param dt: Size of the integration step.
        :param eps: Error per step.
        :param tspan: Simulation time range.
    """

    # number of equations
    ne, num_cols = A.shape

    lstol = 1e-10
    steps = round((tspan[1] - tspan[0])/dt)


    # number of steps used for the stopping rule
    nStopping = 3
    # result array
    #y = np.zeros(shape=(ne, steps+2))
    ORD = []

    y = y0

    # copy the initial condition into the y
    #y[:, [0]] = y0

    time = []
    time[0] = tspan[0]

    t = tspan[0]+dt
    i = 1
    ORD[0] = 0
    #for i in range(1, steps+1, 1):
    #while not math.isclose(t+dt, tspan[1], rel_tol=eps, abs_tol=eps) or (t + dt < tspan[1]):
    while time[i] + 1e-10 < tspan[1]:
        k = 0
        dy = np.zeros(shape=(ne, steps + 1))

        stopping = np.ones(nStopping) * 1e10

        y[:, [i]] = y[:, [i-1]]

        while np.linalg.norm(stopping) > eps:
            if k == 0:
                Ay = np.matmul(A, y[:, [i]])
                dy[:, [k]] = dt*(Ay + b)
            else:
                ADy = np.matmul(A, dy[:, [k-1]])
                dy[:, [k]] = (dt/(k+1))*(ADy)

            stopping[k%nStopping] = np.fabs(max(dy[:, k]))
            k = k + 1
        #np.concatenate((y, np.sum(dy[:, :], axis=1))) # y[:, [i]] = np.sum(dy[:, :])
        y[:, [i]] += np.sum(dy, axis=1, keepdims=True)

        ORD[i] = k
        time[i] = t
        i += 1
        t += dt

    return time, y, ORD

    
import numpy as np
from mpmath import mp
import matplotlib.pyplot as plt
import solvers.mtsm_linear_vpa as mtsm_linear_vpa
import solvers.mtsm_linear as mtsm_linear
import argparse
import scipy.integrate as integrate
import time


def linear(t, y):
    return np.matmul(A, y)+b


parser = argparse.ArgumentParser()
parser.add_argument("dps", help="number of decimal places", type=int)
args = parser.parse_args()

A = None
b = None
y0 = None
y = None
ORD = None

omega = 100

tspan = [0, 10]


if args.dps > 15:
    A = mp.matrix([[0, omega], [-omega, 0]])
    b = mp.matrix([[0], [0]])
    y0 = mp.matrix([0], [1])
    eps = mp.mpf(1e-9)
    dt = mp.mpf(1e-1)

    y, ORD = mtsm_linear_vpa.integrate(A, b, y0, dt, tspan, eps, args.dps)

    sol = integrate.solve_ivp(linear, tspan, y0, method='RK45', dense_output=True, rtol=1e-12, atol=1e-12)
else:
    A = np.array([[0, omega], [-omega, 0]])
    b = np.array([[0], [0]])
    y0 = np.array([[0], [1]])
    eps = 1e-9
    dt = 1e-1

    start = time.process_time()
    sol = integrate.solve_ivp(linear, tspan, y0[:, 0], method='RK45', dense_output=True, vectorized=True, rtol=1e-12, atol=1e-12)
    print(time.process_time() - start)

    start = time.process_time()
    t, y, ORD = mtsm_linear.integrate(A, b, y0, dt, tspan, eps)
    print(time.process_time() - start)


'''
plt.plot(sol.t, sol.y[0])
plt.show()

plt.plot(sol.y[0], sol.y[1], '*')
plt.grid()
plt.show()
'''

plt.plot(t, y[0, :])
plt.grid()
plt.show()

plt.plot(t, ORD ,'*')
plt.grid()
plt.show()

plt.plot(y[0, :], y[1, :], '*')
plt.grid()
plt.show()

 


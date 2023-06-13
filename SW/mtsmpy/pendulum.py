import numpy as np
from mpmath import mp
import matplotlib.pyplot as plt
import solvers.mtsm_linear_vpa as mtsm_linear_vpa
import solvers.mtsm_linear as mtsm_linear
import argparse

# from common.number import Number as number

parser = argparse.ArgumentParser()
parser.add_argument("dps", help="number of decimal places", type=int)
args = parser.parse_args()

A = None
b = None
y0 = None
y = None
ORD = None

''' 
n = number(3.56, args.dps)
print(n)
'''

mp.dps = args.dps

d = 1
m = 1
k = mp.fdiv(d, m)
g = mp.mpf(9.81)
L = 1
a = mp.mpf(m * g * L)

tspan = [0, 10]

A = mp.matrix([[-k, -a], [1, 0]])
b = mp.matrix([[1], [0]])
y0 = mp.matrix([5, 0])
eps = mp.mpf(1e-9)
dt = mp.mpf(1e-1)

y, ORD = mtsm_linear_vpa.integrate(A, b, y0, dt, tspan, eps, args.dps)

t = list(np.arange(tspan[0], tspan[1] + dt, dt, dtype=float))

plt.plot(t, y[0, :])
plt.grid()
plt.show()

plt.plot(t, ORD, '*')
plt.grid()
plt.show()

plt.plot(y[0, :], y[1, :], '*')
plt.grid()
plt.show()

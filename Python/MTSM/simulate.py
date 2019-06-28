import solver_linear as sl
import problem as p
import numpy as np
import scipy as sc

A,b,y = p.specifyProblem()

dt = 0.1
eps = 1e-9
tmax = 6

t,y = sl.solve_linear(A,b,y,dt,eps,tmax)
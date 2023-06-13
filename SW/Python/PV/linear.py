import numpy as np
import scipy as sc

def MTSM_linear(A,b,init,tmax,dt,eps,maxord):
    i = 0 # step index
    n = np.round(tmax/dt)

     


import numpy as np
import scipy as sc
import scipy.integrate as ode
import linear
import ode_linear

## 
settings = {}
settings['dt'] = 0.1
settings['eps'] = 1e-9
settings['tmax'] = 2*np.pi
settings['maxord'] = 64

# problem definition
A = np.matrix([
    [0,1], 
    [-1,0]
])

b = np.array([
    [0],
    [0]
])

init = np.array([0,1])

print(A[0,0])

#init = init.reshape((2, 1))
#print(init)
#print(init.ndim)

#####################
vy = ode.solve_ivp(fun=lambda t, y: ode_linear.ode_linear(t,y,A,b),t_span=[0,settings['tmax']],y0=init,t_eval=np.arange(0,settings['dt'],settings['tmax']))


#linear.MTSM_linear(A,b,init,settings['tmax'],settings['dt'],settings['eps'],settings['maxord'])




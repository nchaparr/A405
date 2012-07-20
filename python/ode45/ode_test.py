from scipy.integrate import ode
import numpy as np
import matplotlib.pyplot as plt

def F(t, x, a, b):
    xp = np.zeros((2,1))
    xp[0] = x[1]
    xp[1] = -t*x[0] - np.exp(t)*x[1] + a*np.sin(b*t)
    return xp

#use ode45 equivalent integrator
r = ode(F).set_integrator('dopri5')
r.set_f_params(3,2)
y0, t0 = [2, 8], 0
r.set_initial_value(y0, t0)

dt = 0.01
t1 = 10

t = np.array(t0)
y = np.array(y0)
                 
while r.successful() and r.t < t1:
    #returns y at the next time step, and sets y as the initial condition 
    r.integrate(r.t+dt)
    #keep track of the values of y at each time step
    y = np.vstack((y, r.y))
    t = np.vstack((t, r.t))
 
plt.figure(1)    
plt.plot(t, y[:,0])
plt.title('height vs. time, no stopping')
plt.show()
    

r.set_initial_value(y0, t0)
t = np.array(t0)
y = np.array(y0)


stopHeight = 5

#stops integration when height = 5
#try to figure out how to stop at height = 5 only if height is increasing/decreasing
while r.successful() and r.t < t1 and r.y[0] < stopHeight:
    r.integrate(r.t+dt)
    y = np.vstack((y, r.y))
    t = np.vstack((t, r.t))


plt.figure(2)    
plt.plot(t, y[:,0])
plt.title('height vs. time, stop when descending past 5')
plt.show()
    

    









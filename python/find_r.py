import numpy as np
import scipy.optimize as optimize
import matplotlib.pyplot as plt
import sys

""""Purpose: 

    To find the radius of a droplet nucleated by aerosol
    given a relative humidity value

    Steps:

    plot relative humidity as calculated by the given function,
    against radius

    determine the brackets, since the function has two values within
    a certain range

    root find over both bracketed intervals

"""

#Constants:
rho_w = 1*10**3 #density of water in Kg per m**3
M_w = 18*10**-3 #molecular mass of water in Kg/mol
k = 1.38*10**-23 #boltzmann's constant in m**2 Kg s**-2 K**-1
sigma_w = 7.28*10**-2 #surface tension for water at 25 C in N/m
Av_no = 6.0221*10**23 #Avagadro's number
Rv = 461.51

#Given:
r_a = 1*10**-8 #radius of aerosol in meters
v_a = (1.0*4/3)*np.pi*r_a**3 #volume of the aerosol partical in m**3
rho_a = 1.77*10**3 #density of hydrated sodium sulfate in Kg/m**3
m_a = v_a*rho_a #mass of the aerosol in Kg
M_a = 140*10**-3 #molecular mass of the aerosol in Kg/mol
I = 3 #distinct ions from aerosol
T = 290 #surrounding temperature in Kelvin

"""Constructing eq 6.8 from W&H"""

def rel_h_shift(rel_h_zero, r):

    v_d = (1.0*4/3)*np.pi*r**3 # volume of the spherical droplet
    rho_d = 1.0*(m_a + (v_d - v_a)*rho_w)/v_d #density of the droplet
    m_d = rho_d*v_d # mass of the droplet
    rho_dw = (v_d - v_a)*rho_w/v_d 
    
    n_w = 1.0*Av_no*(v_d - v_a)*rho_w/M_w #number concentration of  water molecules in droplet
    
    exp_fac = 1.0*2*sigma_w/(rho_d*Rv*T*r)
    fac1 = 1.0*(I*m_a*M_w)/(M_a*(m_d - m_a))

    zero = - rel_h_zero + np.exp(exp_fac)*(1.0/(1 + fac1))
    
    return [n_w, exp_fac, np.exp(exp_fac), (1.0/(1 + fac1))]


brackets = [1.2*10**-7]#list for brackets with first lower bracket
delta = [1000] #list for slopes with dummy first element
r_vals = np.linspace(1.2*10**-7, 10**-5, 1500)#array of radii
rel_h_vals  = np.zeros(len(r_vals))#empty array for relative humidity values to be calculated

for i in range(len(r_vals)):
    print rel_h_shift(0, r_vals[i])  
    rel_h_vals[i] = rel_h_shift(0, r_vals[i])[0]
    #get slopes at each interval
    delta.append(np.abs(1.0*(rel_h_vals[i] - rel_h_vals[i-1])/(r_vals[i]-r_vals[i-1])))
   
deltas = np.array(delta)

"""plot relative humidity against radius"""

fig=plt.figure(1)
fig.clf()
ax1=fig.add_subplot(111)
ax1.plot(r_vals, rel_h_vals)
fig.canvas.draw()
plt.show()

""" 
    Determine the brackets
"""

#get minimum slope, return index.  
est_max_rh = rel_h_vals.max
index = np.where(rel_h_vals == est_max_rh) 
print "index", index

#Is this less than tol?

#If so get max of the two corresponding relative humidity values and append this to the brackets list.   

if np.absolute(min_delta) < tol:
    if rel_h_val[index]>rel_h_val[index - 1]:
        bracket.append(r_val[index])
    else:
        np.append,bracket(r_val[index])

#If not
 
else:
    if min_delta > 0: # if the minimum slope is > than zero, shorten the interval over which it's determined.
        interval = [r_vals[index], 1.0*(r_vals[index] - r_vals[index-1])/2]  #new interval
        count = 0 
        max_iter = 1000
        while count< max_iter:            
            if np.abs(1.0*(rel_h_shift(interval[0]) - rel_h_shift(interval[1]))/(interval[0] - interval[1])) < tol:
                break            
            else:
                if rel_h_shift(interval[0]) > rel_h_shift(interval[1]): #indirectly check direction of slope, and narrow interval accordingly
                    interval[1] = 1.0*(interval[0] - interval[1])/2
                else:
                    interval[0] = 1.0*(interval[0] - interval[1])/2
            count = count+1

    #similarly if slope is negative, except bracket is shortened from the other end
    else:
         interval = [r_vals[index-1], 1.0*(r_vals[index] - r_vals[index-1])/2]   
         count = 0 
         max_iter = 1000
         while count< max_iter:
              if np.absolute(1.0*(rel_h_shift(interval[0]) - rel_h_shift(interval[1]))/(interval[0] - interval[1])) < tol:
                  break
            
              else:
                  if rel_h_shift(interval[0]) > rel_h_shift(interval[1]):
                      interval[1] = 1.0*(interval[0] - interval[1])/2
                  else:
                      interval[0] = 1.0*(interval[0] - interval[1])/2

              count = count+1
            
print "Maximum Iterations exceeded"
brackets.append(max(interval))
        
brackets.append(10**-5)

a = brackets[0]
b = brackets[1]
c = brackets[2]

                #try:
                #r1 = optimize.zeros.brenth(r_find, a, b, rel_h_zero)
                #except
   
                #try:
                #r2 = optimize.zeros.brenth(r_find, b, c, rel_h_zero)
                #except

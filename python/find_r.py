import numpy as np
import scipy.optimize as optimize
import matplotlib.pyplot as plt
import sys

""""
    Purpose: 

    To find the radius of a droplet nucleated by aerosol
    given a relative humidity value

    Steps:

    plot relative humidity as calculated by the given function: rel_h_shift()
    against radius

    set up brackets centered on critical radius which is find using: find_r_crit(), r_crit_zero()    

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
r_a = .5*10**-7 #radius of aerosol in meters
v_a = (1.0*4/3)*np.pi*r_a**3 #volume of the aerosol partical in m**3
rho_a = 1.77*10**3 #density of hydrated amonium sulfate in Kg/m**3
m_a = v_a*rho_a #mass of the aerosol in Kg
print 'mass of aerosol', m_a
M_a = 140*10**-3 #molecular mass of the aerosol in Kg/mol
I = 3 #distinct ions from aerosol
T = 290 #surrounding temperature in Kelvin

#radius range
[a, b] = [1.1*r_a, 10**-6]

"""Constructing eq 6.8 from W&H"""

def rel_h_shift(rel_h_zero, r):
    v_d = (1.0*4/3)*np.pi*r**3 # volume of the spherical droplet
    rho_d = 1.0*(m_a + (v_d - v_a)*rho_w)/v_d #density of the droplet
    m_d = rho_d*v_d # mass of the droplet
    rho_dw = (v_d - v_a)*rho_w/v_d #density of water in droplet    
    n_w = 1.0*Av_no*(v_d - v_a)*rho_w/M_w #number concentration of  water molecules in droplet    
    exp_fac = 1.0*2*sigma_w/(rho_d*Rv*T*r)
    fac1 = 1.0*(I*m_a*M_w)/(M_a*(m_d - m_a))
    zero = - rel_h_zero + np.exp(exp_fac)*(1.0/(1 + fac1))    
    return zero

def rcrit_zero(r):    
    v_d = (1.0*4/3)*np.pi*r**3 # volume of the spherical droplet
    rho_d = 1.0*(m_a + (v_d - v_a)*rho_w)/v_d #density of the droplet
    m_d = rho_d*v_d # mass of the droplet
    rho_dw = (v_d - v_a)*rho_w/v_d #density of water in droplet
    a = 1.0*2*sigma_w/(rho_d*Rv*T)
    b = 1.0*(I*m_a*M_w)/((1.0*4/3)*M_a*np.pi*rho_d)    
    rcrit_zero = np.sqrt(3.0*b/a) - r
    return rcrit_zero


def find_rcrit(a, b):
    r_vals = np.linspace(a, b, 1000)#create an array of radius values
    r_crit_rh = np.array([rel_h_shift(0, r) for r in r_vals])#calculate the relative humidities
    index = np.where(r_crit_rh == np.max(r_crit_rh)) #get the index of the maximum value
    a = r_vals[index]#get the corresponding radius, as an estimate to plug into the root finder
    r_crit = optimize.zeros.newton(rcrit_zero, a)#root find
    return r_crit

"""
   Calculate a range of relative humidity values over a range of radii
"""
r_vals = np.linspace(a, b, 1500)#array of radii
rel_h_vals  = np.zeros(len(r_vals))#empty array for relative humidity values to be calculated
for i in range(len(r_vals)):
    rel_h_vals[i] = rel_h_shift(0, r_vals[i])
        
"""plot relative humidity against radius"""

fig=plt.figure(1)
fig.clf()
ax1=fig.add_subplot(111)
ax1.plot(r_vals, rel_h_vals)
fig.canvas.draw()
plt.show()

"""
   Determine brackets centered around critical radius
"""
r_crit = find_rcrit(a, b)#determine critical radius
brackets = [a, r_crit, b] #list for bracket ends

a = brackets[0]
b = brackets[1]
c = brackets[2]

                #try:
                #r1 = optimize.zeros.brenth(r_find, a, b, rel_h_zero)
                #except
   
                #try:
                #r2 = optimize.zeros.brenth(r_find, b, c, rel_h_zero)
                #except

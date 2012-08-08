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

    root find over both bracketed intervals: do_r_find()
"""

#Constants:
rho_w = 1*10**3 #density of water in Kg per m**3
M_w = 18*10**-3 #molecular mass of water in Kg/mol
k = 1.38*10**-23 #boltzmann's constant in m**2 Kg s**-2 K**-1
sigma_w = 7.28*10**-2 #surface tension for water at 25 C in N/m
Av_no = 6.0221*10**23 #Avagadro's number
Rv = 461.51

#Given:
r_a = .2*10**-7 #radius of aerosol in meters
v_a = (1.0*4/3)*np.pi*r_a**3 #volume of the aerosol partical in m**3
rho_a = 1.77*10**3 #density of hydrated amonium sulfate in Kg/m**3
m_a = v_a*rho_a #mass of the aerosol in Kg
print 'mass of aerosol', m_a
M_a = 140*10**-3 #molecular mass of the aerosol in Kg/mol
I_no = 3 #distinct ions from aerosol
TempK = 290 #surrounding temperature in Kelvin

#radius range
[brackend_l, brackend_r] = [1.1*r_a, 10**-6]


def rel_h_shift(r, rel_h_zero):
    v_d = (1.0*4/3)*np.pi*r**3 # volume of the spherical droplet
    rho_d = 1.0*(m_a + (v_d - v_a)*rho_w)/v_d #density of the droplet
    a = 1.0*2*sigma_w/(rho_d*Rv*TempK)
    b = 1.0*(I_no*m_a*M_w)/((1.0*4/3)*np.pi*M_a*rho_d)
    zero = - rel_h_zero + 1 + 1.0*a/r  - 1.0*b/(r)**3#See excercise 6.12  W&H   
    return zero

def rcrit_zero(r):    
    v_d = (1.0*4/3)*np.pi*r**3 # volume of the spherical droplet
    rho_d = 1.0*(m_a + (v_d - v_a)*rho_w)/v_d #density of the droplet
    m_d = rho_d*v_d # mass of the droplet
    rho_dw = (v_d - v_a)*rho_w/v_d #density of water in droplet
    a = 1.0*2*sigma_w/(rho_d*Rv*TempK)
    b = 1.0*(I_no*m_a*M_w)/((1.0*4/3)*M_a*np.pi*rho_d)    
    rcrit_zero = np.sqrt(3.0*b/a) - r #See excercise 6.12 in W&H
    return rcrit_zero


def find_rcrit(a, b):
    r_vals = np.linspace(a, b, 1000)#create an array of radius values
    r_crit_rh = np.array([rel_h_shift(r, 0) for r in r_vals])#calculate the relative humidities
    index = np.where(r_crit_rh == np.max(r_crit_rh)) #get the index of the maximum value
    a = r_vals[index]#get the corresponding radius, as an estimate to plug into the root finder
    r_crit = optimize.zeros.newton(rcrit_zero, a)#root find
    return r_crit

def do_r_find(rh):
    r_crit = find_rcrit(brackend_l, brackend_r)[0]#determine critical radius
    print "Critical Radius", r_crit
    brackets = [brackend_l, r_crit, brackend_r] #list for bracket ends

    a = brackets[0]
    b = brackets[1]
    c = brackets[2]
    r1 = optimize.zeros.brenth(rel_h_shift, a, b, rh)
    print "root found radius in first bracket", r1

    try:
        r2 = optimize.zeros.brenth(rel_h_shift, b, c, rh)
    except ValueError:
        r2 = 0
        print "only one radius for this relative humidity"
    else:        
        print "root found radius in second bracket", r2
        
    return [r1, r2]
    
if __name__ == "__main__":

    r_vals = np.linspace(brackend_l, brackend_r, 1500)#array of radii
    rel_h_vals  = np.zeros(len(r_vals))#empty array for relative humidity values to be calculated
    for i in range(len(r_vals)):
        rel_h_vals[i] = rel_h_shift(r_vals[i], 0)
        

    fig=plt.figure(1)
    fig.clf()
    ax1=fig.add_subplot(111)
    ax1.plot(r_vals, rel_h_vals)
    plt.ylim(.9, 1.01)
    fig.canvas.draw()
    plt.show()

    rh = float(raw_input('relative humidity'))
    radii = do_r_find(rh)
    



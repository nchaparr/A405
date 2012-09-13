from netCDF4 import Dataset
import site
import sys
site.addsitedir('/nfs/kite/home/nchaparr/A405/repos/a405repo/python//thermlib')
site.addsitedir('/nfs/kite/home/nchaparr/A405/repos/a405repo/python//skew_T')
from scipy.integrate import ode
import matplotlib.pyplot as plt
import numpy as np

""" 
Pulls data from littlerock sounding.

Sets initial values of for parcel variables.

Integrator integrates function (F) for calculating the time derivatives for all variables

"""

from constants import constants as c
from nudge import nudge
from my_thetaep import thetaep
from calcVars import *
from T_thetaep import t_thetaep
from wsat import wsat
from findWvWl import findWvWl
from find_r import *
from aerodist import *
from esat import esat
    
def ode_littlerock():
    filename = 'littlerock.nc'
    print 'reading file: %s\n' %(filename)
    nc_file = Dataset(filename)
    var_names = nc_file.variables.keys()
    #print nc_file.ncattrs()
    #print nc_file.units
    #print nc_file.col_names
    
    sound_var = nc_file.variables[var_names[3]]
    press = sound_var[:,0]
    height = sound_var[:,1]
    temp = sound_var[:,2]
    dewpoint = sound_var[:,3]
    RelH = sound_var[:,4]
    
    #height must have unique values
    newHeight= nudge(height)
    #Tenv and TdEnv interpolators return temp. in deg C, given height in m
    #Press interpolator returns pressure in hPa given height in m
    interpTenv = lambda zVals: np.interp(zVals, newHeight, temp)
    interpTdEnv = lambda zVals: np.interp(zVals, newHeight, dewpoint)
    interpPress = lambda zVals: np.interp(zVals, newHeight, press)
    interpRelH = lambda zVals: np.interp(zVals, newHeight, RelH)
    p880_level = np.where(abs(880 - press) < 2.)
    p900_level = np.where(abs(900 - press) < 2.)

    height0 = 885    

    #Setting initial properties of parcel
    press0 = interpPress(height0)*100
    Tparc0 = 290
    Wt = .013
    ws0 = wsat(Tparc0, press0)
    es0 = esat(Tparc0)
    pressv0 = (1.0*Wt/(Wt + c.eps))*press0
    SS0 = 1.0*pressv0/es0 - 1
    RelH0 = (100*Wt/ws0)

    #if saturated, stop
    
    if ws0 <= Wt:
        sys.exit('Saturated.  Change Wt.')       
    else:
        wv0 = Wt
        
    thetaeVal = thetaep(wv0, Tparc0, press0)
        
    #aerosol properties
     
    r0, r_a, Num_a = drop_props(RelH0)
    
    Num_rad = len(r0)
    print''
    print'Here are my initial Values'
    print''
    yinit = [height0, 0.5, Tparc0, SS0]
    for i in range(len(r0)):
        yinit.append(r0[i])
    print''
    print yinit
    print''
        
    tinit = 0
    tfin = 200
    dt = 1
   
    r = ode(F).set_integrator('dopri5')

    r.set_f_params(Wt, interpTenv, interpTdEnv, interpPress)

    print''
    print'r is about to set initial_value'
    print''

    #r.set_f_params(Wt, rho_a, r_a, r0, Num_a, interpTenv, interpTdEnv, interpPress)
    r.set_initial_value(yinit, tinit)

    print ''
    print'r completed setting initial value'
    print''
    
    y = np.array(yinit)
    t = np.array(tinit)
    
    Press = [press0]
    WSat = [ws0]
    #stop integration when the parcel changes direction, or time runs out

    print''
    print'about to enter while loop'
    print''
    
    while r.successful() and r.t < tfin and r.y[1] > 0:
        
        print''
        print'r.successful(), r.t, r.y',r.successful() ,r.t, r.y
        print''
        #find y at the next time step
        #(r.integrate(t) updates the fields r.y and r.t so that r.y = F(t) and r.t = t 
        #where F is the function being integrated)
        
        r.integrate(r.t+dt)

        print''
        print'have successfully integrated'
        print''
        
        #keep track of y at each time step
        y = np.vstack((y, r.y))
        t = np.vstack((t, r.t))
        P = interpPress(r.y[0])*100
        Press.append(P)
        
        #test for thetaep and point at which parcel becomes saturated
        
        T = r.y[2]
        [wv, wl] = findWvWl(T, Wt, P)
        
        ws = wsat(T, P)
        WSat.append(ws)
       
        #if Wt > ws - .000001 and Wt < ws + .000001:

        #print 'becomes saturated at around:' , r.y[0], 'meters'
               
    wvel = y[:,1]
    Tparc = y[:,2]-273.5
    height = y[:,0]
    SS = y[:,3]
    Press = np.array(Press)
    #radius = y[:,4] # will become a 2d array of 
    radii = y[:, 4:4+len(r0)]
        
    fig1 = plt.figure(1)
    plt.clf()
    plt.ylabel('height above surface (m)')
    ax1=fig1.add_subplot(111)
    plt.plot(WSat, -Press, 'ro')
    labels = ax1.get_xticklabels()
    for label in labels:
        label.set_rotation(30)
    plt.xlabel('Wv')
    plt.xlim(0, .017)

    fig2 = plt.figure(2)
    #will become another figure with a loop, plotting radii
    ax2=fig2.add_subplot(111)
    #plt.plot(radius, -Press, 'o')
    for i in len(r0): 
        plt.plot(radius[:,i], -Press)
    labels = ax2.get_xticklabels()
    for label in labels:
       label.set_rotation(30) 
       #plt.legend(loc = 'lower left')
    plt.xlabel('Radii')
    plt.xlim(.00000001, .00001)
    plt.show()

    fig3 = plt.figure(3)
    plt.clf()
    #plt.ylabel('height above surface (m)')
    #ax3=fig3.add_subplot(121)
    #plt.plot(wvel, -Press, 'k *')
    #labels = ax3.get_xticklabels()
    #for label in labels:
    # label.set_rotation(30)
    #plt.xlabel('Vertical Velocity')
    
    ax3=fig3.add_subplot(111)
    plt.plot(Tparc, -Press, 'o')
    #labels = ax4.get_xticklabels()
    #for label in labels:
    #   label.set_rotation(30) 
    plt.legend(loc = 'lower left')
    plt.xlabel('Temperature')
    plt.show()

def drop_props(RelH0):
    print''
    print'now running drop_props'
    print''
    r_a, Num_a, mass_a, prob_a = Aero_dist()
    rho_a = 1775

    print''
    print'Relative Humitity', RelH0, 1.0*RelH0/100
    print''
    
    #r0 = np.array([do_r_find(1.0*RelH0/100, r, rho_a)[0] for r in r_a])
    
    r0 = np.zeros(len(r_a))
    for i in range(len(r_a)):
        print''
        print 'index, dry radius',i, r_a[i]
        print''
        r0[i] = do_r_find(1.0*RelH0/100, r_a[i], rho_a)[0]

    print''
    print'now leaving drop_props'
    print''

    return r0, r_a, Num_a 
   
#F returns the buoyancy (and height) and rates of change of Temperature, Droplet Radius and Vapour Pressure with time, at a given time step and height
#def F(t, y, Wt, rho_a, r_a, r0, interpTenv, interpTdEnv, interpPress):
def F(t, y, Wt, interpTenv, interpTdEnv, interpPress):
    [rho_a, RelH0, Num_rad] = [1775, 90, 24] 
    print''
    print'now starting F'
    print''
    #yp = np.zeros((5,1))#will be bigger to accomodate radii
    yp = []
    #yp[0] = y[1]
    yp.append(y[1])
    #whats a good way of including all the radii here? pass radii as an array eg y[5:25], receive back as an array within an array from calc_Vars, Vars = calc_Vars(), and unpack to fill yp 
    #will also be passing Num_a
    #yp[1], yp[2], yp[3], yp[4] = calc_Vars(y[0], Wt, y[2], y[1], y[4], y[3], rho_a, r_a, r0, 140*10**-3, 3, interpTenv, interpTdEnv, interpPress)
    Vars = calc_Vars(y[0], Wt, y[2], y[1], y[4:4 + Num_rad], y[3], rho_a, RelH0, 140*10**-3, 3, interpTenv, interpTdEnv, interpPress)
    yp.extend([Vars[0], Vars[1], Vars[2]])

    for i in range(len(Vars[3])): 
        yp.append(Vars[3][i]) 
    print''
    print'F is finished now.'
    print''
    return yp

if __name__ == "__main__":
    ode_littlerock()








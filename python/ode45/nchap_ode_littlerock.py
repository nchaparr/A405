from netCDF4 import Dataset
import site
import sys
site.addsitedir('/nfs/kite/home/nchaparr/A405/repos/a405repo/python//thermlib')
site.addsitedir('/nfs/kite/home/nchaparr/A405/repos/a405repo/python//skew_T')
from scipy.integrate import ode
import matplotlib.pyplot as plt
import numpy as np

""" 

"""

from constants import constants as c
from nudge import nudge
from my_thetaep import thetaep
from calcVars import calcBuoy
from T_thetaep import t_thetaep
from wsat import wsat
from calcVars import calcdT
from calcVars import calcdr
from findWvWl import findWvWl
from find_r import *
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
    pressv0 = 1.0*(ws0/(ws0 + c.eps))*press0
    RelH0 = (100*Wt/ws0)

    #if saturated, stop
    
    if ws0 <= Wt:
        sys.exit('Saturated.  Change Wt.')       
    else:
        wv0 = Wt
        
    thetaeVal = thetaep(wv0, Tparc0, press0)
        
    print 'Initial Temperature, thetaeVal', Tparc0, thetaeVal
    print 'Initial Relative humidity, Wsat, Wt', RelH0, ws0, Wt
      
    #get initial radius
    r0 = do_r_find(1.0*RelH0/100)[0]
    print 'Initial radius', r0       

    yinit = [height0, 0.5, Tparc0, Tparc0, Tparc0, r0, pressv0]  #(intial velocity = 0.5 m/s, initial height in m)
    tinit = 0
    tfin = 100
    dt = .1
    
    #want to integrate F using ode45 (from MATLAB) equivalent integrator
    r = ode(F).set_integrator('dopri5')
    r.set_f_params(Wt, interpTenv, interpTdEnv, interpPress)
    r.set_initial_value(yinit, tinit)
    
    y = np.array(yinit)
    t = np.array(tinit)
    
    #stop integration when the parcel changes direction, or time runs out
    while r.successful() and r.t < tfin and r.y[1] > 0:

        #find y at the next time step
        #(r.integrate(t) updates the fields r.y and r.t so that r.y = F(t) and r.t = t 
        #where F is the function being integrated)
        r.integrate(r.t+dt)

        #keep track of y at each time step
        y = np.vstack((y, r.y))
        t = np.vstack((t, r.t))
        P = interpPress(r.y[0])*100

        #test for thetaep and point at which parcel becomes saturated
        
        T = r.y[2]
        [wv, wl] = findWvWl(T, Wt, P)
        ws = wsat(T, P)
        
       
        if Wt > ws - .000001 and Wt < ws + .000001:
            print 'becomes saturated at around:' , r.y[0], 'meters'
                            
        print "thetaep test: ", thetaep(wv, T, P), thetaeVal
    
    wvel = y[:,1]
    Tparc = y[:,2]
    height = y[:,0]
    Tcheck = y[:,3]
    Tcheck1 = y[:,4]
    pressv = y[:,6]
    radius = y[:,5]
    
    fig = plt.figure(1)
    #ax1=fig.add_subplot(121)
    #ax1.plot(wvel, height, 'o')
    #plt.xlabel('vertical velocity')
    #plt.ylabel('height above surface (m)')
    ax2=fig.add_subplot(121)
    plt.plot(pressv, height, 'k *')
    #plt.plot(Tcheck, height, 'b o')
    #plt.plot(Tcheck1, height, 'r +')
    labels = ax2.get_xticklabels()
    for label in labels:
        label.set_rotation(30)
    ax3=fig.add_subplot(122)
    plt.plot(radius, height, 'o')
    labels = ax3.get_xticklabels()
    for label in labels:
        label.set_rotation(30) 
    #plt.legend(loc = 'lower left')
    plt.xlabel('radius profile')
    plt.show()
        
#F returns the buoyancy (and height)at a given time step and height
def F(t, y, Wt, interpTenv, interpTdEnv, interpPress):
    yp = np.zeros((7,1))
    yp[0] = y[1]
    yp[1] = calcBuoy(y[0], Wt, y[2], interpTenv, interpTdEnv, interpPress)
    yp[2], yp[3], yp[4] = calcdT(y[0], Wt, y[2], interpPress, y[1])
    yp[5], yp[6] = calcdr(y[5], y[2], y[0], interpPress, y[6], 1.77*10**3, .2*10**-7, 140*10**-3, 3, y[1])
    return yp

if __name__ == "__main__":
    ode_littlerock()








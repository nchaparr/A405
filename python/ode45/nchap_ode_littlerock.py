from netCDF4 import Dataset
import site
site.addsitedir('/nfs/kite/home/nchaparr/A405/repos/a405repo/python//thermlib')
site.addsitedir('/nfs/kite/home/nchaparr/A405/repos/a405repo/python//skew_T')
site.addsitedir('/nfs/kite/home/nchaparr/A405/repos/a405repo/python//ode45')
from scipy.integrate import ode
import matplotlib.pyplot as plt
import numpy as np

""" """

from constants import constants as c
from nudge import nudge
from new_thermo import thetaep
from calcVars import calcBuoy
from findTmoist import findTmoist
from wsat import wsat
from calcVars import calcdT 
from findWvWl import findWvWl
    
def ode_littlerock():
    filename = 'littlerock.nc'
    print 'reading file: %s\n' %(filename)
    nc_file = Dataset(filename)
    var_names = nc_file.variables.keys()
    print nc_file.ncattrs()
    print nc_file.units
    print nc_file.col_names
    
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
    p900_level = np.where(abs(900 - press) < 2.)
    p800_level = np.where(abs(800 - press) < 7.)

    thetaeVal=thetaep(dewpoint[p900_level] + c.Tc,temp[p900_level] + c.Tc, press[p900_level]*100.)
    
    height_800=height[p800_level]
    height_900=height[p900_level]

    press0 = interpPress(height_800)*100
    Tcloud0 = findTmoist(thetaeVal, press0)[0]
    RelH0 = interpRelH(height_900)[0]

    Ws0 = wsat(Tcloud0, press0)[0]
    Wt0 = 1.0*RelH0*Ws0/100

    print 'Initial Relative humidity, Wsat and Wt of parcel', RelH0, Ws0, Wt0   #get initial radius
    

    yinit = [height_800, 0.5, Tcloud0, Tcloud0, Tcloud0]  #(intial velocity = 0.5 m/s, initial height in m)
    tinit = 0
    tfin = 2500
    dt = 10
    
    #want to integrate F using ode45 (from MATLAB) equivalent integrator
    r = ode(F).set_integrator('dopri5')
    r.set_f_params(Wt0, thetaeVal, interpTenv, interpTdEnv, interpPress)
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

        T = r.y[2]
        wv = findWvWl(T, Wt0, P)[0]
        e = 1.0*wv*P/(c.eps + wv);
        denom=(1.0*17.67/np.log(e/611.2)) - 1.;
        Td = 1.0*243.5/denom;
        Td = Td + 273.15;
                
        print "thetae test: ", thetaep(Td, T, P), thetaeVal[0]
        
    wvel = y[:,1]
    Tcloud = y[:,2]
    height = y[:,0]
    Tcheck = y[:,3]
    Tcheck1 = y[:,4]
        
    fig = plt.figure(1)
    ax1=fig.add_subplot(121)
    ax1.plot(wvel, height, 'o')
    plt.xlabel('vertical velocity')
    plt.ylabel('height above surface (m)')
    ax2=fig.add_subplot(122)
    plt.plot(Tcloud, height, 'o', label = '1')
    plt.plot(Tcheck, height, '+', label = '2')
    #plt.plot(Tcheck1, height, 'o', label = '3')
    plt.legend(loc = 'lower left')
    plt.xlabel('temperature profiles')
    plt.show()
        
#F returns the buoyancy (and height)at a given time step and height
def F(t, y, Wt0, thetae0, interpTenv, interpTdEnv, interpPress):
    yp = np.zeros((5,1))
    yp[0] = y[1]
    yp[1] = calcBuoy(y[0], Wt0, thetae0, interpTenv, interpTdEnv, interpPress)
    yp[2], yp[3], yp[4] = calcdT(y[0], Wt0, interpPress, thetae0, y[1])
    return yp

if __name__ == "__main__":
    ode_littlerock()








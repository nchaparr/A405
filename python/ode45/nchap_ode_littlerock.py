from netCDF4 import Dataset
import site
site.addsitedir('/nfs/kite/home/nchaparr/A405/repos/a405repo/python//thermlib')
site.addsitedir('/nfs/kite/home/nchaparr/A405/repos/a405repo/python//skew_T')
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
from calcVars import calcdr
from findWvWl import findWvWl
from find_r import * 
    
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
    p880_level = np.where(abs(880 - press) < 2.)
    p900_level = np.where(abs(900 - press) < 2.)

    thetaeVal=thetaep(dewpoint[p880_level] + c.Tc,temp[p880_level] + c.Tc, press[p880_level]*100.)
    
    height_900=height[p900_level]
    height_880=height[p880_level]
    print 'initial height',height_900

    press0 = interpPress(height_900)*100
    Tparc0 = findTmoist(thetaeVal, press0)[0]
    RelH0 = interpRelH(height_880)[0]

    Ws0 = wsat(Tparc0, press0)[0]
    Wt0 = 1.0*RelH0*Ws0/100

    #get initial radius
    r0 = do_r_find(RelH0/100)[0]
    
    print 'Initial Relative humidity, Wsat and Wt of parcel, initial radius', RelH0, Ws0, Wt0, r0       

    yinit = [height_900, 0.5, Tparc0, Tparc0, Tparc0, r0]  #(intial velocity = 0.5 m/s, initial height in m)
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
    Tparc = y[:,2]
    height = y[:,0]
    Tcheck = y[:,3]
    Tcheck1 = y[:,4]
    radius = y[:,5]
    
    fig = plt.figure(1)
    #ax1=fig.add_subplot(121)
    #ax1.plot(wvel, height, 'o')
    #plt.xlabel('vertical velocity')
    #plt.ylabel('height above surface (m)')
    #ax2=fig.add_subplot(121)
    #plt.plot(Tparc, height, 'o', label = '1')
    ax3=fig.add_subplot(111)
    plt.plot(radius, height, 'o')
    labels = ax3.get_xticklabels()
    for label in labels:
        label.set_rotation(30) 
    #plt.legend(loc = 'lower left')
    plt.xlabel('radius profile')
    plt.show()
        
#F returns the buoyancy (and height)at a given time step and height
def F(t, y, Wt0, thetae0, interpTenv, interpTdEnv, interpPress):
    yp = np.zeros((6,1))
    yp[0] = y[1]
    yp[1] = calcBuoy(y[0], Wt0, thetae0, interpTenv, interpTdEnv, interpPress)
    yp[2], yp[3], yp[4] = calcdT(y[0], Wt0, interpPress, thetae0, y[1])
    yp[5] = calcdr(y[5], y[2], y[0], interpPress, 1.77*10**3, .2*10**-7, 140*10**-3, 3)
    return yp

if __name__ == "__main__":
    ode_littlerock()








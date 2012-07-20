from netCDF4 import Dataset
import site
site.addsitedir('C:\Users\Den\mya405\python\\thermlib')
site.addsitedir('C:\Users\Den\mya405\python\\skew_T')
from scipy.integrate import ode
import matplotlib.pyplot as plt
import numpy as np

from constants import constants as c
from nudge import nudge
from new_thermo import thetaep
from calcBuoy import calcBuoy
from findTmoist import findTmoist

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
    
    #height must have unique values
    newHeight= nudge(height)
    #Tenv and TdEnv interpolators return temp. in deg C, given height in m
    #Press interpolator returns pressure in hPa given height in m
    interpTenv = lambda zVals: np.interp(zVals, newHeight, temp)
    interpTdEnv = lambda zVals: np.interp(zVals, newHeight, dewpoint)
    interpPress = lambda zVals: np.interp(zVals, newHeight, press)
    p900_level = np.where(abs(900 - press) < 2.)
    p800_level = np.where(abs(800 - press) < 7.)
    thetaeVal=thetaep(dewpoint[p900_level] + c.Tc,temp[p900_level] + c.Tc,press[p900_level]*100.)
    
    height_800=height[p800_level]
    
    yinit = [0.5, height_800]  #(intial velocity = 0.5 m/s, initial height in m)
    tinit = 0
    tfin = 2500
    dt = 10
    
    #want to integrate F using ode45 (from MATLAB) equivalent integrator
    r = ode(F).set_integrator('dopri5')
    r.set_f_params(thetaeVal, interpTenv, interpTdEnv, interpPress)
    r.set_initial_value(yinit, tinit)
    
    y = np.array(yinit)
    t = np.array(tinit)
    
    #stop integration when the parcel changes direction, or time runs out
    while r.successful() and r.t < tfin and r.y[0] > 0:
        #find y at the next time step
        #(r.integrate(t) updates the fields r.y and r.t so that r.y = F(t) and r.t = t 
        #where F is the function being integrated)
        r.integrate(r.t+dt)
        #keep track of y at each time step
        y = np.vstack((y, r.y))
        t = np.vstack((t, r.t))
        
    wvel = y[:,0]
    height = y[:,1]
    
    plt.figure(1)
    plt.plot(wvel, height)
    plt.xlabel('vertical velocity (m/s)')
    plt.ylabel('height about surface (m)')
    plt.show()
        
#F returns the buoyancy (and height)at a given time step and height
def F(t, y, thetae0, interpTenv, interpTdEnv, interpPress):
    #y[0] is the velocity, y[1] is the height
    yp = np.zeros((2,1))
    #yp[0] is the buoyancy (acceleration), yp[1] is the velocity
    yp[0] = calcBuoy(y[1], thetae0, interpTenv, interpTdEnv, interpPress)
    yp[1] = y[0] 
    return yp

if __name__ == "__main__":
    ode_littlerock()








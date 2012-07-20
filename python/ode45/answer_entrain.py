from netCDF4 import Dataset
import site
site.addsitedir('C:\Users\Den\mya405\python\\thermlib')
site.addsitedir('C:\Users\Den\mya405\python\\skew_T')
from scipy.integrate import ode
import matplotlib.pyplot as plt
import numpy as np
from constants import constants as c
from nudge import nudge
from new_thermo import thetaep, tinvert_thetae, wsat
from calcBuoy import calcBuoy
from findTmoist import findTmoist

def answer_entrain():
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
    envHeight= nudge(height)
    #Tenv and TdEnv interpolators return temp. in deg C, given height in m
    #Press interpolator returns pressure in hPa given height in m
    interpTenv = lambda zVals: np.interp(zVals, envHeight, temp)
    interpTdEnv = lambda zVals: np.interp(zVals, envHeight, dewpoint)
    interpPress = lambda zVals: np.interp(zVals, envHeight, press)
    
    p900_level = np.where(abs(900 - press) < 2.)
    p800_level = np.where(abs(800 - press) < 7.)
    thetaeVal=thetaep(dewpoint[p900_level] + c.Tc,temp[p900_level] + c.Tc,press[p900_level]*100.)
    height_800=height[p800_level]
    wTcloud = wsat(dewpoint[p900_level] + c.Tc, press[p900_level]*100.)
    entrain_rate = 2.e-4
    winit = 0.5 #initial velocity (m/s)
    yinit = [winit, height_800, thetaeVal, wTcloud]  
    tinit = 0
    tfin = 2500
    dt = 10
    
    #want to integrate F using ode45 (from MATLAB) equivalent integrator
    r = ode(F).set_integrator('dopri5')
    r.set_f_params(entrain_rate, interpTenv, interpTdEnv, interpPress)
    r.set_initial_value(yinit, tinit)
    
    y = np.array(yinit)
    t = np.array(tinit)
    
    #stop tracking the parcel when the time runs out, or if the parcel stops moving/is desecnding
    while r.successful() and r.t < tfin and r.y[0] > 0:
        #find y at the next time step
        #(r.integrate(t) updates the fields r.y and r.t so that r.y = F(t) and r.t = t 
        #where F is the function being integrated)
        r.integrate(r.t+dt)
        if r.y[0] <= 0:
            break
        #keep track of y at each time step
        y = np.vstack((y, r.y))
        t = np.vstack((t, r.t))
        
    wvel = y[:,0]
    cloud_height = y[:,1]
    thetae_cloud = y[:,2]
    wT_cloud = y[:,3]
    
    plt.figure(1)
    plt.plot(wvel, cloud_height)
    plt.xlabel('vertical velocity (m/s)')
    plt.ylabel('height above surface (m)')
    plt.gca().set_title('vertical velocity of a cloud parcel vs height,\
 entrainment rate of %4.1e $s^{-1}$' %entrain_rate)
    
    Tcloud = np.zeros(cloud_height.size)
    wvCloud = np.zeros(cloud_height.size)
    wlCloud = np.zeros(cloud_height.size)
    
    for i in range(0, len(cloud_height)):
        the_press = interpPress(cloud_height[i])*100.
        Tcloud[i], wvCloud[i], wlCloud[i] = tinvert_thetae(thetae_cloud[i], 
                                                        wT_cloud[i], the_press)
        
    Tadia= np.zeros(cloud_height.size)
    wvAdia = np.zeros(cloud_height.size)
    wlAdia = np.zeros(cloud_height.size)
    
   
    for i in range(0, len(cloud_height)):
        the_press = interpPress(cloud_height[i])*100.
        Tadia[i], wvAdia[i], wlAdia[i] = tinvert_thetae(thetae_cloud[0], 
                                                      wT_cloud[0], the_press)
    
    plt.figure(2)
    TcloudHandle, = plt.plot(Tcloud - c.Tc, cloud_height, 'r-')
    TenvHandle, = plt.plot(temp, envHeight, 'g-')
    TadiaHandle, = plt.plot(Tadia - c.Tc, cloud_height, 'b-')
    plt.xlabel('temperature (deg C)')
    plt.ylabel('height above surface (m)')
    plt.gca().set_title('temp. of rising cloud parcel vs height,\
 entrainment rate of %4.1e $s^{-1}$' %entrain_rate)
    plt.gca().legend([TcloudHandle, TenvHandle, TadiaHandle],['cloud', 'environment', 'moist adiabat'])
    plt.show()
    
       
#F returns the buoyancy (and height), velocity, rate of change of thetae_cloud 
#(w.r.t. time) and rate of change of total mixing ratio at a given time step and height
def F(t, y, entrain_rate, interpTenv, interpTdEnv, interpPress):
    yp = np.zeros((4,1))
    velocity = y[0]
    height = y[1]
    thetae_cloud = y[2]
    wT_cloud = y[3]
    #yp[0] is the acceleration, in this case the buoyancy 
    yp[0] = calcBuoy(height, thetae_cloud, interpTenv, interpTdEnv, interpPress)
    press = interpPress(height)*100. #Pa
    Tdenv = interpTdEnv(height) + c.Tc #K
    Tenv = interpTenv(height) + c.Tc #K
    wTenv = wsat(Tdenv, press) #kg/kg
    thetaeEnv = thetaep(Tdenv, Tenv, press)
    #yp[1] is the rate of change of height
    yp[1] = velocity
    #yp[2] is the rate of change of thetae_cloud
    yp[2] = entrain_rate*(thetaeEnv - thetae_cloud)
    #yp[3] is the rate of change of wT_cloud
    yp[3] = entrain_rate*(wTenv - wT_cloud)
    return yp

if __name__ == "__main__":
    answer_entrain()








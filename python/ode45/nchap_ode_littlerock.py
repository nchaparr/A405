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
    Wt = .0135
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

    thetaVal = Tparc0*(c.p0/(press0-pressv0))**(c.Rd/c.cpd)
    thetaeVal = thetaep(wv0, Tparc0, press0)
        
    #aerosol properties
     
    r0, r_a, Num_a = drop_props(RelH0, Tparc0)
    
    Num_rad = len(r0)
    
    yinit = [height0, 0.5, Tparc0, SS0]
    for i in range(len(r0)):
        yinit.append(r0[i])
           
    tinit = 0
    tfin = 100
    dt = 1
   
    r = ode(F).set_integrator('dopri5')

    r.set_f_params(Wt, r0, r_a, Num_a, Num_rad, interpTenv, interpTdEnv, interpPress)
   
   
    r.set_initial_value(yinit, tinit)
   
    y = np.array(yinit)
    t = np.array(tinit)
    
    Press = [press0]
    WSat = [ws0]
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
        Press.append(P)
        
        #test for thetaep and point at which parcel becomes saturated
        
        T = r.y[2]
        [wv, wl] = findWvWl(T, Wt, P)
       
        RelH = r.y[3] + 1
        pressv = esat(T)*RelH
        
        wv1 = c.eps*pressv/(P-pressv)
                
        ws = wsat(T, P)
        WSat.append(ws)

        if ws<Wt:
            
            print''
            print "thetae test", thetaeVal, thetaep(wv1, T, press0)
        else:
            print 'thetat test', thetaVal, T*(c.p0/(P-pressv))**(c.Rd/c.cpd)
            print''
            
        print'Checking bulk and micro wv',wv , wv1, 1.0*(wv - wv1)/wv*100, '%' 
        print'' 
    
        if Wt > ws - .000001 and Wt < ws + .000001:
            print 'becomes saturated at around:' , r.y[0], 'meters'
               
    wvel = y[:,1]
    Tparc = y[:,2]-273.5
    height = y[:,0]
    SS = y[:,3]
    Press = np.array(Press)
    radii = y[:, 4:4+len(r0)]

    fig1 = plt.figure(1)
    plt.clf()
    plt.ylabel('height above surface (m)')
    ax1=fig1.add_subplot(111)
    plt.plot(SS, -Press, 'ro')
    labels = ax1.get_xticklabels()
    for label in labels:
        label.set_rotation(30)
    plt.xlabel('Super Saturation')

    ax2 = ax1.twinx()
    s2 = np.sin(2*np.pi*t)
    ax2.plot(SS, height, 'ro' )
    #plt.xlim(0, .017)

    fig2 = plt.figure(2)
    ax3=fig2.add_subplot(111)
    #plt.plot(radius, -Press, 'o')
    for i in range(len(r0)): 

        plt.semilogx(radii[:,i], -Press)
    labels = ax3.get_xticklabels()
    for label in labels:
       label.set_rotation(30) 
       #plt.legend(loc = 'lower left')
    plt.xlabel('Radii')
    plt.xlim(10**-8, 10**-4)
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

    #distributions of radii: 

    pdf = 1.0*Num_a/np.sum(Num_a)
    fig3 = plt.figure(3)
    plt.clf()
    ax3=fig3.add_subplot(111)
    for i in range(len(Press)):
        plt.semilogx(radii[i,:], pdf)
    plt.xlabel('log plot of radii in meters')
    plt.ylabel('pdf of radii')

    fig4 = plt.figure(4)    
    ax4=fig4.add_subplot(111)
    plt.plot(Tparc, -Press, 'o')
    #labels = ax4.get_xticklabels()
    #for label in labels:
    #   label.set_rotation(30) 
    plt.legend(loc = 'lower left')
    plt.xlabel('Temperature')
    plt.show()
    
def drop_props(RelH0, TempK):
   
    r_a, Num_a, mass_a, prob_a = Aero_dist()
    rho_a = 1775
    
    #r0 = np.array([do_r_find(1.0*RelH0/100, r, rho_a)[0] for r in r_a])
    
    r0 = np.zeros(len(r_a))
    for i in range(len(r_a)):
        
        r0[i] = do_r_find(1.0*RelH0/100, r_a[i], rho_a, TempK)[0]

    return r0, r_a, Num_a 
   
#F returns the buoyancy (and height) and rates of change of Temperature, Droplet Radius and Vapour Pressure with time, at a given time step and height
#def F(t, y, Wt, rho_a, r_a, r0, interpTenv, interpTdEnv, interpPress):
def F(t, y, Wt, r0, r_a, Num_a, Num_rad, interpTenv, interpTdEnv, interpPress):
    #[rho_a, RelH0, Num_rad] = [1775, 90, 24] 
   
    yp = np.zeros((4,1))    
    yp[0] = y[1]
    Vars = calc_Vars(y[0], Wt, r0, r_a, Num_a, y[2], y[1], y[4:4 + Num_rad], y[3], rho_a, 140*10**-3, 3, interpTenv, interpTdEnv, interpPress)
   
    yp[1] = Vars[0] 
    yp[2] = Vars[1]
    
    yp[3] = Vars[2]
    for i in range(len(Vars[3])): 
        #yp.append(Vars[3][i])
        yp = np.append(yp, [[Vars[3][i]]], 0)

    yp = np.array((yp))     

   
    return yp

if __name__ == "__main__":
    ode_littlerock()








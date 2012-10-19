from netCDF4 import Dataset
import site
import sys
site.addsitedir('/nfs/kite/home/nchaparr/A405/repos/a405repo/python//thermlib') #specific to home/nchaparr
site.addsitedir('/nfs/kite/home/nchaparr/A405/repos/a405repo/python//skew_T')
from scipy.integrate import ode
import matplotlib.pyplot as plt
import numpy as np
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

"""
Note: This seemse to be a stiff problem, so is very slow
initially.


Initialization:
         Pulls data from littlerock sounding, adjusts 
         (nudge) height for those values within a 
         tolerance of eachother, defines interpolation 
         function handles (lambda).

         Sets initial values of for parcel variables. 

         Calculates moist and dry adiabats.

         Calls function for setting number and dry 
         radius of aerosols.


Integration:
         Sets initial values, solver type and parameters.
         Integrates F().  Does some tests.     

Handling of Results:
         Creates arrays of key variables.  Plot.

Functions Defined:
         drop_props(RelH0, TempK)
         calls aerodist()
         returns aerosol number and dry radius

         F(t, y, Wt, r0, r_a, Num_a, Num_rad, interpTenv, interpTdEnv, interpPress)
         calls calcVars()
         returns array of differentials for integrator
         
"""
def ode_littlerock():
    filename = 'littlerock.nc'
    print 'reading file: %s\n' %(filename)
    nc_file = Dataset(filename)
    var_names = nc_file.variables.keys()
    
    sound_var = nc_file.variables[var_names[3]]
    press = sound_var[:,0]
    height = sound_var[:,1]
    temp = sound_var[:,2]
    dewpoint = sound_var[:,3]
    RelH = sound_var[:,4]
    
    
    newHeight= nudge(height)
    
    interpTenv = lambda zVals: np.interp(zVals, newHeight, temp)
    interpTdEnv = lambda zVals: np.interp(zVals, newHeight, dewpoint)
    interpPress = lambda zVals: np.interp(zVals, newHeight, press)
    interpRelH = lambda zVals: np.interp(zVals, newHeight, RelH)
    p880_level = np.where(abs(880 - press) < 2.)
    p900_level = np.where(abs(900 - press) < 2.)

    height0 = 885    
    press0 = interpPress(height0)*100
    Tparc0 = 290
    Wt = .0135
    ws0 = wsat(Tparc0, press0)
    es0 = esat(Tparc0)
    pressv0 = (1.0*Wt/(Wt + c.eps))*press0 #from w&h 3.63
    SS0 = 1.0*pressv0/es0 - 1
    RelH0 = SS0 + 1
    
    thetaVal = Tparc0*(c.p0/(press0))**(c.Rd/c.cpd) #w&h 3.54
    thetaeVal = thetaep(ws0, Tparc0, press0)
        
    r0, r_a, Num_a = drop_props(RelH0, Tparc0)
    Num_rad = len(r0)
    
    yinit = [height0, 0.5, Tparc0, SS0]
    for i in range(len(r0)):
        yinit.append(r0[i])
    tinit = 0
    tfin = 100
    dt = .001
   
    r = ode(F).set_integrator('vode', method = 'BDF', order=15, nsteps = 300000)

    r.set_f_params(Wt, r0, r_a, Num_a, Num_rad, interpTenv, interpTdEnv, interpPress)
      
    r.set_initial_value(yinit, tinit)
   
    y = np.array(yinit)
    t = np.array(tinit)
    
    Press = [press0]
    WSat = [ws0]
        
    while r.successful() and r.t < tfin and r.y[1] > 0:
        r.integrate(r.t+dt)
        y = np.vstack((y, r.y))
        t = np.vstack((t, r.t))

        #for tests
        P = interpPress(r.y[0])*100
        Press.append(P)
        T = r.y[2]
        [wv, wl] = findWvWl(T, Wt, P)
        RelH = r.y[3] + 1
        pressv = esat(T)*RelH
        wv1 = c.eps*pressv/(P-pressv)
        ws = wsat(T, P)
        WSat.append(ws)

        #Miscellaneous print statements and tests
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
    plt.ylabel('Height')
    ax1=fig1.add_subplot(111)
    plt.plot(SS, height, 'ro')
    labels = ax1.get_xticklabels()
    for label in labels:
        label.set_rotation(30)
    plt.xlabel('Super Saturation')
   
    fig2 = plt.figure(2)
    ax2=fig2.add_subplot(111)
    for i in range(len(r0)): 
        plt.semilogx(radii[:,i], height)
    labels = ax2.get_xticklabels()
    for label in labels:
       label.set_rotation(30) 
    plt.xlabel('Radii (m)')
    plt.ylabel('Height (m)')
    plt.show()

    pdf = 1.0*Num_a/np.sum(Num_a)
    fig3 = plt.figure(3)
    plt.clf()
    ax3=fig3.add_subplot(111)
    for i in range(len(Press)):
        plt.semilogx(radii[i,:], pdf, 'o' )
    plt.xlabel('Log of radius')
    plt.ylabel('dr(D)/d(lnr)')

    fig4 = plt.figure(4)    
    ax4=fig4.add_subplot(111)
    plt.plot(Tparc, height, 'o')
    plt.xlabel('Temperature (C)')
    plt.ylabel('Height (m)')    
    plt.show()
    
def drop_props(RelH0, TempK):
    r_a, Num_a, mass_a, prob_a = Aero_dist()
    rho_a = 1775
    r0 = np.zeros(len(r_a))
    for i in range(len(r_a)):        
        r0[i] = do_r_find(RelH0, r_a[i], rho_a, TempK)[0]
    return r0, r_a, Num_a 
   
def F(t, y, Wt, r0, r_a, Num_a, Num_rad, interpTenv, interpTdEnv, interpPress):
    yp = np.zeros((4,1))    
    yp[0] = y[1]
    Vars = calc_Vars(y[0], Wt, r0, r_a, Num_a, y[2], y[1], y[4:4 + Num_rad], y[3], rho_a, 140*10**-3, 3, interpTenv, interpTdEnv, interpPress)
    yp[1] = Vars[0] 
    yp[2] = Vars[1]
    yp[3] = Vars[2]
    for i in range(len(Vars[3])): 
        yp = np.append(yp, [[Vars[3][i]]], 0)
    yp = np.array((yp))     
    return yp

if __name__ == "__main__":
    ode_littlerock()








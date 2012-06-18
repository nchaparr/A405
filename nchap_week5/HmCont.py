import numpy as np
from scipy.interpolate import UnivariateSpline
from scipy.interpolate import interp1d
import scipy.optimize
import matplotlib.pyplot as plt
import matplotlib.pylab as pyl
from mccla import mccla
from thermo import qsat

"""This script creates a contour plot of constant Thetav(Virtual Potential Temperature) levels on a an Hm (Moist Static Energy) vs Wtotal (total water) grid. 

First pressure and height values are pulled from fortran script.  These are subjected to an interpolation so that height can be obtained for any given pressure.

Several functions are defined for calculating variables, the most central being that which calculates the Thetav values given Hm and Wtot.  A rootfinder is usedon the function TfromHm to find the correct Temperature and Wv (Water Vapour content).

The number of grid points is give by the user, the grid of Hm and Wt is created and points from Phils convect.m script are plotted at two pressure levels.
"""
#Constants
C2K = 273.15
Lv = 2.5*10**(-6)
g = 9.8
Rd = 287
Cp = 1004
P0 = 10**5

# get pressure and height for one of the mcclatchey standard
# atmospheres
# use f2py to produce mccla.so from mcclatchey.f#
#f2py --overwrite-signature -m mccla -h mccla.pyf mcclatchey.f
#f2py -c mccla.pyf mcclatchey.f

npts=33  #number of points
height,press,temp,rvapmix,o3den,airden = mccla('midsummer',npts)

Readme_variables=['height: height [m]','press: pressure [pa]',\
                  'temp: temperature [K]',\
                  'rvapmix: vapor mixing ratio [kg/kg]',\
                  'o3den: ozone density [kg/m^3]',\
                  'airden: air density [kg/m^3]']
out_dict={'height':height,'press':press,'temp':temp,'rvapmix':rvapmix,\
           'o3den':o3den,'airden':airden,'Readme_variables':Readme_variables}

interp_height = UnivariateSpline(height, np.log(press))# interpolate height against pressure

interp_press = UnivariateSpline(np.log(press), height)

plot_interp = raw_input('Plot interpolation? yes or no:')
if plot_interp == 'yes':
    #Plot to test interpolation
    fig=plt.figure()
    fig.clf()
    ax1=fig.add_subplot(111)
    ax1.plot(height,np.log(press))
    ax1.plot(height,interp_height(height), 'o')
    ax1.plot(interp_press(np.log(press)), np.log(press), '+')
    fig.canvas.draw()
    plt.show()

#Fuctions
def Hmcalc(wt, T, P):
    z = interp_press(np.log(P))
    Hm = Cp*T + Lv*wt + g*z
    return Hm

def Thetavcalc(wv, wl, T, P):
    Theta =  T*(1.0*P0/P)**.286
    Thetav = Theta*(1 + .61*wv - wl)
    return Thetav

def Tvcalc(wv, T):
    Tv = T*(1.0*(wv + .622)/(.622*(1+wv)))
    return Tv

def rhocalc(P, T):
    rho = 1.0*P/(Rd*T)
    return rho

def TdfromHm(Tguess, Hmtarget, z, P):
    Wvguess = Wsatcalc(Tguess, P)
    zero = Hmtarget - (Cp*Tguess + Lv*Wvguess + g*z)
    return zero

def TfromHm(Tguess, wt,  Hmtarget, z, P):
    wsat = Wsatcalc(Tguess, P)
    if wsat <= wt:
        zero = Hmtarget - (Cp*Tguess + Lv*wsat + g*z)
    else:
        zero = Hmtarget - (Cp*Tguess + Lv*wt + g*z)
    return zero

def Wsatcalc(T, P):
    Tc = T - 273.15
    Esat = 611.2*np.exp(1.0*(17.67*Tc)/(Tc + 243.5))
    Wsat = 1.0*(.622*Esat)/(P - Esat) 
    return Wsat

def Getpoints(wt, Hmtarget, z):
    logP  = interp_height(z) 
    P = np.exp(logP)  
    T = scipy.optimize.zeros.brenth(TfromHm, 200, 400, (wt, Hmtarget, z, P))
    wsat = Wsatcalc(T, P)
    
    if wsat<=wt:
        wv = wsat
    else:
        wv = wt
        
    wl = wt - wv
        
    Tv = T*(1 + .622*wv - wl)
    rho = 1.0*P/(Rd*Tv)

    Thetav = Thetavcalc(wv, wl, T, P)
    
    if wt >= wsat:
        tag = 'saturated'
    else:
        tag = 'unsaturated'

    Td = scipy.optimize.zeros.brenth(TdfromHm, 200, 400, (Hmtarget, z, P))
    Tdv = Td*(1 + .622*wt)
    rhosat = 1.0*P/(Rd*Tdv)
    
    return np.array([Hmtarget, rho, Thetav, tag])

# Temperature and total water content from two pressure levels in Phils convec.m
P = np.array([90000,85000])
Temp = [293.15, 296.6907, 300.5701, 288.4194, 291.900, 295.7144]
wtot = [.012, .0093, .0071 ,.012, .0093, .0071]
Hm  = [Hmcalc(wtot[i], Temp[i], 90000) for i in range(3)]
Hm1  = [Hmcalc(wtot[i+3], Temp[i+3], 85000) for i in range(3)]
for i in range(3):
    Hm.append(Hm1[i])

#get number of points for contour plot grid
numpoints = int(raw_input('number of grid points?:'))
for i in range(2):
#create grid of Hm and Wt
    Press = P[i]
    logPress  = np.log(Press)
    z = interp_press(logPress)
    Hmtarget = np.linspace(2.9*10**5, 3.2*10**5, num = numpoints)
    wt = np.linspace(0.00, 0.015, num = numpoints)

#get corresponding thetav and wsat values
#create empty arrays
    Thetav = np.zeros((numpoints, numpoints))
    wsats = np.zeros(numpoints)
#fill arrays  
    for j in range(numpoints):
        T = scipy.optimize.zeros.brenth(TfromHm, 200, 400, (wt[j], Hmtarget[j], z, Press))
        wsats[j] = Wsatcalc(T, Press)
        for k in range(numpoints):
            Thetav[j, k] = Getpoints(wt[j], Hmtarget[k], z)[2]

#create grid for contour plot
    X, Y = np.meshgrid(Hmtarget, wt)
    Thetavlevels = np.mean(Thetav, axis=0)

# Create a contour plot with points from convec.m
    plt.figure(i)
    ConstThetav = plt.contour(X, Y, Thetav, levels = Thetavlevels, colors = 'k')
    plt.clabel(ConstThetav, inline=1, fontsize=10)
    print Hm[0+3*i:3+3*i], wtot[0+3*i:3+3*i]
    plt.plot(Hm[0+3*i:3+3*i], wtot[0+3*i:3+3*i], 'o')
    plt.title("Constant ThetaV Levels " + str(Press) + "Pa")
    plt.show()

    
    

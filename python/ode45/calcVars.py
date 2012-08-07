
import site
site.addsitedir('/nfs/kite/home/nchaparr/A405/repos/a405repo/python/thermlib')
from constants import constants as c
from findTmoist import findTmoist
from wsat import wsat
from findWvWl import findWvWl

def calcBuoy(height, Wt0, thetae0, interpTenv, interpTdEnv, interpPress):

    #input: height (m), thetae0 (K), plus function handles for
    #T,Td, press soundings
    #output: Bout = buoyant acceleration in m/s^2
    #neglect liquid water loading in the virtual temperature
    
    Press=interpPress(height)*100 #Pa
    Tcloud=findTmoist(thetae0,Press) #K
    wvcloud= findWvWl(Tcloud, Wt0, Press)[0]; #kg/kg
    Tvcloud=Tcloud*(1. + c.eps*wvcloud)
    Tenv=interpTenv(height) + c.Tc
    Tdenv=interpTdEnv(height) + c.Tc
    wvenv=wsat(Tdenv,Press); #kg/kg
    Tvenv=Tenv*(1. + c.eps*wvenv)
    TvDiff=Tvcloud - Tvenv
    dW = c.g0*(TvDiff/Tvenv)    
    return dW 

def calcdT(height, Wt0, interpPress, thetae0, Wvel):
    Press=interpPress(height)*100.
    Tcloud=findTmoist(thetae0,Press)
    wvcloud=findWvWl(Tcloud, Wt0, Press)[0]
    wl = Wt0 - wvcloud
    Tvcloud=Tcloud*(1. + c.eps*wvcloud)
    rho= 1.0*Press/(c.Rd*Tvcloud)
    Ws = wsat(Tcloud, Press)#assumed saturated?

    rho= 1.0*Press/(c.Rd*Tvcloud) #density of the cloud

    Pressv = 1.0*(Ws/(Ws + c.eps))*Press
    
    dWs=(Ws + Ws**2)*(1.0*c.lv0/(c.Rv*Tcloud**2))
    
    dT = (1.0*c.Rd*Tcloud/(c.cpd*Press))*Wvel*(-c.g0*rho)*(1.0/(1 - (c.lv0*Ws)/(c.cpd*Tcloud) + (1.0*c.lv0/c.cpd)*dWs))

    a = 1.0*(1+Wt0)/(1+wvcloud*(1.0*c.cpv/c.cpd))
    b = 1+1.0*c.lv0/(c.Rd*Tcloud)
    C = 1.0*wl*c.cl/(c.cpd+wvcloud*c.cpv)
    d = 1.0*c.lv0**2*wvcloud*(1 + 1.0*wvcloud/c.eps)/(c.Rv*Tcloud**2*(c.cpd + wvcloud*c.cpv))

    check1dT = -1.0*(1.0*c.g0*Wvel/c.cpd)*a*b/(1+C+d)
    checkdT = -1.0*(c.g0*Wvel)/(c.cpd + c.lv0*dWs)
    return dT, checkdT, check1dT

def calcdWdT(tempK, pressPa):
    deriv=c.lv0/(c.Rv*tempK**2.)  #1/es des/dT
    thewsat=wsat(tempK,pressPa)
    #chain rule the derivitive of epsilon*esat/(p - esat)
    out= (thewsat[0]  + 1.0*thewsat[0]**2/c.eps)*deriv
    return out

    
    
    

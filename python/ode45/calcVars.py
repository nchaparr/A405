import site
site.addsitedir('/nfs/kite/home/nchaparr/A405/repos/a405repo/python/thermlib')
from constants import constants as c
from findTmoist import findTmoist
from wsat import wsat
from findWvWl import findWvWl
from esat import *
from T_thetaep import t_thetaep

def calcBuoy(height, Wt, Tparc, interpTenv, interpTdEnv, interpPress):

    #input: height (m), thetae0 (K), plus function handles for
    #T,Td, press soundings
    #output: Bout = buoyant acceleration in m/s^2
    #neglect liquid water loading in the virtual temperature
    
    Press = interpPress(height)*100 #Pa
    [wvparc, wlparc] = findWvWl(Tparc, Wt, Press) #kg/kg
    Tvparc = Tparc*(1. + c.eps*wvparc)
    Tenv = interpTenv(height) + c.Tc
    Tdenv = interpTdEnv(height) + c.Tc
    wvenv = wsat(Tdenv,Press) #kg/kg
    Tvenv = Tenv*(1. + c.eps*wvenv)
    TvDiff = Tvparc - Tvenv
    dW = c.g0*(TvDiff/Tvenv)    
    return dW 

def calcdT(height, Wt, Tparc, interpPress, Wvel):
    Press=interpPress(height)*100.
    [wvparc, wlparc] = findWvWl(Tparc, Wt, Press)
    Tvparc=Tparc*(1. + c.eps*wvparc)
    rho = 1.0*Press/(c.Rd*Tvparc)
    Ws = wsat(Tparc, Press)

    rho= 1.0*Press/(c.Rd*Tvparc) #density of the parcel

    Pressv = 1.0*(wvparc/(wvparc + c.eps))*Press
    es = esat(Tparc)
    dWsdT=(Ws + Ws**2)*(1.0*c.lv0/(c.Rv*Tparc**2))
    dWsdP = -1.0*(c.eps*es)/(Press - es)**2

    Pressd = Press - Pressv
    
    dT = (1.0*c.Rd*Tparc/(c.cpd*Pressd) - 1.0*(c.lv0*dWsdP)/(c.cpd))*Wvel*(-c.g0*rho)*(1.0/(1 - (c.lv0*Ws)/(c.cpd*Tparc) + (1.0*c.lv0/c.cpd)*dWsdT))

    a = 1.0*(1+Wt)/(1+wvparc*(1.0*c.cpv/c.cpd))
    b = 1+1.0*c.lv0/(c.Rd*Tparc)
    C = 1.0*wlparc*c.cl/(c.cpd+wvparc*c.cpv)
    d = 1.0*c.lv0**2*wvparc*(1 + 1.0*wvparc/c.eps)/(c.Rv*Tparc**2*(c.cpd + wvparc*c.cpv))

    check1dT = -1.0*(1.0*c.g0*Wvel/c.cpd)*a*b/(1+C+d)
    checkdT = -1.0*((c.g0 + c.lv0*dWsdP*(-rho*c.g0))*Wvel)/(c.cpd + c.lv0*dWsdT)
    
    return dT, checkdT, check1dT

def calcdr(rad, Tparc, height, interpPress, pressv, rho_a, r_a, M_a, I_no, wvel):
    rho_d = 1000 #can i assume this?
    Rg = 8.3143
    Md = 28.97
    rho_l =  1000 #can i assume this?
    sigma_w = .076
    M_w = 18

    v_a = (1.0*4/3)*np.pi*r_a**3 #volume of the aerosol partical in m**3
    m_a = v_a*rho_a #mass of the aerosol in Kg
       
    a = 1.0*2*sigma_w/(rho_d*c.Rv*Tparc)#W&H ex6.12
    b = 1.0*(I_no*m_a*M_w)/((1.0*4/3)*np.pi*M_a*rho_d) 
    
    e_s = esat(Tparc)
    
    e_drop = e_s*(1 + 1.*a/rad  - 1.0*b/rad**3)#W&H ex6.12

    Press=interpPress(height)*100
    
    Ws = wsat(Tparc, Press)

    wvparc = .622*(pressv/(Press - pressv)) 
    Tvparc = Tparc*(1. + c.eps*wvparc)
    rho = 1.0*Press/(c.Rd*Tvparc)
    
    #Pressv = 1.0*(Ws/(Ws + c.eps))*Press

    Dv = 1000.0*2.21*10**-5/Press#Curry&Webster p144 
    dr = (1.0/rad)*(1.0*Dv/(rho_l*c.Rv*Tparc))*(pressv - e_drop)
    dpressv = -(1.0*Press/c.eps)*(rad**2)*(dr) - c.g0*wvel*rho
        
    return dr, dpressv
    
    

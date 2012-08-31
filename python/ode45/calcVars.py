import site
site.addsitedir('/nfs/kite/home/nchaparr/A405/repos/a405repo/python/thermlib')
from constants import constants as c
from findTmoist import findTmoist
from wsat import wsat
from findWvWl import findWvWl
from esat import *
from T_thetaep import t_thetaep

"""
For use in 

Input:
    height (m)
    Wt  = totaly water 
    Tparc  = parcel temperature (K)
    Wvel = vertical velocity (m/s)
    rad = droplet radius (m)
    pressv = vapour pressure (Pa)
    rho_a, r_a, M_a, I_no = density, radius, mass and vant hoff factor for inucleating aerosol. 
    plus function handles for T,Td, press soundings

output: 
    Bout = buoyant acceleration in m/s^2
     = rate of change of temperature
     = rate of change of vapour pressure
     = rate of change of droplet radius
"""

def calc_Vars(height, Wt, Tparc, Wvel, rad, pressv, rho_a, r_a, r0, M_a, I_no, interpTenv, interpTdEnv, interpPress):

    #Should include liquid water loading in the virtual temperature
    #Should test for saturation and replace Ws with Wv
    
    Press = interpPress(height)*100 #Pa
    print 'e', pressv
    wvparc = c.eps*(pressv)/(Press - pressv)
    wlparc = Wt - wvparc
    
    Tvparc = Tparc*(1. + c.eps*wvparc)
    Tenv = interpTenv(height) + c.Tc
    Tdenv = interpTdEnv(height) + c.Tc
    wvenv = wsat(Tdenv,Press) #kg/kg
    Tvenv = Tenv*(1. + c.eps*wvenv)
    TvDiff = Tvparc - Tvenv
    Ws = wsat(Tparc, Press)
    es = esat(Tparc)
    rho= 1.0*Press/(c.Rd*Tvparc) #density of the parcel
    Pressd = Press - pressv

    dW = c.g0*(TvDiff/Tvenv) #acceleration or buoyancy
        

    dWsdT=(Ws + Ws**2)*(1.0*c.lv0/(c.Rv*Tparc**2))#see if there's another derivable dT for the unsaturated adiabat, OR just assume the parcel is never unsaturated
    dWsdP = -1.0*(c.eps*es)/(Press - es)**2

    if Wt < Ws:
        dT = -(1.0*c.Rd*Tparc/(c.cpd*Pressd))*(c.g0*rho*Wvel)
    else:   
        dT = (1.0*c.Rd*Tparc/(c.cpd*Pressd) - 1.0*(c.lv0*dWsdP)/(c.cpd))*Wvel*(-c.g0*rho)*(1.0/(1 - (c.lv0*Ws)/(c.cpd*Tparc) + (1.0*c.lv0/c.cpd)*dWsdT))
    
    a = 1.0*(1+Wt)/(1+wvparc*(1.0*c.cpv/c.cpd))
    b = 1+1.0*c.lv0/(c.Rd*Tparc)
    C = 1.0*wlparc*c.cl/(c.cpd+wvparc*c.cpv)
    d = 1.0*c.lv0**2*wvparc*(1 + 1.0*wvparc/c.eps)/(c.Rv*Tparc**2*(c.cpd + wvparc*c.cpv))

    check1dT = -1.0*(1.0*c.g0*Wvel/c.cpd)*a*b/(1+C+d)
    checkdT = -1.0*((c.g0 + c.lv0*dWsdP*(-rho*c.g0))*Wvel)/(c.cpd + c.lv0*dWsdT)

    rho_d = 1000 #kg/m**3 can i assume this?
    Rg = 8.3143
    rho_l =  1000 #can i assume this?
    sigma_w = .076
    M_w = 18*10**-3 #kg/mol
    v_a = (1.0*4/3)*np.pi*r_a**3 #volume of the aerosol partical in m**3
    m_a = v_a*rho_a #mass of the aerosol in Kg
    N_a = 10**9 #number/m**3

     #Test on total water according to adiabats
    testwv, testwl = findWvWl(Tparc, Wt, Press)    

    A = 1.0*2*sigma_w/(rho_d*c.Rv*Tparc)#W&H ex6.12
    B = 1.0*(I_no*m_a*M_w)/((1.0*4/3)*np.pi*M_a*rho_d) 
    e_drop = es*(1 + 1.*A/rad  - 1.0*B/rad**3)#W&H ex6.12
     
    Dv = 1000.0*(2.21*10**-5)/Press #Diffusion constant, Curry&Webster p144 
    #rad will be an array
    if rad > r0:
        dr = (1.0/rad)*(1.0*Dv/(rho_l*c.Rv*Tparc))*(pressv - e_drop)
    else:
        dr = 0            
    #will have three arrays here: Num_a, rad, and dr, so will need to do two np.multiplies() to get a product array for first term here 
    dpressv =  - ((N_a*4*np.pi*rad**2*rho_l)*(dr) + (c.eps*pressv*c.g0*Wvel*rho/(Press-pressv)**2))/(c.eps*Press/(Press - pressv)**2)
    #print 'dr, dpressv', dr
    #print 'e_drop, A, B', e_drop, A*Tparc, B*(M_a/(m_a*I_no)), (M_a/(m_a*I_no))**-1

    print 'Liquid Water: total water, micro vapor, micro liquid , bulk', Wt, wvparc, wlparc, testwl
    print '' 
    print 'drop mass, aerosol mass, drop radius, aerosol radius, initial frop radius, dr', (1.0*4/3)*np.pi*(rad**3)*rho_l, (1.0*4/3)*np.pi*(r_a**3)*rho_a, rad, r_a, r0, dr  
    print ''
    print 'increase in water on droplets', N_a*(1.0*4/3)*np.pi*(rad**3 - r0**3)*rho_l
    print ''
    print 'Water Vapour', wvparc    
    
    return dW, dT, checkdT, check1dT, dr, dpressv 


    
    
    

import site
site.addsitedir('/nfs/kite/home/nchaparr/A405/repos/a405repo/python/thermlib')
from constants import constants as c
from findTmoist import findTmoist
from wsat import wsat
from findWvWl import findWvWl
from esat import *
from T_thetaep import t_thetaep
from nchap_ode_littlerock import drop_props
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

#def calc_Vars(height, Wt, Tparc, Wvel, rad, SS, rho_a, r_a, r0, M_a, I_no, interpTenv, interpTdEnv, interpPress):
def calc_Vars(height, Wt, Tparc, Wvel, rad, SS, rho_a, RelH0, M_a, I_no, interpTenv, interpTdEnv, interpPress):
    print''
    print'now CalcVars is running...'
    print''
    #Should include liquid water loading in the virtual temperature
    #Should test for saturation and replace Ws with Wv
    print''
    print''
    print''
    r0, r_a, Num_a = drop_props(RelH0)
    Press = interpPress(height)*100 #Pa
    es = esat(Tparc)
    pressv = (SS + 1)*es
    #print 'e', pressv
    wvparc = c.eps*(pressv)/(Press - pressv)
    wlparc = Wt - wvparc
    testwv,testwl = findWvWl(Tparc, Wt, Press)
    Tvparc = Tparc*(1. + c.eps*testwv)
    Tenv = interpTenv(height) + c.Tc
    Tdenv = interpTdEnv(height) + c.Tc
    wvenv = wsat(Tdenv,Press) #kg/kg
    Tvenv = Tenv*(1. + c.eps*wvenv)
    TvDiff = Tvparc - Tvenv
    Ws = wsat(Tparc, Press)
    
    rho= 1.0*Press/(c.Rd*Tvparc) #density of the parcel
    Pressd = Press - pressv

    dW = c.g0*(TvDiff/Tvenv) #acceleration or buoyancy
    
    dWsdT=(Ws + Ws**2)*(1.0*c.lv0/(c.Rv*Tparc**2))#see if there's another derivable dT for the unsaturated adiabat, OR just assume the parcel is never unsaturated
    dWsdP = -1.0*(c.eps*es)/(Press - es)**2
    
    rho_d = 1000 #kg/m**3 can i assume this?
    Rg = 8.3143
    rho_l =  1000 #can i assume this?
    sigma_w = .076
    M_w = 18*10**-3 #kg/mol
    v_a = (1.0*4/3)*np.pi*r_a**3 #volume of the aerosol partical in m**3
    m_a = v_a*rho_a #mass of the aerosol in Kg
    N_a = 10**9 #number/m**3

     #Test on total water according to adiabats
     
    dWsdT=(Ws + Ws**2)*(1.0*c.lv0/(c.Rv*Tparc**2))
    
    dWsdP = -1.0*(c.eps*es)/(Press - es)**2

    rho_d = 1000 #kg/m**3 can i assume this?
    Rg = 8.3143
    rho_l =  1000 #can i assume this?
    sigma_w = .076
    M_w = 18*10**-3 #kg/mol
    v_a = (1.0*4/3)*np.pi*r_a**3 #volume of the aerosol partical in m**3
    m_a = v_a*rho_a #mass of the aerosol in Kg
    N_a = 10**9 #number/m**3

     #Test on total water according to adiabats
    
    rhovs = es/(c.Rv*Tparc) 
    Dv = 1000.0*(2.21*10**-5)/Press #Diffusion constant, Curry&Webster p144 
    #rad will be an array
  
    v_a = []
    dr = []
    dwl = []
    print ''
    print 'entering loop to calculate dr'
    print ''
    for i in range(len(rad)):
        print''
        print'', r_a[i]
        print''
        v_a.append((1.0*4/3)*np.pi*r_a[i]**3)
        if rad[i] <= r0[i] and (1.0/rad[i])*(1.0*Dv*rhovs/rho_l)*(SS) < 0:
            drad = 0  
        else:
            drad = (1.0/rad[i])*(1.0*Dv*rhovs/rho_l)*(SS)
        dr.append(drad)
        dwl.append(Num_a[i]
                   *rho_l*4.0*np.pi*(rad[i]**2)*drad)
    cumdwl = np.cumsum(np.array(dwl))
    dwv = -cumdwl  
 
    #if rad <= r0 and (1.0/rad)*(1.0*Dv*rhovs/rho_l)*(SS) < 0:        dr = 0
    #else:
    #dr = (1.0/rad)*(1.0*Dv*rhovs/rho_l)*(SS)

         #print 'initial radius, radius, dr', r0, rad, dr
        
    dwl = N_a*rho_l*4.0*np.pi*(rad**2)*dr                
    dwv = -dwl
    des = 1.0*(c.lv0*es)/(c.Rv*(Tparc**2)) 

    #will have three arrays here: Num_a, rad, and dr, so will need to do two np.multiplies() to get a product array for first term here

    if Wt < Ws:
        dT = (-1.0*c.g0/(c.cpd))*(Wvel)
        #dT = -(c.Rd*Tparc/(c.cpd*Press))*c.g0*rho*Wvel
    else:   
        dT = -(1.0*c.lv0)/(c.cpd)*dwv - (1.0*c.g0)/(c.cpd)*Wvel 
        #dT = ((1.0*c.Rd*Tparc)/(c.cpd*Pressd) - (1.0*c.lv0/c.cpd)*dWsdP)/(1 - (c.lv0*Ws/(c.cpd*Tparc)) + (c.lv0*dWsdT)/(c.cpd))*Wvel*(-rho*c.g0)
        
    term1 = (-c.eps*es/Press**2)*(-c.g0*Press*Wvel/(c.Rd*Tparc)) 
    term2 = (c.eps/Press)*(c.eps*es*c.lv0/(c.Rd*Tparc**2))*dT

    dSS =  (Press/(c.eps*es))*(dwv - (1+SS)*(term1 + term2))
    print''
    print'Now CalcVars is finished.'
    print''               
    return dW, dT, dSS, dr 


    
    
    

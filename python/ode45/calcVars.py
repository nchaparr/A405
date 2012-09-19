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
    dW = buoyant acceleration in m/s^2
    dT = rate of change of temperature
    dr = rate of change of 
    dSS = rate of change of 

    """

#def calc_Vars(height, Wt, Tparc, Wvel, rad, SS, rho_a, r_a, r0, M_a, I_no, interpTenv, interpTdEnv, interpPress):
def calc_Vars(height, Wt, r0, r_a, Num_a, Tparc, Wvel, rad, SS, rho_a, M_a, I_no, interpTenv, interpTdEnv, interpPress):

    Press = interpPress(height)*100 #Pa
    es = esat(Tparc)
    pressv = (SS + 1)*es
   
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
    
    dWsdT=(Ws + Ws**2)*(1.0*c.lv0/(c.Rv*Tparc**2))
    
    dWsdP = -1.0*(c.eps*es)/(Press - es)**2

    rho_w = 1000
    rho_d = 1000 #kg/m**3 can i assume this?
    Rg = 8.3143
    rho_l =  1000 #can i assume this?
    sigma_w = .076
    M_w = 18*10**-3 #kg/mol

    #Test on total water according to adiabats
     
    dWsdT=(Ws + Ws**2)*(1.0*c.lv0/(c.Rv*Tparc**2))
    
    dWsdP = -1.0*(c.eps*es)/(Press - es)**2

    rhovs = es/(c.Rv*Tparc) 
    
    Dv = 100000.0*(2.21*10**-5)/Press #Diffusion constant, Curry&Webster p144 
    
    #v_a = []
    dr = []
    dwl = []
    for i in range(len(rad)):
        #v_a.append((1.0*4/3)*np.pi*r_a[i]**3)
        a, b = aero_prms(sigma_w, I_no, M_a, M_w, rho_a, r_a[i], rad[i], rho_w, Tparc)
        #if rad[i] <= r0[i] and (1.0/rad[i])*(1.0*Dv*rhovs/rho_l)*(SS - 1.0*a/rad[i] + 1.0*b/rad[i]**3) < 0:
            drad = 0  
            #else:
        drad = (1.0/rad[i])*(1.0*Dv*rhovs/rho_l)*(SS - 1.0*a/rad[i] + 1.0*b/rad[i]**3)
        print 'drad', drad     
        dr.append(drad)
        dwl.append(Num_a[i]*rho_l*4.0*np.pi*(rad[i]**2)*drad)
    cumdwl = np.sum(np.array(dwl))
    dwv = -cumdwl 
                    
    des = 1.0*(c.lv0*es)/(c.Rv*(Tparc**2)) 
   

    if Wt < Ws:
        #dT = (-1.0*c.g0/(c.cpd))*(Wvel)
        #dT = -(1.0*c.Rd*Tparc/(c.cpd*Pressd))*c.g0*rho*Wvel
        dT = -(1.0*c.Rd*Tparc/(c.cpd*Pressd))*(c.g0*rho*Wvel)
    else:   
        #dT = -(1.0*c.lv0)/(c.cpd)*dwv - (1.0*c.g0)/(c.cpd)*Wvel 
        #dT = ((1.0*c.Rd*Tparc)/(c.cpd*Pressd) - (1.0*c.lv0/c.cpd)*dWsdP)/(1 - (c.lv0*Ws/(c.cpd*Tparc)) + (c.lv0*dWsdT)/(c.cpd))*Wvel*(-rho*c.g0)

        dT = (1.0*(c.Rd*Tparc)/(c.cpd*Pressd) - 1.0*(c.lv0*dWsdP)/(c.cpd))*Wvel*(-c.g0*rho)*(1.0/(1 - (c.lv0*Ws)/(c.cpd*Tparc) + (1.0*c.lv0/c.cpd)*dWsdT))
        
    term1 = (-c.eps*es/Press**2)*(-c.g0*Press*Wvel/(c.Rd*Tparc)) 
    term2 = (c.eps/Press)*(des)*dT

    dSS =  (Press/(c.eps*es))*(dwv - (1+SS)*(term1 + term2))
    print'dSS', dSS    
    return dW, dT, dSS, dr 

def aero_prms(sigma_w, I_no, M_a, M_w, rho_a, r_a, r, rho_w, TempK):
    
    v_d = 1.0*4/3*np.pi*(r**3)
    v_a = 1.0*4/3*np.pi*(r_a**3)
    m_a = rho_a*v_a
    rho_d = 1.0*(m_a + (v_d - v_a)*rho_w)/v_d
        
    a = 1.0*2*sigma_w/(rho_w*c.Rv*TempK)
    b = 1.0*(I_no*m_a*M_w)/((1.0*4/3)*M_a*np.pi*rho_w)    

    #print 'test on a and b', a*TempK, 1.0*b/I_no*(1.0*M_a/m_a)

    return a, b

    
    
    

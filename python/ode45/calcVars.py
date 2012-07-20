
import site
site.addsitedir('C:\Users\Den\mya405\python\\thermlib')
from constants import constants as c
from findTmoist import findTmoist
from new_thermo import wsat

def calcBuoy(height, thetae0, interpTenv, interpTdEnv, interpPress):

    #input: height (m), thetae0 (K), plus function handles for
    #T,Td, press soundings
    #output: Bout = buoyant acceleration in m/s^2
    #neglect liquid water loading in the virtual temperature
    
    press=interpPress(height)*100.#%Pa
    Tcloud=findTmoist(thetae0,press) #K
    wvcloud=wsat(Tcloud,press); #kg/kg
    Tvcloud=Tcloud*(1. + c.eps*wvcloud)
    Tenv=interpTenv(height) + c.Tc
    Tdenv=interpTdEnv(height) + c.Tc
    wvenv=wsat(Tdenv,press); #kg/kg
    Tvenv=Tenv*(1. + c.eps*wvenv)
    TvDiff=Tvcloud - Tvenv
    #print '%10.3f %10.3f %10.3f\n' %(press*0.01,height,TvDiff)
    return c.g0*(TvDiff/Tvenv)

# 
# W&H Exercise 3.9 A parcel of air with an initial temperature
# of 15 °C and dew point 2 °C is lifted adiabatically
# from the 1000-hPa level. Determine its LCL
# and temperature at that level. If the air parcel is
# lifted a further 200 hPa above its LCL, what is its
# final temperature and how much liquid water is condensed
# during this rise?
import site
site.addsitedir('C:\Users\Den\mya405\python\\thermlib')
from constants import constants as c
from new_thermo import wsat, thetaep, tinvert_thetae
from findLCL0 import findLCL0
from findTmoist import findTmoist


press0=1.e5
Temp0=15 + c.Tc
Td0=2 + c.Tc
#
#step 1: find the lcl
#
wv0=wsat(Td0,press0)
plcl, Tlcl =findLCL0(wv0,press0,Temp0)

print 'found Plcl=%8.2f (hPa) and Tlcl=%8.2f (deg C)\n' %(plcl*1.e-2, Tlcl - c.Tc)

the_thetae=thetaep(Tlcl,Tlcl,plcl)

# step 2  raise air 200 hPa above plcl along pseudo adiabat,
# find wsat at that temperature, compare to wv0 to find the amount
# of liquid water condensed

pnew = plcl - 200.e2

newTemp, newWv, newWl = tinvert_thetae(the_thetae, wv0, pnew)

#check:
#ws = wsat(newTemp, pnew)
#wl = wv0 - ws

out_msg = '\nat new pressure level %8.2f hPa\n\
          the temperature is  %8.2f deg C\n\
          the vapor mixing ratio is %8.2f g/kg\n\
          and there are %8.4f g/kg of liquid water\n\n'

print out_msg %(pnew*1.e-2, newTemp - c.Tc, newWv*1.e3, newWl*1e3)

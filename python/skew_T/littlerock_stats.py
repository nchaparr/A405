import site
site.addsitedir('C:\Users\Den\mya405\python\\thermlib')
import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from constants import constants as c
from new_thermo import convertTempToSkew, thetaep, wsat, tinvert_thetae
from convecSkew import convecSkew
from calcAdiabat import calcAdiabat
from calcTvDiff import calcTvDiff
from nudge import nudge
    
filename='littlerock.nc';
print 'reading file: %s\n' %filename
nc_file = Dataset(filename)
#nc_file.variables returns a dictionary with variable names as keys and
#variable instances as values
var_names = nc_file.variables.keys()
print var_names
print nc_file.ncattrs()   
print nc_file.units
print nc_file.col_names


#grab the March 2 12Z sounding
sound_var = nc_file.variables[var_names[3]]
print 'found sounding: \n%s\n' %(sound_var)
    
press= sound_var[:,0]
temp= sound_var[:,2]
dewpoint= sound_var[:,3]
direct= sound_var[:,6]
speed= sound_var[:,7]

plt.figure(1)
plt.semilogy(temp, press)
plt.semilogy(dewpoint, press)
plt.gca().invert_yaxis()
plt.gca().set_ybound((200, 1000))
plt.gca().set_title('littlerock sounding, %s' %var_names[3])
plt.ylabel('pressure (hPa)')
plt.xlabel('temperature (deg C)')
plt.show()

plt.figure(2)
skew, ax2 = convecSkew(2)
xtemp=convertTempToSkew(temp,press,skew)
xdew=convertTempToSkew(dewpoint,press,skew)
plt.semilogy(xtemp,press,'g-',linewidth=4)
plt.semilogy(xdew,press,'b-',linewidth=4)

#use 900 hPa sounding level for adiabat
#array.argmin() finds the index of the min. value of array 
p900_level = np.abs(900 - press).argmin()

thetaeVal=thetaep(dewpoint[p900_level] + c.Tc,temp[p900_level] + c.Tc,press[p900_level]*100.)
pressVals,tempVals =calcAdiabat(press[p900_level]*100.,thetaeVal,200.e2)
xTemp=convertTempToSkew(tempVals - c.Tc,pressVals*1.e-2,skew)
p900_adiabat = plt.semilogy(xTemp,pressVals*1.e-2,'r-',linewidth=4)
xleft=convertTempToSkew(-20,1.e3,skew)
xright=convertTempToSkew(25.,1.e3,skew)

#
#interpolator fails if two pressure values are the same
#
newPress = nudge(press)
#independent variable used to interpolate must be in increasing order 
#pVals must be in hPa
interpTenv = lambda pVals: np.interp(pVals, newPress[::-1], temp[::-1])
interpTdenv = lambda pVals: np.interp(pVals, newPress[::-1], dewpoint[::-1])
interpDirec = lambda pVals: np.interp(pVals, newPress[::-1], direct[::-1])
interpSpeed = lambda pVals: np.interp(pVals, newPress[::-1], speed[::-1])
trytemp = interpTenv(pressVals*1.e-2)
xTemp=convertTempToSkew(trytemp,pressVals*1.e-2,skew)
plt.semilogy(xTemp,pressVals*1.e-2,'k.',markersize=6)
labels = range(100, 1100, 100)
ax2.set_yticks(labels)
ax2.set_yticklabels(labels)
ax2.set_ybound((200, 1000))
ax2.set_xbound((xleft, xright))
plt.gca().set_title('littlerock sounding, %s' %var_names[3])
#plt.gca().legend([p900_adiabat], ['900 hPa moist adiabat'])
plt.show()

calcTvDiffHandle = lambda pVals: calcTvDiff(pVals, thetaeVal, interpTenv, interpTdenv)
presslevs = np.linspace(200, press[0], 100)*1e2
#start integrating from first sounding level
presslevs = presslevs[::-1]
Tvdiff = [calcTvDiffHandle(p) for p in presslevs]
    
plt.figure(3)
plt.plot(Tvdiff, presslevs/100)
plt.title('virtual temperature difference vs. pressure (hPa)')
plt.gca().invert_yaxis()
plt.show()

cumCAPE = -c.Rd*np.cumsum(Tvdiff[1:]*np.diff(np.log(presslevs)))

plt.figure(4)
plt.plot(cumCAPE, presslevs[1:]/100)
plt.title('cumulative CAPE (J/kg) vs. pressure (hPa)')
plt.gca().invert_yaxis()
plt.show()

#equate kinetic and potential energy to get maximum
#updraft speed
   
plt.figure(5)
maxvel=np.sqrt(2*cumCAPE);
plt.plot(maxvel, presslevs[1:]*0.01,'k-');
plt.title('maximum updraft (m/s) vs. pressure (hPa)');
plt.gca().invert_yaxis()
plt.show()

#
# find storm indices
#
#  lifted index 
thetaeVal=thetaep(dewpoint[0] + c.Tc,temp[0] + c.Tc,press[0]*100.)
wT=wsat(dewpoint[0] + c.Tc,press[0]*100.)
Tadia_500,wv,wl=tinvert_thetae(thetaeVal,wT,500.e2)
Temp_500=interpTenv(500.) + c.Tc
lifted_index= Temp_500 - Tadia_500
# total totals = vertical totals plus cross totals
Temp_850=interpTenv(850.) + c.Tc
dew_850 = interpTdenv(850.) + c.Tc
TT_index=Temp_850 + dew_850 - 2*Temp_500
#  Sholwater
thetaeVal=thetaep(dew_850,Temp_850,850.*100.)
wT=wsat(dew_850,850*100.)
Tadia_500,wv,wl=tinvert_thetae(thetaeVal,wT,500.e2)
sholwater=Temp_500 - Tadia_500
#  SWEAT
speed_850=interpSpeed(850.)
speed_500=interpSpeed(500.)
dir_850=interpDirec(850.)
dir_500=interpDirec(500.)
angle=(dir_500 - dir_850);
shear=125.*(np.sin(angle*(np.pi/180)) + 0.2)
term1 = 12.*(dew_850 - c.Tc)
term2 = 20*(TT_index - 49)
if term1 < 0:
    term1 = 0
if TT_index < 49:
    term2 = 0
if not (dir_850 >= 130 and dir_850 <= 250):
    shear = 0
if not (dir_500 >= 150 and dir_500 <= 310):
    shear = 0
if dir_500 - dir_850 <= 0:
    shear = 0
if not (speed_850 and speed_500 >= 15):
    shear = 0
sweat=term1 + term2 + 2.*speed_850 + speed_500 + shear
out_mesg='\nCAPE is %10.4f (J/kg)\n\
\nLifted index= %10.4f (K)\n\
\nTotal Totals=%10.4f (K)\n\
\nSholwater=%10.4f (K)\n\
\nSWEAT =%10.4f'

print out_mesg %(cumCAPE[-1],lifted_index,TT_index,sholwater,sweat)







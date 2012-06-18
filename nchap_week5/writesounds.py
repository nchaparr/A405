# get temperature and height for one of the mcclatchey standard
# atmospheres
# use f2py to produce mccla.so from mcclatchey.f
#
#f2py --overwrite-signature -m mccla -h mccla.pyf mcclatchey.f
#f2py -c mccla.pyf mcclatchey.f
from scipy.interpolate import interp1d
#from scipy.io import savemat
from mccla import mccla
from thermo import qsat
np=33  #number of points
height,press,temp,rvapmix,o3den,airden = mccla('midsummer',np)
import matplotlib.pyplot as plt
Readme_variables=['height: height [m]','press: pressure [pa]',\
                  'temp: temperature [K]',\
                  'rvapmix: vapor mixing ratio [kg/kg]',\
                  'o3den: ozone density [kg/m^3]',\
                  'airden: air density [kg/m^3]']
#Readme_variables=['press: pressure [pa]',\
 #                 'temp: temperature [K]']
out_dict={'height':height,'press':press,'temp':temp,'rvapmix':rvapmix,\
           'o3den':o3den,'airden':airden,'Readme_variables':Readme_variables}
#out_dict={'press':press,'temp':temp,'Readme_variables':Readme_variables}
fig=plt.figure(1)
fig.clf()
ax1=fig.add_subplot(111)
ax1.plot(height,press)
#ax1.semilogy(height,press)
ax1.invert_yaxis()
plt.xlim(0, 2000)
fig.canvas.draw()
plt.show()



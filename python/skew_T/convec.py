#calculate thermodynamic variables
#for a convectively unstable layer
import site
site.addsitedir('C:\Users\Den\mya405\python\\thermlib')
from constants import constants as c
from convecSkew import convecSkew
from new_thermo import thetaep, wsat, tinvert_thetae, convertTempToSkew, LCLfind, Tdfind
from findTmoist import findTmoist
import numpy as np
import matplotlib.pyplot as plt


Tbot=20.
Ttop=17.5
Tdewbot=15
Tdewtop=5.5

Pbot=900.
Ptop=800.

#calculate the temperature and dewpoint sounding
#assuming linear profiles
slope=(Ttop - Tbot)/(Ptop - Pbot)
press=np.arange(Pbot, Ptop-10, -10)
Temp=(Tbot + slope*(press - Pbot));
slope=(Tdewtop - Tdewbot)/(Ptop - Pbot)
Tdew = (Tdewbot + slope*(press - Pbot))

#how big is the pressure vector?
numPoints= np.size(press);


#figure 1: cloud base at 900 hPa
plt.figure(1)
skew, ax = convecSkew(1)
#zoom the axis to focus on layer
skewLimits=convertTempToSkew([5,30],1.e3,skew)
plt.axis([skewLimits[0],skewLimits[1],1000,600])
xplot1=convertTempToSkew(Temp,press,skew)
#plot() returns a list of handles for each line plotted
Thandle, = plt.plot(xplot1,press,'k-', linewidth=2.5)
xplot2=convertTempToSkew(Tdew,press,skew)
TdHandle, = plt.plot(xplot2,press,'b--', linewidth=2.5)
plt.title('convectively unstable sounding: base at 900 hPa')

plt.show()
#print -dpdf initial_sound.pdf

#put on the top and bottom LCLs and the thetae sounding
Tlcl=np.zeros(numPoints)
pLCL=np.zeros(numPoints)
theTheta=np.zeros(numPoints)
theThetae=np.zeros(numPoints)
Tpseudo=np.zeros(numPoints)
wtotal=np.zeros(numPoints)
for i in range(0, numPoints):
  wtotal[i]=wsat(Tdew[i] + c.Tc,press[i]*100.);
  Tlcl[i],pLCL[i]=LCLfind(Tdew[i] + c.Tc,Temp[i]+c.Tc,press[i]*100.)
  theThetae[i]=thetaep(Tdew[i] + c.Tc,Temp[i] + c.Tc,press[i]*100.)
  #find the temperature along the pseudo adiabat at press[i]
  Tpseudo[i]=findTmoist(theThetae[i],press[i]*100.)
  #no liquid water in sounding
xplot=convertTempToSkew(Tlcl[0] - c.Tc,pLCL[0]*0.01,skew);
bot,=plt.plot(xplot,pLCL[0]*0.01,'ro',markersize=12, markerfacecolor ='r')
xplot=convertTempToSkew(Tlcl[-1] - c.Tc,pLCL[-1]*0.01,skew)
top,=plt.plot(xplot,pLCL[-1]*0.01,'bd',markersize=12,markerfacecolor='b')
#print -dpdf initial_lcls.pdf
xplot=convertTempToSkew(Tpseudo - c.Tc,press,skew)
thetaEhandle,=plt.plot(xplot,press,'c-', linewidth=2.5)
ax.legend([Thandle, TdHandle, bot, top, thetaEhandle], ['Temp (deg C)','Dewpoint (deg C)',
       'LCL bot (835 hPa)','LCL top (768 hPa)','$\\theta_e$'])
plt.title('convectively unstable sounding: base at 900 hPa')
plt.show()
#print -dpdf base900_thetae.pdf


#figure 2: lift cloud base by 50 hPa to 850 hPa 
for i in range(0,numPoints):
  press[i]=press[i] - 50
  #find the temperature along the pseudoadiabats at press[i]
  Tpseudo[i]=findTmoist(theThetae[i],press[i]*100.)
  #find the actual temperature and dewpoint
  Temp[i],wv,wl=tinvert_thetae(theThetae[i],wtotal[i],press[i]*100.)
  Tdew[i]=Tdfind(wv,press[i]*100.)

plt.figure(2)
skew, ax = convecSkew(2)
skewLimits=convertTempToSkew([5,30],1.e3,skew)
plt.axis([skewLimits[0],skewLimits[1],1000,600])
xplot1=convertTempToSkew(Temp - c.Tc,press,skew)
Thandle,=plt.plot(xplot1,press,'k-', linewidth=2.5)
xplot2=convertTempToSkew(Tdew -c.Tc,press,skew)
TdHandle,=plt.plot(xplot2,press,'b--', linewidth=2.5)
xplot=convertTempToSkew(Tpseudo - c.Tc,press,skew)
thetaEhandle,=plt.plot(xplot,press,'c-', linewidth=2.5)
xplot=convertTempToSkew(Tlcl[0] - c.Tc,pLCL[0]*0.01,skew);
bot,=plt.plot(xplot,pLCL[0]*0.01,'ro',markersize=12, markerfacecolor ='r')
xplot=convertTempToSkew(Tlcl[-1] - c.Tc,pLCL[-1]*0.01,skew)
top,=plt.plot(xplot,pLCL[-1]*0.01,'bd',markersize=12,markerfacecolor='b')
xplot=convertTempToSkew(Tpseudo - c.Tc,press,skew)
thetaEhandle,=plt.plot(xplot,press,'c-', linewidth=2.5)
ax.legend([Thandle, TdHandle, bot, top, thetaEhandle], ['Temp (deg C)','Dewpoint (deg C)',
       'LCL bot (835 hPa)','LCL top (768 hPa)','$\\theta_e$'])
plt.title('convectively unstable sounding: base at 850 hPa')
plt.show()
#print -dpdf base_850.pdf


#figure 3: lift cloud base by 1470 Pa to 835.3 hPa
for i in range(0,numPoints):
  press[i]=press[i] - 14.7
  #find the temperature along the pseudoadiabats at press[i]
  Tpseudo[i]=findTmoist(theThetae[i],press[i]*100.)
  #find the actual temperature and dewpoint
  Temp[i],wv,wl=tinvert_thetae(theThetae[i],wtotal[i],press[i]*100.)
  Tdew[i]=Tdfind(wv,press[i]*100.)

plt.figure(3)
skew, ax = convecSkew(3)
skewLimits=convertTempToSkew([5,30],1.e3,skew)
plt.axis([skewLimits[0],skewLimits[1],1000,600])
xplot1=convertTempToSkew(Temp - c.Tc,press,skew)
Thandle,=plt.plot(xplot1,press,'k-', linewidth=2.5)
xplot2=convertTempToSkew(Tdew -c.Tc,press,skew)
TdHandle,=plt.plot(xplot2,press,'b--', linewidth=2.5)
xplot=convertTempToSkew(Tpseudo - c.Tc,press,skew)
thetaEhandle,=plt.plot(xplot,press,'c-', linewidth=2.5)
xplot=convertTempToSkew(Tlcl[0] - c.Tc,pLCL[0]*0.01,skew);
bot,=plt.plot(xplot,pLCL[0]*0.01,'ro',markersize=12, markerfacecolor ='r')
xplot=convertTempToSkew(Tlcl[-1] - c.Tc,pLCL[-1]*0.01,skew)
top,=plt.plot(xplot,pLCL[-1]*0.01,'bd',markersize=12,markerfacecolor='b')
xplot=convertTempToSkew(Tpseudo - c.Tc,press,skew)
thetaEhandle,=plt.plot(xplot,press,'c-', linewidth=2.5)
ax.legend([Thandle, TdHandle, bot, top, thetaEhandle], ['Temp (deg C)','Dewpoint (deg C)',
       'LCL bot (835 hPa)','LCL top (768 hPa)','$\\theta_e$'])
plt.title('convectively unstable sounding: base at 835.3 hPa')
plt.show()
#print -dpdf base_835.pdf


#figure 4: lift cloud base by 10.30 hPa to 825 hPa
for i in range(0,numPoints):
  press[i]=press[i] - 10.3
  #find the temperature along the pseudoadiabats at press[i]
  Tpseudo[i]=findTmoist(theThetae[i],press[i]*100.)
  #find the actual temperature and dewpoint
  Temp[i],wv,wl=tinvert_thetae(theThetae[i],wtotal[i],press[i]*100.)
  Tdew[i]=Tdfind(wv,press[i]*100.)

plt.figure(4)
skew, ax = convecSkew(4)
skewLimits=convertTempToSkew([5,30],1.e3,skew)
plt.axis([skewLimits[0],skewLimits[1],1000,600])
xplot1=convertTempToSkew(Temp - c.Tc,press,skew)
Thandle,=plt.plot(xplot1,press,'k-', linewidth=2.5)
xplot2=convertTempToSkew(Tdew -c.Tc,press,skew)
TdHandle,=plt.plot(xplot2,press,'b--', linewidth=2.5)
xplot=convertTempToSkew(Tpseudo - c.Tc,press,skew)
thetaEhandle,=plt.plot(xplot,press,'c-', linewidth=2.5)
xplot=convertTempToSkew(Tlcl[0] - c.Tc,pLCL[0]*0.01,skew);
bot,=plt.plot(xplot,pLCL[0]*0.01,'ro',markersize=12, markerfacecolor ='r')
xplot=convertTempToSkew(Tlcl[-1] - c.Tc,pLCL[-1]*0.01,skew)
top,=plt.plot(xplot,pLCL[-1]*0.01,'bd',markersize=12,markerfacecolor='b')
xplot=convertTempToSkew(Tpseudo - c.Tc,press,skew)
thetaEhandle,=plt.plot(xplot,press,'c-', linewidth=2.5)
ax.legend([Thandle, TdHandle, bot, top, thetaEhandle], ['Temp (deg C)','Dewpoint (deg C)',
       'LCL bot (835 hPa)','LCL top (768 hPa)','$\\theta_e$'])
plt.title('convectively unstable sounding: base at 825 hPa')
plt.show()
#print -dpdf base_825.pdf


#figure 5: lift cloud base by 25 hPa to 800 hPa
for i in range(0,numPoints):
  press[i]=press[i] - 25
  #find the temperature along the pseudoadiabats at press[i]
  Tpseudo[i]=findTmoist(theThetae[i],press[i]*100.)
  #find the actual temperature and dewpoint
  Temp[i],wv,wl=tinvert_thetae(theThetae[i],wtotal[i],press[i]*100.)
  Tdew[i]=Tdfind(wv,press[i]*100.)

plt.figure(5)
skew, ax = convecSkew(5)
skewLimits=convertTempToSkew([5,30],1.e3,skew)
plt.axis([skewLimits[0],skewLimits[1],1000,600])
xplot1=convertTempToSkew(Temp - c.Tc,press,skew)
Thandle,=plt.plot(xplot1,press,'k-', linewidth=2.5)
xplot2=convertTempToSkew(Tdew -c.Tc,press,skew)
TdHandle,=plt.plot(xplot2,press,'b--', linewidth=2.5)
xplot=convertTempToSkew(Tpseudo - c.Tc,press,skew)
thetaEhandle,=plt.plot(xplot,press,'c-', linewidth=2.5)
xplot=convertTempToSkew(Tlcl[0] - c.Tc,pLCL[0]*0.01,skew);
bot,=plt.plot(xplot,pLCL[0]*0.01,'ro',markersize=12, markerfacecolor ='r')
xplot=convertTempToSkew(Tlcl[-1] - c.Tc,pLCL[-1]*0.01,skew)
top,=plt.plot(xplot,pLCL[-1]*0.01,'bd',markersize=12,markerfacecolor='b')
xplot=convertTempToSkew(Tpseudo - c.Tc,press,skew)
thetaEhandle,=plt.plot(xplot,press,'c-', linewidth=2.5)
ax.legend([Thandle, TdHandle, bot, top, thetaEhandle], ['Temp (deg C)','Dewpoint (deg C)',
       'LCL bot (835 hPa)','LCL top (768 hPa)','$\\theta_e$'])
plt.title('convectively unstable sounding: base at 800 hPa')
plt.show()
#print -dpdf base_800.pdf


#figure 6: lift cloud base by 32.25 hPa to 768 hPa
for i in range(0,numPoints):
  press[i]=press[i] - 32.25
  #find the temperature along the pseudoadiabats at press[i]
  Tpseudo[i]=findTmoist(theThetae[i],press[i]*100.)
  #find the actual temperature and dewpoint
  Temp[i],wv,wl=tinvert_thetae(theThetae[i],wtotal[i],press[i]*100.)
  Tdew[i]=Tdfind(wv,press[i]*100.)

plt.figure(6)
skew, ax = convecSkew(6)
skewLimits=convertTempToSkew([5,30],1.e3,skew)
plt.axis([skewLimits[0],skewLimits[1],1000,600])
xplot1=convertTempToSkew(Temp - c.Tc,press,skew)
Thandle,=plt.plot(xplot1,press,'k-', linewidth=2.5)
xplot2=convertTempToSkew(Tdew -c.Tc,press,skew)
TdHandle,=plt.plot(xplot2,press,'b--', linewidth=2.5)
xplot=convertTempToSkew(Tpseudo - c.Tc,press,skew)
thetaEhandle,=plt.plot(xplot,press,'c-', linewidth=2.5)
xplot=convertTempToSkew(Tlcl[0] - c.Tc,pLCL[0]*0.01,skew);
bot,=plt.plot(xplot,pLCL[0]*0.01,'ro',markersize=12, markerfacecolor ='r')
xplot=convertTempToSkew(Tlcl[-1] - c.Tc,pLCL[-1]*0.01,skew)
top,=plt.plot(xplot,pLCL[-1]*0.01,'bd',markersize=12,markerfacecolor='b')
xplot=convertTempToSkew(Tpseudo - c.Tc,press,skew)
thetaEhandle,=plt.plot(xplot,press,'c-', linewidth=2.5)
ax.legend([Thandle, TdHandle, bot, top, thetaEhandle], ['Temp (deg C)','Dewpoint (deg C)',
       'LCL bot (835 hPa)','LCL top (768 hPa)','$\\theta_e$'])
plt.title('convectively unstable sounding: base at 768 hPa')
plt.show()
#print -dpdf base_768.pdf
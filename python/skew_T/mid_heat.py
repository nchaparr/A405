import site
site.addsitedir('C:\Users\Den\mya405\python\\thermlib')
import matplotlib.pyplot as plt
import numpy as np
from constants import constants
from new_thermo import Tdfind, thetaep, tinvert_thetae, convertTempToSkew
from convecSkew import convecSkew

c=constants()
wtA=14.e-3
pressA=900.e2
tempA=25 + c.Tc

TdA=Tdfind(wtA,pressA)
thetaeA=thetaep(TdA,tempA,pressA)
wtB=wtA
pressB=700.e2
TdB=Tdfind(wtB,pressB)
thetaeB=thetaeA
tempB,wvB,wlB=tinvert_thetae(thetaeB, wtB, pressB)

wtC=wtA
pressC=900.e2
TdC=Tdfind(wtC,pressC)
tempC=tempB
thetaeC=thetaep(TdC,tempC,pressC)

plt.figure(1)
skew, ax1 = convecSkew(1)
xtempA=convertTempToSkew(tempA - c.Tc,pressA*0.01,skew)
xtempB=convertTempToSkew(tempB - c.Tc,pressB*0.01,skew)
xtempC=convertTempToSkew(tempC - c.Tc,pressC*0.01,skew)
plt.text(xtempA,pressA*0.01,'A', fontweight='bold',fontsize= 22, color='b')
plt.text(xtempB,pressB*0.01,'B', fontweight='bold',fontsize= 22,color='b')
plt.text(xtempC,pressC*0.01,'C', fontweight='bold',fontsize= 22,color='b')

pressLevs=np.linspace(700,900,60)*100.
lineAB = np.zeros(pressLevs.size)
rhoAB = np.zeros(pressLevs.size)

#adiabatic expansion from A to B
for i in range(0, len(pressLevs)):
    thePress=pressLevs[i]
    temp,wv,wl=tinvert_thetae(thetaeA, wtA, thePress)
    lineAB[i]=temp
    rho=thePress/(c.Rd*temp)
    rhoAB[i]=rho

#isothermal compression from B to C
rhoBC=pressLevs/(c.Rd*tempB)

#isobaric heating from C to A 
tempCA=np.linspace(tempC,tempA,100)
rhoCA=pressA/(c.Rd*tempCA)
press900Vec=np.ones(rhoCA.size)*pressA


xtemp=convertTempToSkew(lineAB - c.Tc,pressLevs*0.01,skew)
plt.semilogy(xtemp,pressLevs*0.01,'k-.',linewidth=2)
plt.semilogy([xtempB,xtempC],[700.,900.],'b-.',linewidth=2);
plt.semilogy([xtempC,xtempA],[900.,900.],'r-.',linewidth=2)
xleft=convertTempToSkew(10,1000.,skew)
xright=convertTempToSkew(30,1000.,skew)
labels = np.array(range(100, 1100, 100))
ax1.set_yticks(labels)
ax1.set_yticklabels(labels)
ax1.set_ybound((650, 1000))
ax1.set_xbound((xleft, xright))
plt.title('heat engine problem')
plt.show()
#print -depsc prob1.eps

plt.figure(2)
plt.plot(1./rhoCA,press900Vec*1.e-2, 'r')
plt.ylim([700,1000.])
plt.plot(1./rhoAB,pressLevs*1.e-2, 'k')
plt.plot(1./rhoBC,pressLevs*1.e-2, 'b')
plt.xlabel('specific volume ($m^3 kg^{-1}$)')
plt.ylabel('pressure (hPa)')
plt.title('volume - pressure plot')
plt.show()

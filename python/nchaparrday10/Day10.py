import numpy as np

eps = .622
Rd = 287

print 'excercise 3.39'
dT = 5
Cp = 4.2*10**3
Lv = 2.5*10**6

waterpercentage = dT*1.0*Cp/Lv*100
print "mass of water as a percentage of mass of person ", waterpercentage

print 'excercise 3.40'

Rv = 461.51
esat = 2300
e1 = .6*esat
T = 20 + 273.15
V1 = 20
V2 = 4

rho1 = e1/(Rv*T)
print rho1

rho2 = (1.0*V2/V1)*rho1
print rho2

e2 = rho2*Rv*T
print e2

econd = esat - e2
print econd

rhocond = econd/(Rv*T)
print rhocond

masscond = rhocond*V2
print masscond

print 'excercise 4.41'

q = .0196
T = 30 + 273.15
w = 1.0*q/(1-q)
Tv = T*1.0*(w + .622)/(.622*(1+w))
pressure = 1014*10**2

print 'Tv', Tv

#pressure = Rd*density*Tv

density = 1.0*pressure/(Rd*Tv)
print'density', density

print 'excercise 4.42'

p = 97500
T = 15 + 273.15
w = 1.8*10**(-3)

q = 1.0*w/(eps+w)*p
print "Vapour pressure", q
Tv = T*1.0*(w + eps)/( eps*(1 + w))
print 'Tv', Tv

print 'excercise 4.43'

wsat = 8.7*10**(-3)
T = 18 + 273.15
Tw = 12 + 273.15
Lv = 2.25*10**6
Cp = 1004
Cpw = 1952

#(wsat - w)*Lv = (T - Tw)*(Cp + *Cpw*w)

w = (Cp*(Tw - T) + wsat*Lv)/(Cpw*(T-Tw) + Lv)

print "w", w


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
e1 = 0.6*esat
print 'initial e', e1
T = 20 + 273.15
V1 = 20*10**(-3)
V2 = 4*10**(-3)

rho1 = e1/(Rv*T)
print 'initial density of water vapour', rho1

mass1 = rho1*V1
print 'initial mass of water vapour', mass1

rho2 = 1.0*esat/(Rv*T)
mass2 = rho2*V2
print 'mass after decreasing volume, and increasing pressure above wsat', mass2

print 'mass condensed out', mass1 - mass2

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


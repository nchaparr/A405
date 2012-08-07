import numpy as np

""" 
    Estimates the number and mass of aerosols within
    a series of bins.  Functions for nd and md are as derived in 
    Problem set 4
"""

def nd(D):#number of aerosols for a certain diameter
    rho = 1775
    nd = (rho**-.25)*(1.0941*10**3)*(D)**(-7.0/4)
    return nd

def md(D):#mass of aerosol with a certain diameter
    rho = 1775
    md = (rho**.75)*(.5728*10**3)*(D**(5.0/4.0))
    return md
#get a series of 6 increments of the logD over a given range
logDs = np.linspace(np.log10(10**-8), np.log10(10**-7), 6)
print 'log(D)', logDs

#get the corresponding diameters
Ds = [10**x for x in logDs]
print "diameters", Ds

#now get 5 corresponding diameter differences, between the above 
dD = [0]
for i in range(len(Ds) - 1):
    dD.append(Ds[i+1] - Ds[i])
dD = np.array(dD)

#calculated number and mass in each bin.  See (12) in Aerosol notes
nDs = np.array([nd(D) for D in Ds])
mDs = np.array([md(D) for D in Ds])
num = np.multiply(nDs, dD)
mass = np.multiply(mDs, dD)

print "number and mass concentrations for 5 bins", num, mass




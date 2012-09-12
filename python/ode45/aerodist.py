import numpy as np
import matplotlib.pyplot as plt 

"""
For use in

Gives the number and mass of aerosol of a particular diameter as determined by the lognormal distribution

"""
"""
LogDist(): 
   out=logDist(themass,Amp,meanval,logsd)
   return a vector with the mass (kg) in each bin specified by
   themass, using the lognormal distribution given in the
   aerosol notes.  Input: themass (kg), meanval(kg), logsd
   unitless)
   utput:  kg/logbinwidth, that is, out(i)*(log(themass(i+1)) - log(themass(i)))   
   is the aerosol mass, in kg in bin(i), for aerosols with logmss
   between log(themass(i)) and log(themass(i+1))
"""
  
def logDist(themass, Amp, meanval, logsd):
    theexp=np.exp(-0.5*(((np.log(themass) - np.log(meanval))/logsd))**2.);
    out=Amp/(np.sqrt(2*np.pi)*logsd)*theexp;
    print 'array of masses', themass
    print 'probabilities', out
    return out 

def Aero_dist():
# setup the underlying  lognormal mass distribution

#distribution parameters
    rhoaero = 1775
    Mtot =.25e-9 #0.25 mug/m^3
    logsd = 1.7
    meandiam = 0.2*10**-6
    meanmass = 1.0*rhoaero*4/3*np.pi*(meandiam/2.)**3.

#mass bins from 0.01 to 0.5 microns
#start with log10(mass)
    logdrydiam = np.linspace(-2, .1, num = 24)
#convert to meters
    drydiam = 10.0**logdrydiam*1.0e-6
    dryrad = 1.0*drydiam/2
    
#convert aerosol volume to kg
    themass = rhoaero*(1.0*4/3)*np.pi*(dryrad)**3.0

#get the number distribtuion for themass vector
    out = logDist(themass, Mtot, meanmass, logsd)
    veclength = len(out)

#now divide this by the bin mass, which gives
#the number of aerosols in each bin
    thenum = np.divide(1.0*out, themass)

    end = len(themass)
    print 'end', end
# now work backwards to check our number distribution
    dlogMass = np.subtract(np.log(themass[1:end]), np.log(themass[0:end-1]))
    print 'check total mass:'
    binmass = np.array([x*y for x in out[0:end-1]for y in dlogMass])
    totmass=np.sum(binmass)
    print totmass
    print 'check total number'
    end = len(thenum)
    binnum = np.array([x*y for x in thenum[1:end] for y in dlogMass])
    totnum = np.sum(binnum)
    print totnum

    return dryrad, thenum, themass, out
#plot the distribution

if __name__ == "__main__":

    dryrad, thenum, themass, out = Aero_dist()
    fig=plt.figure(1)
    fig.clf()
    ax1 = fig.add_subplot(111)
    #plt.semilogx(themass,out)
    plt.plot(themass, out)
    #plt.title('lognormal mass distribution')
    plt.xlim(-10**-17, 10**-15)
    #plt.ylabel('m x mu(m) (kg/m^3/per log binwidth)')
    plt.show()

   

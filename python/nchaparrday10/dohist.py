from netCDF3 import Dataset
from scipy.interpolate import UnivariateSpline
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.mlab as mlab

"""This script:
   takes mean pressure, pressure perturbation, temperature, and water content variables
   from an nc file
   calculates virtual temperature and density as well as their respective means and
   perturbations
   and histogram-plots two slices of each perturbation array and divides them by the corresponding mean. 
"""

#Gas Constant
Rd = 287

#function for getting the variable to be plotted as a histogram

def getslice(h, var, meanvec):
    ind = np.where(height==h)[0][0]
    varslice = var[ind,:,:]
    mean = meanvec[ind]
    varslice = np.ravel(1.0*varslice/mean)
    return varslice

#beginning the process of abstracting data from the nc file
filename='subset.nc';
nc_file=Dataset(filename, 'r')
var_names=nc_file.variables.keys()

print "variable names: ",var_names
print "global attributes: ",nc_file.ncattrs()
print 'z', nc_file.variables['z'][...]

#absracting data from nc file as arrays

height = nc_file.variables['z'][...]
T = nc_file.variables['TABS'][...]
qv = nc_file.variables['QV'][...]
ql = nc_file.variables['QN'][...]
pmean = nc_file.variables['p'][...]
ppert = nc_file.variables['PP'][...]

#Calculating virtual temperature and density means and perturbations

Tv = np.zeros(ppert.shape)
Tv = np.multiply(np.add(0.622*qv, -ql) + 1, T)
Tvmean = np.mean(np.mean(Tv, axis=1), axis=1)

Tvpert = np.zeros(ppert.shape)
p = np.zeros(ppert.shape)

for i in range(ppert.shape[0]):
    p[i,:,:] = ppert[i,:,:] + pmean[i]

rho = np.multiply(p, 1.0/(Rd*Tv))
rhomean = np.mean(np.mean(rho, axis = 1), axis = 1)

rhopert = np.zeros(ppert.shape)
for i in range(ppert.shape[0]):
    rhopert[i,:,:] = rho[i,:,:] - rhomean[i]
    Tvpert[i,:,:] = rho[i,:,:] - Tvmean[i]    

#calculating desired variables and plotting the histograms

sliceh = [810, 1510]
varlist = [[rhopert, rhomean], [Tvpert,Tvmean], [ppert, p]]
titlelist = ["Density Perturbation divided by Mean", "Tv Perturbation divided by Mean", "Pressure Perturbation divided by Mean"]
count = 0
for i in range(len(sliceh)):
    for j in range(len(varlist)):
            var =  np.ravel(getslice(sliceh[i], varlist[j][0], varlist[j][1]))
            print var
            plt.figure(count)
            n, bins, patches = plt.hist(var, 50, normed=1, histtype='stepfilled')
            plt.setp(patches, 'facecolor', 'g', 'alpha', 0.75)
            plt.title(titlelist[j] + ' at ' + str(sliceh[i]) + ' meters' )
            plt.show()
            count = count+1

 


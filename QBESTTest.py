# -*- coding: utf-8 -*-
"""
Created on Thu Mar 22 20:53:53 2018

@author: alexa
"""

# -*- coding: utf-8 -*-
"""
Created on Tue Mar 13 16:21:11 2018

@author: ahwch
"""

from scipy import genfromtxt
from Bootstrap import Bootstrap
import numpy as np
from math import sqrt
from numpy import matlib
from copy import deepcopy
from scipy.stats import norm
from scipy.spatial.distance import mahalanobis


tnspec = genfromtxt("tnspec.csv",delimiter=",").T
newspec = genfromtxt("newspec.csv",delimiter=",")

'''
BTRAIN, CNTER, BSAMP = Bootstrap(tnspec,100)

BSAMP = BSAMP[0]
'''

b = np.size(BTRAIN,0)

qdist = b

RADFRAC = .1

SENSITIV = .01

residuals = np.subtract(newspec,CNTER.T)

SO2 = np.sqrt(np.sum((np.sum(residuals**2,axis=0))))

SOR = np.sqrt(np.sum((BTRAIN-matlib.repmat(CNTER,b,1))**2,axis=1))
 
S2R = np.sqrt(np.sum((BTRAIN-matlib.repmat(newspec,b,1))**2,axis=1))

sub = (SO2+SOR+S2R)/2

area = np.sqrt((sub*(sub-SO2)*(sub-SOR)*(sub-S2R)))

radial = (2*area)/SO2

project = np.sqrt(SOR**2-radial**2)

locsarray = (SO2**2)+(SOR**2)

locs = np.where(locsarray < S2R**2)

project[locs] = project[locs]*-1

qdist = deepcopy(project)

qrr = np.sort(radial)

rindex = int(b * RADFRAC)-1

radii=qrr[rindex]

rrloc = np.where(radial > radii)

qdist[rrloc] = 0

rrloc = np.where( qdist != 0)

qdist = qdist[rrloc]
 
qdist = np.sort(qdist)


lindex = np.floor(.16*len(qdist))

uindex = np.floor(.84*len(qdist))

if len(qdist) < 50:
    print("Need more replicates in hypercylinder")

qdisttest = deepcopy(qdist)

sd = np.std(qdist)*sqrt(len(tnspec))
    
sds = np.sqrt(np.sum(CNTER-newspec)**2)/sd

alpha = norm.cdf(-1,loc=0,scale=1)

za = norm.ppf(alpha,loc=0,scale=1)

tcenter = np.median(tnspec,axis=0)

csO2 = deepcopy(SO2)

csOR= sqrt(np.sum((tcenter-CNTER)**2))

cs2R = sqrt(np.sum((tcenter-newspec)**2))

csub = (csO2+csOR+cs2R)/2

carea = np.sqrt((csub*(csub-csO2)*(csub-csOR)*(csub-cs2R)))

cradial = 2 * carea/csO2

cproject = sqrt(csOR**2 - cradial**2)

if((SO2**2 + csOR**2) > cs2R**2):
    cproject = -cproject
    
n = len(qdist)

if(np.floor(n/2) == n/2):
    div = int(n/2) - 1
    div2 = int(div + 1)
    
    md = (qdist[div] + qdist[div2])/2
    
else:
    tocall = np.floor((n/2)+.5)
    tocall = int(tocall)
    md = qdist[tocall]
    
    
cproject = cproject*SENSITIV+md

fdist = qdist-cproject

index = np.arange(1,len(fdist)+1)

if(cproject > np.max(qdist)):
    zelement = len(qdist)-1
    
if(cproject < np.min(qdist)):
    zelement = 1
    
else:
    rootloc = np.where(abs(fdist) == min(abs(fdist)))
    zelement = rootloc[0]
    
z0 = 2.5333/len(qdist)

if(abs(2*z0) > abs(za)):
    print("Decrease skew sensitivity")
    
sensitiv = abs(SENSITIV)

lowind = np.floor(norm.cdf((2*z0+za),loc=0,scale=1)*len(qdist))
upind = np.floor(norm.cdf((2*z0-za),loc=0,scale=1)*len(qdist))

if(lowind < 2):
    print ("WARNING TOO FEW REPLICATES")

if(upind > len(qdist)-2):
    print ("WARNING TOO FEW REPLICATES")
    
if(lowind < 1):
    lowind = 1
    
if(upind > len(qdist)):
    upind = len(qdist)
    
lowind = int(lowind)
upind = int(upind)
    
lowlim = qdist[lowind-1 ]
uplim = qdist[upind-1]

euc = np.sqrt(np.sum(np.power((CNTER-newspec),2)))

fac = abs(norm.ppf(alpha))

erd = np.sqrt(len(tnspec))

if(abs(2*z0)> abs(fac)):
    print("WARNING, SKEW CORRECTION exceeds replicates")
    
sdskew = euc/((uplim/fac)*erd)


    

    









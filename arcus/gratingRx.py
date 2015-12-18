import numpy as np
from numpy import sin,cos,exp,sqrt,pi,tan

x,y = np.meshgrid(np.linspace(-75./2,75./2,1000),\
                  np.linspace(11500-96./2,11500+96./2,1000))

#Compute d spacing for 8.4 m case, 12 m case, and infinite case
rho84 = sqrt(x**2 + (y-3600.)**2)
d84 = 160./(11500.-3600)*rho84

rho12 = sqrt(x**2 + y**2)
d12 = 160./11500.*rho12

#Discretize 8.4 m d spacing
d84[np.logical_and(d84<160.5,d84>160.)] = 160.25
d84[np.logical_and(d84<161.,d84>160.5)] = 160.75
d84[np.logical_and(d84<160.,d84>159.5)] = 159.75
d84[np.logical_and(d84<159.5,d84>159.)] = 159.25

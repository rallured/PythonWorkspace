from numpy import *
import matplotlib.pyplot as plt
import matplotlib
from os import *

def angulardependence(theta,T):
    theta = double(theta)
    wk = double(.038) #Fluorescence yield of K shell photons
    l1 = double(31.9) #Attenuation length of 5893 eV photon in microns
    l2 = double(9.197) #Attenuation length of 1486 eV photon in microns

    #return (wk/(4*np.pi))*(l2*np.cos(theta))/ \
    #       (l2*np.cos(theta)+l1)*np.sin(theta)

##    const = (wk/(4*pi))*l2*cos(theta)*sin(theta)/(l1-l2*cos(theta))
##    tdep = exp(-T/l1)-exp(-T/(l2*cos(theta)))
##    jacob = sin(theta)
##    return const*tdep*jacob
##    return (l1*wk*sin(theta)*(l1 - exp(-(T/l1))*l1 + \
##                              (-1 + exp(-(T/(cos(theta)*l2))))*\
##            l2*cos(theta))) / (l2*4*pi*(l1 - l2*cos(theta)))
    return (sin(theta)/4/pi)*((exp(-(T/l1)) - exp(-((T/cos(theta))/l2)))*
   l1*wk*cos(theta))/(l1 - l2*cos(theta))

def blah(t1,t2,dphi,thick):
    #Create inclination angle array
    theta = arange(t1,t2+.001,.001)
    theta = double(theta*pi/180)
    dtheta = theta[1]-theta[0]

    rates = zeros(size(thick))

    for i in range(size(thick)):
        angdep = angulardependence(theta,thick[i])
        rates[i] = 1/sum(angdep*dtheta*dphi)

    return rates

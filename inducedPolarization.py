import numpy as np

###CXC Ir optical constants
##cxcIr = np.transpose(np.genfromtxt('/home/rallured'
##                        '/Dropbox/AXRO/WSTracing/chandraConstants.txt')[2:])

def readNK(fn):
    """
    Read in an optical constant file from Windt's library
    delta = 1-n
    beta = k
    """
    d = np.transpose(np.genfromtxt(fn,skip_header=8))
    ener = 1.24/(d[0]/10.)
    delta = 1 - d[1]
    beta = d[2]
    return [ener,delta,beta]

def computeReflectivities(ang,energy,rough,density,constants):
    """Return reflectivity with 0.5 nm RMS roughness
    calculated using Fresnel equations, Chandra
    optical constants, and "Strehl" factor (Nevot-Croce/Debye-Waller)
    ang in radians
    Energy supplied in eV
    Roughness in RMS nm
    Density as fraction of bulk density
    constants is a list returned by readNK
    """
    #Get proper optical constants
    if np.size(energy) == 1:
        ind = np.argmin(abs(constants[0]-energy/1000.))
        b = constants[2][ind]*density
        d = constants[1][ind]*density
    else:
        b,d = np.zeros(len(energy)),np.zeros(len(energy))
        for i in range(len(energy)):
            ind = np.argmin(abs(constants[0]-energy[i]/1000.))
            b[i] = constants[2][ind]*density
            d[i] = constants[1][ind]*density
    n = 1 - d + 1j*b
    n2 = abs(n)**2
    #Compute reflectivity in each polarization plane
    #Return mean value
    Rp = abs(n*n*np.sin(ang)-np.sqrt(n*n-np.cos(ang)**2))**2/\
         abs(n*n*np.sin(ang)+np.sqrt(n*n-np.cos(ang)**2))**2
    Rs = abs(np.sin(ang)-np.sqrt(n*n-np.cos(ang)**2))**2/\
         abs(np.sin(ang)+np.sqrt(n*n-np.cos(ang)**2))**2
    R = np.mean([Rp,Rs],axis=0)
    wave = 1240./energy #wavelength in nm
    k = 2*np.pi/wave
    strehl = np.exp(-4*k**2*np.sin(ang)**2*rough**2)

    return R*strehl,Rp*strehl,Rs*strehl

def computeInducedPolarization(ang,energy,rough,density,constants):
    R,Rp,Rs = computeReflectivities(ang,energy,rough,density,constants)
    return (Rs-Rp)/(Rp+Rs)

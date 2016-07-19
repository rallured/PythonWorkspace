import numpy as np

def Rp(n1,n2,ang):
    """Compute the amplitude reflectivity for p polarization
    as a function of complex refractive index and
    incidence angle from surface
    n1 is incident medium
    n2 is transmitted medium"""
    return -(n2*n2*np.sin(ang)-n1*np.sqrt(n2*n2-n1*n1*np.cos(ang)**2)) /\
           (n2*n2*np.sin(ang)+n1*np.sqrt(n2*n2-n1*n1*np.cos(ang)**2))

def Rs(n1,n2,ang):
    """Compute the amplitude reflectivity for s polarization
    as a function of complex refractive index and
    incidence angle from surface
    n1 is incident medium
    n2 is transmitted medium"""
    return (n1*np.sin(ang)-np.sqrt(n2*n2-n1*n1*np.cos(ang)**2)) /\
           (n1*np.sin(ang)+np.sqrt(n2*n2-n1*n1*np.cos(ang)**2))

def Tp(n1,n2,ang):
    """Compute the amplitude transmission for p polarization
    as a function of complex refractive index and
    incidence angle from surface
    n1 is incident medium
    n2 is transmitted medium"""
    return 2*n1*n2*np.sin(ang)/\
           (n2*n2*np.sin(ang)+n1*np.sqrt(n2*n2-n1*n1*np.cos(ang)**2))

def Ts(n1,n2,ang):
    """Compute the amplitude transmission for s polarization
    as a function of complex refractive index and
    incidence angle from surface
    n1 is incident medium
    n2 is transmitted medium"""
    return 2*n1*np.sin(ang)/\
           (n1*np.sin(ang)+np.sqrt(n2*n2-n1*n1*np.cos(ang)**2))

def Ra(n1,n2,ang):
    """Compute the amplitude reflectivity for unpolarized
    light as a function of complex refractive index
    and incidence angle from surface
    n1 is incident medium
    n2 is transmitted medium"""
    return np.mean([Rs(n1,n2,ang),Rp(n1,n2,ang)],axis=0)

def Ta(n1,n2,ang):
    """Compute the amplitude transmission for unpolarized
    light as a function of complex refractive index
    and incidence angle from surface
    n1 is incident medium
    n2 is transmitted medium"""
    return np.mean([Ts(n1,n2,ang),Tp(n1,n2,ang)],axis=0)

def phaseFactor(n,wave,alpha,thick):
    """Compute a phase factor for Parratt's recursive method
    for a given layer. n is complex index of refraction, wave
    is the wavelength in mm, thickness is in mm, alpha is
    graze angle in radians
    """
    k = 2*np.pi/wave*n
    return np.exp(1j*k*np.sin(alpha)*thick)


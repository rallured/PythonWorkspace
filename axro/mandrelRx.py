import numpy as np
import matplotlib.pyplot as plt
import traces.conicsolve as conic
from numpy.polynomial.legendre import legval

def parabola():
    """
    Compute parameters and sag table for primary parabola
    Produce plots of profile and sag profile (1st order removed)
    """
    z = np.linspace(8475.-100,8475.+100.,1000)
    r = conic.primrad(z,220.,8400.)
    fig = plt.figure(figsize=[20.7,8])
    fig.add_subplot(121)
    plt.plot(z,r)
    plt.title('Primary Parabolic Profile')
    plt.xlabel('Axial Height (mm)')
    plt.ylabel('Radius (mm)')
    plt.grid()
    fig.add_subplot(122)
    plt.plot(z,r-np.polyval(np.polyfit(z,r,1),z))
    plt.title('Primary Parabolic Profile (Cone Removed)')
    plt.xlabel('Axial Height (mm)')
    plt.ylabel('Radius (mm)')
    plt.grid()

    return z,r
    
def ellipse():
    """
    As above, but with ellipsoid
    """
    psieff = 1.1844342872222389
    S = 91110.
    p,a,b,e,f = conic.ellipsoidFunction(S,psieff,220.,8400.)

    z = np.linspace(8475.-100,8475.+100.,1000)
    r = conic.ellipsoidRad(S,psieff,220.,8400.,z)
    fig = plt.figure(figsize=[20.7,8])
    fig.add_subplot(121)
    plt.plot(z,r)
    plt.title('Primary Ellipsoid Profile')
    plt.xlabel('Axial Height (mm)')
    plt.ylabel('Radius (mm)')
    plt.grid()
    fig.add_subplot(122)
    plt.plot(z,r-np.polyval(np.polyfit(z,r,1),z))
    plt.title('Primary Ellipsoid Profile (Cone Removed)')
    plt.xlabel('Axial Height (mm)')
    plt.ylabel('Radius (mm)')
    plt.grid()

    return z,r

def hyperbola():
    """
    As above, but with hyperbola
    """
    z = np.linspace(8325.-100,8325.+100.,1000)
    r = conic.secrad(z,220.,8400.)
    fig = plt.figure(figsize=[20.7,8])
    fig.add_subplot(121)
    plt.plot(z,r)
    plt.title('Secondary Hyperbolic Profile')
    plt.xlabel('Axial Height (mm)')
    plt.ylabel('Radius (mm)')
    plt.grid()
    fig.add_subplot(122)
    plt.plot(z,r-np.polyval(np.polyfit(z,r,1),z))
    plt.title('Secondary Hyperbolic Profile (Cone Removed)')
    plt.xlabel('Axial Height (mm)')
    plt.ylabel('Radius (mm)')
    plt.grid()

    return z,r

def doublehyperbola():
    """
    Hyperbola plus sag of primary
    """
    #extra sag required for secondary
    primsag = conic.primsag(8600.,220.,8400.)
    #Initial parabla
    z = np.linspace(8325.-100,8325.+100.,1000)
    r = conic.secrad(z,220.,8400.)
    #Sag perturbation
    rpert = legval(np.linspace(-1.,1.,1000),\
                   [0,0,-primsag/1.5])
    r = r+rpert
    #Make plot
    fig = plt.figure(figsize=[20.7,8])
    fig.add_subplot(121)
    plt.plot(z,r)
    plt.title('Secondary Hyperbolic2 Profile')
    plt.xlabel('Axial Height (mm)')
    plt.ylabel('Radius (mm)')
    plt.grid()
    fig.add_subplot(122)
    plt.plot(z,r-np.polyval(np.polyfit(z,r,1),z))
    plt.title('Secondary Hyperbolic2 Profile (Cone Removed)')
    plt.xlabel('Axial Height (mm)')
    plt.ylabel('Radius (mm)')
    plt.grid()

    return z,r

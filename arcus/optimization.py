import numpy as np

def lineDetection(wave,Fl,Fc,A,R,dE=100.,S=5.):
    """Compute observation time required to detect a line
    at a signal-to-noise ratio S.
    Inputs are:
    S - required signal-to-noise (default 5)
    wave - line wavelength in angstroms
    Fl - line flux in ph/cm2/sec
    Fc - continuum+bg flux in ph/cm2/sec/angstrom
    A - vector of effective areas in cm2 for relevant orders
    R - vector of resolutions (lambda/fwhm) for relevant orders
    dE - intrinsic line width (likely thermal width) in km/sec
    """
    #Convert intrinsic width to delta wavelength in angstroms
    dE = dE*1000./3.e8 * wave
    print 'Intrinsic Width: %.2e' % dE
    #Determine line widths for each order
    linewidth = wave/R
    print 'Instrumental Width: %.2e' % linewidth
    #Add intrinsic width to instrumental resolution in quadrature
    linewidth = np.sqrt(linewidth**2 + dE**2)
    #Adjust line flux based on Gaussian FWHM approximation
    Fl = Fl * .76
    #Create vector of denominators in merit function
    denom = Fl*A/(1+linewidth*(Fc/Fl))
    #Compute merit function
    merit = S**2/np.sum(denom)
    #Number of line counts
    nl = np.sum(Fl*A*merit)
    nc = np.sum(Fc*A*linewidth*merit)
    print 'Line counts: ' + str(nl)
    print 'Continuum counts: ' + str(nc)

    return merit

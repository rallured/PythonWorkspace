import numpy as np

#This module contains Fourier analysis routine

def components(d,win=1):
    """Want to return Fourier components with optional window
    Application note: These components are dependent on sampling!
    This means you can *not* interpolate these components onto other
    frequency grids!
    """
    #Handle window
    if win is not 1:
        if np.size(np.shape(d)) is 1:
            win = win(np.size(d))/np.sqrt(np.mean(win(np.size(d))**2))
        else:
            win1 = win(np.shape(d)[0])
            win2 = win(np.shape(d)[1])
            win = np.outer(win1,win2)
            win = win/np.sqrt(np.mean(win**2))

    #Compute Fourier components
    return np.fft.fftn(d*win)/np.size(d)

def continuousComponents(d,dx,win=1):
    """Want to return Fourier components with optional window
    Divide by frequency interval to convert to continuous FFT
    These components can be safely interpolated onto other frequency
    grids. Multiply by new frequency interval to get to numpy format
    FFT. Frequency units *must* be the same in each case.
    """
    #Handle window
    if win is not 1:
        if np.size(np.shape(d)) is 1:
            win = win(np.size(d))/np.sqrt(np.mean(win(np.size(d))**2))
        else:
            win1 = win(np.shape(d)[0])
            win2 = win(np.shape(d)[1])
            win = np.outer(win1,win2)
            win = win/np.sqrt(np.mean(win**2))

    #Compute Fourier components
    return np.fft.fftn(d*win)*dx**2
    

def freqgrid(d,dx=1.):
    """Return a frequency grid to match FFT components
    """
    freqx = np.fft.fftfreq(np.shape(d)[1],d=dx)
    freqy = np.fft.fftfreq(np.shape(d)[0],d=dx)
    freqx,freqy = np.meshgrid(freqx,freqy)
    return freqx,freqy

def ellipsoidalHighFrequencyCutoff(d,fxmax,fymax,dx=1.,win=1):
    """A simple low-pass filter with a high frequency cutoff.
    The cutoff boundary is an ellipsoid in frequency space.
    All frequency components with (fx/fxmax)**2+(fy/fymax)**2 > 1.
    are eliminated.
    fxmax refers to the second index, fymax refers to the first index
    This is consistent with indices in imshow
    """
    #FFT components in numpy format
    fftcomp = components(d,win=win)*np.size(d)

    #Get frequencies
    freqx,freqy = freqgrid(d,dx=dx)

    #Get indices of frequencies violating cutoff
    ind = (freqx/fxmax)**2+(freqy/fymax)**2 > 1.
    fftcomp[ind] = 0.

    #Invert the FFT and return the filtered image
    return fft.ifftn(fftcomp)

def realPSD(d,win=1,dx=1.):
    """This function returns the PSD of a real function
    Gets rid of zero frequency and puts all power in positive frequencies
    Returns only positive frequencies
    """
    #Get Fourier components
    c = components(d,win=win)

    #Reform into PSD
    if np.size(np.shape(c)) is 2:
        f = [np.fft.fftfreq(np.shape(c)[0],d=dx)[:np.shape(c)[0]/2],\
                   np.fft.fftfreq(np.shape(c)[1],d=dx)[:np.shape(c)[1]/2]]
        c = c[:np.shape(c)[0]/2,:np.shape(c)[1]/2]
        c[0,0] = 0.
        #Handle normalization
        c = 2*c
        c[0,:] = c[0,:]/np.sqrt(2.)
        c[:,0] = c[:,0]/np.sqrt(2.)
        
    elif np.size(np.shape(c)) is 1:
        f = np.fft.fftfreq(np.size(c),d=dx)
        f = f[:np.size(c)/2]
        c = c[:np.size(c)/2]
        c[0] = 0.
        c = c*np.sqrt(2.)

    return f,np.abs(c)**2

def lowpass(d,dx,fcut):
    """Apply a low pass filter to a 1 or 2 dimensional array.
    Supply the bin size and the cutoff frequency in the same units.
    """
    #Get shape of array
    sh = np.shape(d)
    #Take FFT and form frequency arrays
    f = np.fft.fftn(d)
    fx = np.fft.fftfreq(sh[0],d=dx)
    fy = np.fft.fftfreq(sh[1],d=dx)
    fa = np.meshgrid(fy,fx)
    fr = np.sqrt(fa[0]**2+fa[1]**2)
    #Apply cutoff
    f[fr>fcut] = 0.
    #Inverse FFT
    filtered = np.fft.ifftn(f)
    return filtered


                     

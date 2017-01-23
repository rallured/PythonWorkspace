import numpy as np
import matplotlib.pyplot as plt
import pdb

def computeResolution(length,dx,period,N,sh=np.pi/5,pad=0,win=None,\
                      phaseerr=None,xr=10.):
    """
    Compute the FFT of a 1D grating. Determine resolution f/df.
    Zero padding might need to be used in order to resolve the
    frequency peak.
    """
    #Construct perfect grating
    x = np.arange(0.,length,dx)
##    y1 = np.cos(2*np.pi*x/period)

    #Construct phase function
    if phaseerr is None:
        phaseerr = np.random.uniform(high=sh,size=N)
    phaseerr = np.tile(phaseerr,(np.ceil(len(x)/N),1))
    phaseerr = np.transpose(phaseerr).flatten()
    phaseerr = phaseerr[:len(x)]

    #Construct error grating
##    y2 = np.cos(2*np.pi*x/period+phaseerr)

    #Construct phase error function
    phaseerr = np.exp(1j*phaseerr)
    const = np.repeat(1.,len(phaseerr))

    #Take the FFTs
    if pad < len(x):
        pad=len(x)
    pad = int(pad)

##    if win is not None:
##        win = win(len(x))
##    else:
##        win = 1.
        
    freq = np.fft.fftfreq(pad,d=dx)
##    ff1 = np.fft.fft(y1,n=pad)
##    ff2 = np.fft.fft(y2,n=pad)
    ff3 = np.fft.fft(phaseerr,n=pad)
    ff4 = np.fft.fft(const,n=pad)

    #Convert to intensity and shift
    freq = np.fft.fftshift(freq)
    ff3 = np.abs(np.fft.fftshift(ff3))**2
    ff4 = np.abs(np.fft.fftshift(ff4))**2
    ff3,ff4 = ff3/np.sum(ff3),ff4/np.sum(ff4)

    #Compute CDFs
    cdf3 = np.cumsum(ff3)
    cdf4 = np.cumsum(ff4)
    hew3 = freq[np.argmin(np.abs(cdf3-.75))] - \
           freq[np.argmin(np.abs(cdf3-.25))]
    hew4 = freq[np.argmin(np.abs(cdf4-.75))] - \
           freq[np.argmin(np.abs(cdf4-.25))]

    ind = abs(freq)<xr
    
    return freq[ind],ff3[ind],ff4[ind],hew3,hew4

def varyingPeriod(length,dx,N,freq=None,pad=0):
    #Construct perfect grating
    x = np.arange(0.,length,dx)
##    y1 = np.cos(2*np.pi*x/period)

    #Construct phase function
    if freq is None:
        freq = np.random.uniform(low=6000.,high=6300.,size=N)
    freq = np.tile(freq,(np.ceil(len(x)/N),1))
    freq = np.transpose(freq).flatten()
    freq = freq[:len(x)]

    #Construct profiles
    y1 = np.cos(2*np.pi*x*freq)
    y2 = np.cos(2*np.pi*x*np.mean(freq))

    #Take the FFTs
    if pad < len(x):
        pad=len(x)
    pad = int(pad)

    ff3 = np.fft.fft(y1,n=pad)
    ff4 = np.fft.fft(y2,n=pad)

    freq = np.fft.fftfreq(pad,d=dx)
    freq = np.fft.fftshift(freq)

    ff3 = np.abs(np.fft.fftshift(ff3))**2
    ff4 = np.abs(np.fft.fftshift(ff4))**2
    ff3,ff4 = ff3/np.sum(ff3),ff4/np.sum(ff4)

    return freq,ff3,ff4,x,y1,y2

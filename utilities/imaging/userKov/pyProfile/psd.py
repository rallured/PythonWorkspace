from profile import line
from pyGeneralRoutines import span

def psd(x,y,retall=False):
    """return frequencies and PSD of a profile in mm^3, with x in mm and y in um.
    If retall is set, also the subtracted leveling profile
    PSD is squared FFT divided by step size. Line through ends is removed. Profile is assumed real, so only positive freqs are doubled and used. Return value has ceil(N/2)
    elements (always odd)."""
    N = len(y)
    L=span(x,True)
    x=y-line(x,y)-y.mean()
    yfft  = np.fft.fft(y)
    # let's take only the positive frequencies and normalize the amplitude
    freqs = np.fft.fftfreq(N,L/N)
    freqs = freqs[:np.ceil(N/2)]
    yfft  = yfft[:np.ceil(N/2)]
    psd  = np.abs(yfft)**2/N*L*2*1e-6
    if retall:
        return freqs,psd,line,np.angle(yfft)
    else:
        return freqs,psd

def profile_from_psd(f,psd,phase=None):
    """build a profile from PSD, if phase (in rad) is not passed use random values."""
    
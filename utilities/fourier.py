from numpy import *

#This module contains Fourier analysis routine


#Want to return Fourier components with optional window
def components(d,win=1):
    #Handle window
    if win is not 1:
        if size(shape(d)) is 1:
            win = win(size(d))/sqrt(mean(win(size(d))**2))
        else:
            win1 = win(shape(d)[0])
            win2 = win(shape(d)[1])
            win = outer(win1,win2)
            win = win/sqrt(mean(win**2))

    #Compute Fourier components
    return fft.fftn(d*win)/size(d)

#This function returns the PSD of a real function
#Gets rid of zero frequency and puts all power in positive frequencies
#Returns only positive frequencies
def realPSD(d,win=1,dx=1.):
    #Get Fourier components
    c = components(d,win=win)

    #Reform into PSD
    if size(shape(c)) is 2:
        f = [fft.fftfreq(shape(c)[0],d=dx)[:shape(c)[0]/2],\
                   fft.fftfreq(shape(c)[1],d=dx)[:shape(c)[1]/2]]
        c = c[:shape(c)[0]/2,:shape(c)[1]/2]
        c[0,0] = 0.
    elif size(shape(c)) is 1:
        f = fft.fftfreq(size(c),d=dx)
        f = f[:size(c)/2]
        c = c[:size(c)/2]
        c[0] = 0.

    return f,abs(c*2)**2

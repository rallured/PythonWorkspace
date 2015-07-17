from pyGeneralRoutines.span import span
import os

def line(x,y):
    """return line through end points of x,y."""
    L=span(x,size=1)
    return (x-x[0])*(y[-1]-y[0])/L+y[0]

def movingaverage(values,window):
    weigths = np.repeat(1.0, window)/window
    #including valid will REQUIRE there to be enough datapoints.
    #for example, if you take out valid, it will start @ point one,
    #not having any prior points, so itll be 1+0+0 = 1 /3 = .3333
    smas = np.convolve(values, weigths, 'same')
    return smas # as a numpy array

'''
def autotilt(x,y):
    """Transform a profile by tilting in a way that two contact points are on same horizontal line, return removed line."""
    if delta[0]==0 and delta[-1]==0
    L=span(x,size=1)
    line=x*(delta[-1]+delta[0])/L-delta[0]
    delta=delta-line  #start from endpoints
    i,j=delta.argsort()[:2]   #index of two maxima
'''

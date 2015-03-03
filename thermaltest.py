from numpy import *

def angunc(x,dx,h,dh,t,dt):
    tterm = x/h+t

    dpdt = .5*(1/(1+tterm**2)-1)
    dpdh = .5*(1/(1+tterm**2))*(x/h**2)
    dpdx = .5*(1/(1+tterm**2))*(1/h)

    return (dpdt*dt)**2+(dpdh*dh)**2+(dpdx*dx)**2

def disp(t,dphi,h):
    return h*(tan(t+2*dphi)-tan(t))

import numpy as np
from numpy import pi,sqrt,sin,cos,tan,exp
import pdb

def blazeYaw(inc,wave,m,d):
    """Determine required yaw angle for blaze wavelength
    at order with groove period d after setting incidence
    angle inc."""
    return np.arcsin(m*wave/2/d/cos(inc))

def blazeAngle(inc,wave,m,d):
    """Determine the blaze angle for wavelength wave at order m
    with an incidence angle inc and groove period d
    """
    psi = blazeYaw(inc,wave,m,d)
    beta1 = cos(inc)*cos(psi)
    alpha1 = cos(inc)*sin(psi)-m*wave/d
    return np.arcsin(alpha1/cos(np.arcsin(beta1)))

def blazeAngle2(inc,wave,m,d):
    a = m*wave/2/d
    return np.arcsin(a/sqrt(1-cos(inc)**2+a**2))
    
def eta(phi0,theta0):
    return np.arcsin(np.cos(phi0)*np.cos(theta0))

def yaw(phi0,theta0):
    return np.arctan(np.sin(theta0)/np.tan(phi0))

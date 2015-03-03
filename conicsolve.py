from numpy import *
from matplotlib.pyplot import *
import pdb

#Return radius of mirror at arbitrary z coordinate
def primrad(z,r0,z0):
    alpha = .25*arctan(r0/z0)
    thetah = 3*alpha
    thetap = alpha
    p = z0*tan(4*alpha)*tan(thetap)
    d = z0*tan(4*alpha)*tan(4*alpha-thetah)
    e = cos(4*alpha)*(1+tan(4*alpha)*tan(thetah))

    return sqrt(p**2+2*p*z+(4*e**2*p*d)/(e**2-1))

def secrad(z,r0,z0):
    alpha = .25*arctan(r0/z0)
    thetah = 3*alpha
    thetap = alpha
    p = z0*tan(4*alpha)*tan(thetap)
    d = z0*tan(4*alpha)*tan(4*alpha-thetah)
    e = cos(4*alpha)*(1+tan(4*alpha)*tan(thetah))

    return sqrt(e**2*(d+z)**2-z**2)

#Wolter parameters
def woltparam(r0,z0):
    alpha = .25*arctan(r0/z0)
    thetah = 3*alpha
    thetap = alpha
    p = z0*tan(4*alpha)*tan(thetap)
    d = z0*tan(4*alpha)*tan(4*alpha-thetah)
    e = cos(4*alpha)*(1+tan(4*alpha)*tan(thetah))

    return (alpha,p,d,e)

#Return distance to primary focus
def primfocus(r0,z0):
    alpha,p,d,e = woltparam(r0,z0)
    return z0 + 2*e**2*d/(e**2-1)

#Test Mathematica raytrace output
def mathraytrace(r0,z0,r):
    alpha,p,d,e = woltparam(r0,z0)
    a = p**2 + 4 * e**2 * p * d / (e**2 - 1)
    b = 2 * p
    l = e**2 * d**2
    m = 2 * e**2 * d
    n = e**2 - 1

    z1 = (r**2 - a) / b
    x2 = (b**2 *(2*b*m - 4*a*n + b**2*n)*r - 4*(2*b*m - 4*a*n + b**2*n)*r**3 +\
        2*b*r * sqrt(b**4*(m**2 - 4*l*n) + \
        4*(2*b*(-8*a*m + b*(8*l + (2*b - m)*m)) + ((-4*a + b**2)**2 + \
        8*b**2*l)*n)*r**2 + 16*(m**2 - 4*l*n)*r**4))/(b**4*n - \
       8*b**2*(2 + n)*r**2 + 16*n*r**4) #Mathematica solve for x2
    z2 = (-m+sqrt(m**2-4*n*(l-x2**2)))/2/n #Quad equation solve for z2
    pdb.set_trace()

#What r0 to achieve rgoal at zmax?
def rGoal_to_rMax(rgoal,z0,zmax):
    rguess = linspace(rgoal-2.,rgoal,10000)
    rnext = primrad(zmax,rguess,z0)
    return rguess[argmin(abs(rgoal-rnext))]

#Determine set of primary prescriptions to intercept beam
def primaryintercept(rmax,rmin,z0,zmin,zmax):
    rnew = rGoal_to_rMax(rmax,z0,zmax)
    print rnew
    while rnew > rmin:
        rmax = primrad(zmin,rnew,z0)
        rnew = rGoal_to_rMax(rmax,z0,zmax)
        print rnew

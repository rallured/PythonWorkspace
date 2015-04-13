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

#W-S Parameters
def wsPrimFunction(r0,z0,psi,x,y,z):
    a,p,d,e = woltparam(r0,z0)
    betas = 4*a
    ff = z0/np.cos(betas)
    g = ff / psi
    k = tan(betas/2)**2
    pdb.set_trace()

    beta = arcsin(sqrt(x**2+y**2)/ff)
    ind = beta<betas
    beta[ind] = betas
    kterm= (1/k)*tan(beta/2)**2 -1
    F = -z - ff*sin(betas/2)**2 + \
        ff**2*sin(beta)**2/(4*ff*sin(betas/2)**2) + \
        g*cos(beta/2)**4*kterm**(1-k)
    Fb = ff**2*sin(beta)*cos(beta)/(2*ff*sin(betas/2)**2) - \
           2*g*cos(beta/2)**3*sin(beta/2)*(kterm)**(1-k) + \
           g*(1-k)*cos(beta/2)*sin(beta/2)*(kterm)**(-k)*(1/k)
    dbdx = x/sqrt(1-(x**2+y**2)/ff**2)/ff/sqrt(x**2+y**2)
    dbdy = y/sqrt(1-(x**2+y**2)/ff**2)/ff/sqrt(x**2+y**2)
    Fx = Fb * dbdx
    Fy = Fb * dbdy
    Fx[ind] = 0.
    Fy[ind] = 0.
    Fz = -1
    return F, Fx, Fy

def wsSecFunction(r0,z0,psi,x,y,z):
    """Changed
    """
    a,p,d,e = woltparam(r0,z0)
    betas = 4*a
    ff = z0/np.cos(betas)
    g = ff / psi
    k = tan(betas/2)**2
    pdb.set_trace()

    beta = arctan2(sqrt(x**2+y**2),z)
    ind = beta<betas
    beta[ind] = betas
    kterm= (1/k)*tan(beta/2)**2 -1
    pdb.set_trace()
    a = (1-cos(beta))/(1-cos(betas))/ff + \
        (1+cos(beta))/(2*g)*kterm**(1+k)
    F = -z + cos(beta)/a
    #Add correction term to beta<betas indices
    dadbs = sin(betas)/ff/(1-cos(betas))+\
            (k+1)*(cos(betas)+1)*tan(betas/2)/cos(betas/2)**2/2/g/k
    dbdzs = -sin(betas)**2/sqrt(x[ind]**2+y[ind]**2)
    gam = (-ff*sin(betas)-ff**2*cos(betas)*dadbs)*dbdzs
    F[ind] = F[ind] + gam*(z[ind]-sqrt(x[ind]**2+y[ind]**2)/tan(betas))
    dadb = sin(beta)/ff/(1-cos(betas)) - \
           sin(beta)/(2*g)*(kterm)**(1+k) + \
           (k+1)*(cos(beta)+1)*tan(beta/2)*(kterm**k)/2/g/k/(cos(beta/2)**2)
    Fb = -sin(beta)/a - cos(beta)/a**2*dadb
##    dbdx = x/z/sqrt(x**2+y**2)/sqrt(1-(x**2+y**2)/z**2)
##    dbdy = y/z/sqrt(x**2+y**2)/sqrt(1-(x**2+y**2)/z**2)
##    dbdz = -sqrt(x**2+y**2)/z**2/sqrt(1-(x**2+y**2)/z**2)
    dbdx = x*z/(x**2+y**2+z**2)/sqrt(x**2+y**2)
    dbdy = y*z/(x**2+y**2+z**2)/sqrt(x**2+y**2)
    dbdz = -sqrt(x**2+y**2)/(x**2+y**2+z**2)
    Fx = Fb * dbdx
    Fy = Fb * dbdy
    Fz = -1. + Fb*dbdz
    Fx[ind] = -2./tan(betas)*x[ind]/sqrt(x[ind]**2+y[ind]**2)
    Fy[ind] = -2./tan(betas)*y[ind]/sqrt(x[ind]**2+y[ind]**2)
    Fz[ind] = gam - 1.
    return F, Fx, Fy, Fz

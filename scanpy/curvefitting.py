##Curve fitting routines
import numpy as np
import scipy as sp
import scipy.optimize
#Par = [amplitude,center,sigma,const] So you know what the different
#fit parameters are

#This is just the gaussian function.  You give it the independent
#variable x and the function parameters and it returns the dependent
#gaussian values.  If nterms=4, it returns a guassian plus a constant
def gaussian(par,x):
    if len(par)==3:
        return par[0]*np.exp((-(x-par[1])**2)/(2*par[2]**2))
    else:
        return par[3]+par[0]*np.exp((-(x-par[1])**2)/(2*par[2]**2))

#This takes the independent variable x, the measured values y,
#and the fit parameters and figures out what your residuals are
def gausserrfunct(par,x,y):
    return y - gaussian(par,x)

#This computes the Jacobian.  It only uses par and x, y is given
#as an input due to the way optimize.leastsq works.
def gaussderiv(par,x,y):
    gauss = gaussian(par,x) #Returns the value of the gaussian
    d0 = gauss/par[0] #Returns the partial with respect to par[0]
    d1 = gauss*(x-par[1])/par[2]**2 #Returns the partial with respect to par[1]
    d2 = gauss*((x-par[1])**2)/(par[2]**3) #Returns the partial with respect to par[2]
    if len(par)==3:
        return np.array([-d0,-d1,-d2])
    else:
        return np.array([-d0,-d1,-d2,-np.ones(len(x))]) #Partial with respect to constant is 1

#Estimates the fit parameters based on the independent variable x and the
#measured values y
def gaussguess(x,y,nterms=3):
    yguess = y
    if nterms==4:
        const = np.mean([y[0:4],y[len(y)-4:len(y)]]) #Estimates value of constant
        yguess = y-const #Subtracts constant off to get the rest of the fit parameters
    x0 = x[np.where(yguess==max(yguess))] #Estimates central x value
    a = np.array(max(yguess)) #Estimates amplitude
    diff = np.abs(yguess-(max(yguess)/2)) #Returns vector of differences from max
    fwhm = np.abs(x0-x[np.where(diff == min(diff))])*2 #Calculates FWHM
    sigma = fwhm/(2*np.sqrt(2*np.log(2))) #Calculates sigma from FWHM
    if nterms==4:
        return [a,x0,sigma,const]
    else:
        return [a,x0,sigma]

#This function does the actual fitting
def gaussfit(x,y,nterms=3):
    #Find initial guess
    par = gaussguess(x,y,nterms=nterms)
    #Calls fitting routine - passes args to gaussderiv and gausserrfunct - ftol and xtol should be played with if the fits are not good
    return sp.optimize.leastsq(gausserrfunct,par,Dfun=gaussderiv,col_deriv=1,ftol=1e-10,xtol=1e-10,args=(x,y))[0] #0th element is the fit parameters

from numpy import pi, inf, sqrt, cos, exp, log, array, arange, isscalar, asarray
import pyfits
import matplotlib.pyplot as plt
from scipy.optimize import leastsq
#from scipy.optimize import curve_fit


def _general_function(params, xdata, ydata, function):
    return function(xdata, *params) - ydata

def _weighted_general_function(params, xdata, ydata, function, weights):
    return weights * (function(xdata, *params) - ydata)

def curve_fit(f, xdata, ydata, p0=None, sigma=None, **kw):
    """
    Use non-linear least squares to fit a function, f, to data.

    Assumes ``ydata = f(xdata, *params) + eps``

    Parameters
    ----------
    f : callable
        The model function, f(x, ...).  It must take the independent
        variable as the first argument and the parameters to fit as
        separate remaining arguments.
    xdata : An N-length sequence or an (k,N)-shaped array
        for functions with k predictors.
        The independent variable where the data is measured.
    ydata : N-length sequence
        The dependent data --- nominally f(xdata, ...)
    p0 : None, scalar, or M-length sequence
        Initial guess for the parameters.  If None, then the initial
        values will all be 1 (if the number of parameters for the function
        can be determined using introspection, otherwise a ValueError
        is raised).
    sigma : None or N-length sequence
        If not None, it represents the standard-deviation of ydata.
        This vector, if given, will be used as weights in the
        least-squares problem.


    Returns
    -------
    popt : array
        Optimal values for the parameters so that the sum of the squared error
        of ``f(xdata, *popt) - ydata`` is minimized
    pcov : 2d array
        The estimated covariance of popt.  The diagonals provide the variance
        of the parameter estimate.

    Notes
    -----
    The algorithm uses the Levenburg-Marquardt algorithm:
    scipy.optimize.leastsq. Additional keyword arguments are passed directly
    to that algorithm.

    Examples
    --------
    >>> import numpy as np
    >>> from scipy.optimize import curve_fit
    >>> def func(x, a, b, c):
    ...     return a*np.exp(-b*x) + c

    >>> x = np.linspace(0,4,50)
    >>> y = func(x, 2.5, 1.3, 0.5)
    >>> yn = y + 0.2*np.random.normal(size=len(x))

    >>> popt, pcov = curve_fit(func, x, yn)

    """
    if p0 is None or isscalar(p0):
        # determine number of parameters by inspecting the function
        import inspect
        args, varargs, varkw, defaults = inspect.getargspec(f)
        if len(args) < 2:
            msg = "Unable to determine number of fit parameters."
            raise ValueError(msg)
        if p0 is None:
            p0 = 1.0
        p0 = [p0]*(len(args)-1)

    args = (xdata, ydata, f)
    if sigma is None:
        func = _general_function
    else:
        func = _weighted_general_function
        args += (1.0/asarray(sigma),)
    res = leastsq(func, p0, args=args, full_output=1, **kw)
    (popt, pcov, infodict, errmsg, ier) = res

    if ier != 1:
        #msg = "Optimal parameters not found: " + errmsg
        #raise RuntimeError(msg)
        print "Optimal parameters not found: " + errmsg

    if (len(ydata) > len(p0)) and pcov is not None:
        s_sq = (func(popt, *args)**2).sum()/(len(ydata)-len(p0))
        pcov = pcov * s_sq
    else:
        pcov = inf

    return popt, pcov


def wavg(value, err):
	a = sum(value/err)
	b =  sum(1/err)
	return (a/b)

def roi_stats(x, c, roi):
    """Returns statistics of a region of interest in the counts histogram
    """
    q = (x >= roi[0]) & (x <= roi[1])
    total = sum(c[q])
    centroid = sum(x[q]*c[q])/total
    sigma = sqrt(sum(c[q]*(x[q]-centroid)**2)/total)
    maximum = max(c[q])
    return [total, centroid, sigma, maximum]
 
def get_rate(filename, roi, plot=False):
# read a fits file and get the count rate in the region of interest
    h = pyfits.open(filename)
    rotation = h[1].header['ROTATION']
    accumtime = h[1].header['EXPOSURE']
    cols = h[1].columns
    tbdata = h[1].data # data is in first extension
    chan = tbdata.field('CHANNEL')
    counts = tbdata.field('COUNTS')
    stats = roi_stats(chan, counts, roi)
    rate = stats[0]/accumtime
    rate_err = sqrt(stats[0])/accumtime
    centroid = stats[1]
    centroid_err = stats[2]/sqrt(stats[0])
    # plot = False
    if plot:
        # print "Plotting MCA data"
        plt.clf()
        plt.plot(chan, counts)
        plt.xlabel('Channel')
        plt.ylabel('Number events/channel')
        plt.show() # draw() #
        for i in range(2):
            plt.plot([roi[i], roi[i]], [0, 0.2*stats[3]])
        plt.plot([stats[1], stats[1]], [0, stats[3]])
        plt.draw()
        plt.show()
    # print filename, rotation, rate, rate_err
    return rotation, rate, rate_err, centroid, centroid_err

def func_cos2t(x, a0, a2, phi2):
    return a0*(1 + exp(a2)*cos(2*(x-phi2)*pi/180))
fitfunc = lambda p, x: p[0]*(1 + exp(p[1])*cos(2*(x-p[2])*pi/180))

def func_2cos(x, a0, a2, phi2, a1, phi1):
    return a0*(1 + exp(a2)*cos(2*(x-phi2)*pi/180.0) + exp(a1)*cos((x-phi1)*pi/180.0))
fitfunc2 = lambda p, x: p[0]*(1 + exp(p[1])*cos(2*(x-p[2])*pi/180) + exp(p[3])*cos((x-p[4])*pi/180))


roi = array([238, 270])

#filebase, angles, roi = "100406/100406a", arange(0.0, 225+1, 15), array([230, 280])
#filebase, angles = "100402/100402a", arange(15.0, 225+1, 15)
#filebase = "100402b"

#filebase, angles, roi = "100412/100412a", arange(0.0, 360+1, 10), array([238, 273])
#angles = angles[angles <> 170]
#filebase, angles, roi = "100412b", arange(0.0, 360+1, 10), array([238, 273])
#angles = angles[angles >= 180]
#filebase, angles, roi = "100412d", arange(0.0, 350+1, 10), array([238, 273])
#filebase, angles, roi = "100412e", arange(0.0, 360+1, 10), array([238, 273])

#datadir = '100413/'
#filebase, angles, roi = "100413c", arange(0.0, 360+1, 10), array([237, 273])
#angles = angles[angles <> 60] # run at 60 degrees has very high rate
#filebase, angles, roi = "100413d", arange(0.0, 360+1, 10), array([237, 273])

#datadir = '100414/'
filebase, angles, roi = "100414b", arange(0.0, 360+1, 10), array([237, 273])
#filebase, angles, roi = "100414d", arange(0.0, 360+1, 10), array([237, 273])
#angles = angles[(angles > 25) & (angles <> 200)]
#filebase, angles, roi = "100414e", arange(0.0, 360+1, 10), array([237, 273])

datadir = ''
#filebase, angles, roi = "100415a", arange(0.0, 360+1, 10), array([237, 273])
filebase, angles, roi = "100415a", arange(0.0, 360+1, 10), array([460, 540])
#filebase, angles, roi = "100415a", arange(0.0, 360+1, 10), array([235, 274])
#angles = angles[(angles <> 10)& (angles <> 30)& (angles <> 120)& (angles <> 350)]
#filebase, angles, roi = "100416/100416a", np.arange(0.0, 360+1, 10), np.array([235, 278])
#angles = angles[(angles <> 30)& (angles <> 90)& (angles <> 150)& (angles <> 180)& (angles <> 320)& (angles <> 350)]
#filebase, angles, roi = "100416/100416b", np.arange(0.0, 360+1, 10), np.array([235, 278])
#angles = angles[(angles <> 50)& (angles <> 60)& (angles <> 200)& (angles <> 230)]
#filebase, angles, roi = "100416/100416c", arange(0.0, 360+1, 15), array([235, 278])
#angles = angles[(angles <> 195)& (angles <> 225)]
#filebase, angles, roi = "100419/100419a", arange(0.0, 360+1, 15), array([235, 278])
#angles = angles[(angles <> 195)& (angles <> 270)]
filebase, angles, roi = "100419/100419b", arange(0.0, 360+1, 15), array([235, 278])
angles = angles[(angles <> 45)& (angles <> 360)]
#roi = array([460, 544])
#filebase, angles, roi = "100419/100419c", arange(0.0, 360+1, 15), array([235, 278])
#angles = angles[(angles <> 75) & (angles <> 210) & (angles <> 255)]
#filebase, angles, roi = "100419/100419d", arange(0.0, 360+1, 10), array([235, 278])
#angles = angles[(angles <> 50)&(angles <> 100)&(angles <> 200)&(angles <> 270)&(angles <> 280)]
#filebase, angles, roi = "100419/100419e", arange(0.0, 360+1, 10), array([235, 278])
#angles = angles[(angles <> 200)]
#filebase, angles, roi = "100419/100419f", arange(0.0, 360+1, 10), array([235, 278])
#angles = angles[(angles <> 10)&(angles <> 40)&(angles <> 150)]
#filebase, angles, roi = "100419/100419g", arange(0.0, 360+1, 10), array([235, 278])
#angles = angles[(angles <> 40)&(angles <> 50)&(angles <> 160)]
#filebase, angles, roi = "100419/100419h", arange(0.0, 360+1, 15), array([235, 278])
#filebase, angles, roi = "100420/100420a", arange(0.0, 360+1, 15), array([235, 278])
#angles = angles[(angles <> 30)&(angles <> 285)]
#filebase, angles, roi = "100420/100420b", arange(0.0, 360+1, 10), array([235, 278])
#angles = angles[(angles <> 200)]
#roi = array([460, 544])
#filebase, angles, roi = "100420/100420c", arange(0.0, 360+1, 10), array([235, 278])
#angles = angles[(angles <> 270)&(angles <> 320)&(angles <> 330)]
#filebase, angles, roi = "100420/100420e", arange(0.0, 360+1, 10), array([235, 278])
#angles = angles[(angles <> 70)&(angles <> 160)&(angles <> 270)]
#filebase, angles, roi = "100420/100420f", arange(0.0, 360+1, 10), array([235, 278])
#angles = angles[(angles <> 0)&(angles <> 170)&(angles <> 190)]
#filebase, angles, roi = "100421/100421a", arange(0.0, 360+1, 15), array([230, 286])
#angles = angles[(angles <> 15)&(angles <> 30)&(angles <> 180)&(angles <> 330)]
#filebase, angles, roi = "100421/100421b", arange(0.0, 360+1, 15), array([230, 286])
#angles = angles[(angles <> 75)&(angles <> 105)&(angles <> 120)&(angles <> 240)&(angles <> 315)]
#filebase, angles, roi = "100421/100421c", arange(0.0, 360+1, 10), array([230, 286])
#angles = angles[(angles <> 50)&(angles <> 80)&(angles <> 110)&(angles <> 130)&(angles <> 330)]
#filebase, angles, roi = "100421/100421d", arange(0.0, 360+1, 10), array([230, 286])
#filebase, angles, roi = "100421/100421e", arange(0.0, 360+1, 10), array([230, 286])
filebase, angles, roi = "100421/100421f", arange(0.0, 360+1, 15), array([230, 286])
#roi = array([460, 550])
filebase, angles, roi = "100421/100421g", arange(0.0, 360+1, 15), array([230, 286])
filebase, angles, roi = "100421/100421h", arange(0.0, 360+1, 15), array([230, 286])
filebase, angles, roi = "100430/100430a", arange(0.0, 360+1, 15), array([230, 286])
filebase, angles, roi = "100430/100430b", arange(0.0, 360+1, 15), array([230, 286])
#roi = array([460, 560])
filebase, angles, roi = "100430/100430c", arange(0.0, 360+1, 15), array([230, 280])
filebase, angles, roi = "100430/100430d", arange(0.0, 360+1, 15), array([230, 280])
filebase, angles, roi = "100430/100430e", arange(0.0, 360+1, 15), array([230, 280])
filebase, angles, roi = "100430/100430f", arange(0.0, 360+1, 15), array([230, 285])
filebase, angles, roi = "100430/100430g", arange(0.0, 360+1, 15), array([230, 285])
filebase, angles, roi = "100503/100503a", arange(0.0, 360+1, 15), array([230, 285])
filebase, angles, roi = "100503/100503b", arange(0.0, 360+1, 15), array([230, 285])
filebase, angles, roi = "100503/100503c", arange(0.0, 360+1, 15), array([230, 285])
filebase, angles, roi = "100503/100503d", arange(0.0, 360+1, 10), array([230, 285])
filebase, angles, roi = "100503/100503e", arange(0.0, 360+1, 10), array([230, 285])
filebase, angles, roi = "100503/100503f", arange(0.0, 360+1, 10), array([230, 285])
filebase, angles, roi = "100503/100503g", arange(0.0, 360+1, 10), array([230, 285])
filebase, angles, roi = "100503/100504a", arange(0.0, 360+1, 10), array([230, 285])

filebase,angles, roi = "100611/100611b", arange(0.0, 360+1, 15), array([110, 145])
#roi = array([230,285])
#roi = array([470, 540])

print "Processing ", filebase, " with roi=", roi
rate = 0.0*angles
rate_err = 0.0*rate
cent = 0.0*rate
cent_err = 0.0*rate
for i in range(len(angles)):
    pos = angles[i]
    spos = "%03d" % pos
    filename = datadir+filebase+"_r"+spos+".pha"
    r = get_rate(filename, roi, plot=True)
    if r[0] <> pos:
        print "Rotation angle does not match filename"
        pos = r[0]
    rate[i] = r[1]
    rate_err[i] = r[2]
    cent[i] = r[3]
    cent_err[i] = r[4]
    print angles[i], rate[i], rate_err[i], cent[i], cent_err[i]
#plt.clf()
#plt.plot(angles, rate)
#plt.show()


# find average, remove high points
avg = wavg(rate, rate_err)
print 'threshold for removal = ', 1.3*avg

## remove points more than 30% above average
#q = (rate < 1.3*avg)
#angles = angles[q]
#rate = rate[q]
#rate_err = rate_err[q]
#cent = cent[q]
#cent_err = cent_err[q]

# first, fit with a constant
avg = wavg(rate, rate_err)
chisq0 = sum(((rate-avg)/rate_err)**2)
print 'Average (c/s) = %.2f' % (avg)
dof0 = len(angles)-1
print 'Chisq/Dof = %.1f/%d' % (chisq0, dof0)

# fit to a*(1 + b*cos(2*(phi-phi0))
pinit = [avg, log(0.01*avg), 45.0]
p,t = curve_fit(func_cos2t, angles, rate, p0=pinit, sigma=rate_err)
p2t = p
perr = array([t[0,0], t[1,1], t[2,2]])
# define our fitting functions
#fitfunc = lambda p, x: p[0] + p[1]*cos(2*(x-p[2])*3.14159/180.0)
#errfunc = lambda p, x, y, err: (y - fitfunc(p, x)) / err
#pinit = [250.0, 200.0, 90.0]
#out = leastsq(errfunc, pinit, args=(angles, rate, rate_err), full_output=1)
#p7 = out[0]
#t = out[1]
#perr7 = array([t[0,0], t[1,1], t[2,2]])
print 'Average (c/s) = %.2f +/- %.2f' % (p[0], perr[0])
#mod = 100.0*p[1]/p[0]
#moderr = 100.0*sqrt((perr[1]/p[0])**2+(perr[0]*p[1]/p[0])**2)
mod = 100.0*exp(p[1])
moderr = 100.0*exp(p[1])*perr[1]
print 'Modulation (percent) = %.2f +/- %.2f' % (mod, moderr)
ang = p[2]
angerr = perr[2]
print 'Angle (degrees) = %.1f +/- %.1f ' % (ang, angerr)
chisq = sum(((fitfunc(p, angles)-rate)/rate_err)**2)
dof = len(angles)-3
print 'Chisq/Dof = %.1f/%d' % (chisq, dof)

# fit to a*(1 + b*cos(2*(phi-phi0) + c*cos(phi-phi0)
#pinit = [avg, log(0.02*avg), 150.0, log(0.02*avg), 180.0]
pinit = [p[0], p[1], p[2], log(0.02*avg), 180.0]
print 'pinit=', pinit
p,t = curve_fit(func_2cos, angles, rate, p0=pinit, sigma=rate_err)
if len(t) == 1:
    perr = array([-1.0, -1.0, -1.0, -1.0, -1.0])
else:
    perr = array([t[0,0], t[1,1], t[2,2], t[3,3], t[4,4]])
avg2 = p[0]
avg2err = perr[0]    
print 'Constant (c/s) = %.2f +/- %.2f' % (avg2, avg2err)
mod2 = 100.0*exp(p[1])
mod2err = 100.0*exp(p[1])*perr[1]
print '2*phi Modulation (percent) = %.2f +/- %.2f' % (mod2, mod2err)
ang2 = p[2]
ang2err = perr[2]
print '2*phi Angle (degrees) = %.1f +/- %.1f ' % (ang2, ang2err)
mod1 = 100.0*exp(p[3])
mod1err = 100.0*exp(p[3])*perr[3]
print '2*phi Modulation (percent) = %.2f +/- %.2f' % (mod1, mod1err)
ang1 = p[4]
ang1err = perr[4]
print '2*phi Angle (degrees) = %.1f +/- %.1f ' % (ang1, ang1err)
chisq2 = sum(((fitfunc2(p, angles)-rate)/rate_err)**2)
dof2 = len(angles)-5
print 'Chisq/Dof = %.1f/%d' % (chisq2, dof2)

pang = arange(0.0, 360.0, 1.0)
plt.clf()
plt.errorbar(angles, rate, rate_err, fmt='o')
plt.plot(pang, fitfunc(p2t, pang), 'b--')
plt.plot(pang, fitfunc2(p, pang), 'b')
ylo = min(rate-rate_err)-0.3*(max(rate)-min(rate))
yhi= max(rate+rate_err)+0.2*(max(rate)-min(rate))
plt.axis([-5, 365, ylo, yhi])
plt.title(filebase+'  Varian VF-50J-Rh/S 71053-6W')
plt.xlabel('Rotation angle (degrees)')
plt.ylabel('Count rate (events/second)')
s1 = 'Average (c/s) = %.2f +/- %.2f' % (avg2, avg2err)
s2 = '2*phi mod (%s) = %.2f +/- %.2f' % ('%', mod2, mod2err)
s3 = 'Angle (deg) = %.1f +/- %.1f ' % (ang2, ang2err)
s4 = '1*phi mod (%s) = %.2f +/- %.2f' % ('%', mod1, mod1err)
s5 = 'Angle (deg) = %.1f +/- %.1f ' % (ang1, ang1err)
s6 = 'Chisq/Dof: c= %.1f/%d  c+cos= %.1f/%d  c+cos+cos= %.1f/%d' \
     % (chisq0, dof0, chisq, dof, chisq2, dof2)
plt.text(5, ylo, s1+'\n'+s2+' '+s3+'\n'+s4+' '+s5+'\n'+s6+'\n')

#clo = min(cent)
#chi= max(cent)
#plt.plot(angles, ylo+(yhi-ylo)*(cent-clo)/(chi-clo))
plt2 = plt.twinx()
cent_avg = sum(cent)/len(cent)
#plt2.plot(angles, cent-cent_avg, 'r')
plt2.errorbar(angles, 1000*(cent-cent_avg)/cent_avg, 1000*cent_err/cent_avg, fmt='+r')
plt2.set_ylabel(r'PHA $\Delta$/average ($\times$1000)')

plt.show()

""" code to read Amptek MCA files and fit peak with Gaussian
use fitmca to fit the data in one MCA file
 fitmca(datadir, datafile, plot=True)
   datadir = string with directory containing file, '' is ok for local directory
   datafile = string with name of data file
   plot=True for plot, =False for no plot

use fitall to fit the data in all the MCA files in a directory
 fitall(datadir, plotfile)
   datadir = string with directory containing file, '' is ok for local directory
   plotfile = name of pdffile with plots for all MCA files
            = 'Screen' plots to screen
            = 'None' no plots

"""

from numpy import pi, inf, sqrt, cos, exp, log, array, arange, isscalar, asarray
import time, os, pdb
import numpy as np
#import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.optimize import leastsq
from matplotlib.backends.backend_pdf import PdfPages

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
        msg = "Optimal parameters not found: " + errmsg
        raise RuntimeError(msg)

    if (len(ydata) > len(p0)) and pcov is not None:
        s_sq = (func(popt, *args)**2).sum()/(len(ydata)-len(p0))
        pcov = pcov * s_sq
    else:
        pcov = inf

    return popt, pcov

def readmca(filename):
# read data from an Amptek MCA file and return the
# fractional day of year for the start of the run,
# the live time in seconds,
# the MCA data histogram
  # Amptek MCA save file is in text format in multiple lines, read all lines
  f = open(filename, 'rU').readlines()
  # find lines beginning with '<<' these are tags
  t = [i for i, e in enumerate(f) if e.find('<<') == 0]
  # histogram of data lies between last two tags
  # should check that these lines are <<DATA>> and <<END>>
  mcahist = [float(f[i]) for i in range(t[-2]+1,t[-1]-1)]
  # write the Amptek header into fits header keywords
  for i in range(t[0]+1, t[1]-1):
    s = f[i]
    fs = f[i].find(' ')
    kname = s[0:fs]
    kvalue = s[fs+3:-1]
    if len(kvalue) > 0:
      if kname == 'LIVE_TIME':
        livetime = float(kvalue) ; 
      if kname == 'START_TIME': # has form 09/20/2011 18:16:45
        starttime = time.strptime(kvalue, '%m/%d/%Y %H:%M:%S')
        dayofyear = starttime.tm_yday + starttime.tm_hour/24.0 \
                    + starttime.tm_min/1440.0 + starttime.tm_yday/86400.0
  return dayofyear, livetime, mcahist

def roi_stats(x, c, roi):
    """Returns statistics of a region of interest in the counts histogram
    """
    q = (x >= roi[0]) & (x <= roi[1])
    total = sum(c[q])
    centroid = sum(x[q]*c[q])/total
    sigma = sqrt(sum(c[q]*(x[q]-centroid)**2)/total)
    return [total, centroid, sigma]
 
def func_gaussc(x, norm, cent, width, cons):
    return norm*exp(-(x-cent)**2/(2*width**2)) + cons
gaussc = lambda p, x: p[0]*exp(-(x-p[1])**2/(2*p[2]**2)) + p[3]

def func_gausszexp(x, norm, cent, width, nnorm, nwid):
    return norm*exp(-(x-cent)**2/(2*width**2)) + nnorm*exp(-x/nwid)
gausszexp = lambda p, x: p[0]*exp(-(x-p[1])**2/(2*p[2]**2)) + p[3]*exp(-x/p[4])

def func_twogauss(x, norm1, cent1, width1, norm2, cent2, width2, cons):
    return norm1*exp(-(x-cent1)**2/(2*width1**2)) + cons + \
           norm2*exp(-(x-cent2)**2/(2*width2**2))
twogauss = lambda p, x: p[0]*exp(-(x-p[1])**2/(2*p[2]**2)) + p[6] + \
           p[3]*exp(-(x-p[4])**2/(2*p[5]**2))

def func_specshape(x, a1, c1, sig1):
    a2 = .43262*a1
    c2 = .89461*c1
    sig2 = 2.9388754*sig1
    return twogauss([a1,c1,sig1,a2,c2,sig2,0],x)

def fitmca(datadir, filename, plot=False):
# read in an Amptek MCA data file then fit the data with a Gaussian
  # read in the data
  (dayofyear, livetime, mcahist) = readmca(datadir+filename)
  counts = array(mcahist) # convert data to numpy format
  nchan = len(mcahist) # number of channels in MCA data
  chan = arange(nchan) # array of channel numbers
  # find the mean and standard deviation of the data
  (t, c, s) = roi_stats(chan, counts, [0,255])
  #print t, c, s
  # fit a gaussian+constant to the data, use Gehrels variance function for error
  pinit = [t, c, s, 0.0] # use calculated statistics as first guess for fit
  p,t = curve_fit(func_gaussc, chan, counts, p0=pinit, sigma=1+sqrt(counts+0.75))
  #print 'fit1:', p
  # check the fit
  if(p[1] < p[2]*2.35): # if centroid < FWHM then 
    counts[0:s] = 0 # drop channels < stdev(data) and fit again
    (t, c, s) = roi_stats(chan, counts, [0,255])
    pinit = [t, c, s, 0.0] # use calculated statistics as first guess for fit
    p,t = curve_fit(func_gaussc, chan, counts, p0=pinit, sigma=1+sqrt(counts+0.75))
  # fit again using channels within +/- FWHM of centroid
  counts[0:p[1]-p[2]*2.35] = 0
  counts[p[1]+p[2]*2.35:255] = 0
  p,t = curve_fit(func_gaussc, chan, counts, p0=p, sigma=1+sqrt(counts+0.75))
  #print 'fit2:', p
  perr = array([t[0,0], t[1,1], t[2,2]])
  gcounts = p[0]*p[2]*sqrt(2*pi)
  gcountserr = sqrt(2*pi)*sqrt((perr[0]*p[2])**2+(p[0]*perr[2])**2)
  gcent = p[1]
  gcenterr = perr[1]
  gfwhm = p[2]*2.35
  gfwhmerr = perr[2]*2.35
  if plot:
    plt.clf()
    plt.plot(chan, mcahist)
    plt.xlabel('Channel')
    plt.ylabel('Number events/channel')
    plt.title(filename)
    plt.plot(chan, gaussc(p, chan))
    #plt.plot(chan, gausszexp(p, chan))
    s1 = 'Centroid = %.1f' % gcent
    s2 = 'FWHM = %.1f' % gfwhm
    s3 = 'FWHM/Cent = %.2f' % (gfwhm/gcent)
    s4 = 'Counts = %i' % gcounts
    plt.text(nchan*0.75, max(counts)*0.8, s1+'\n'+s2+'\n'+s3+'\n'+s4)
  return dayofyear, livetime, gcent, gcenterr, gfwhm, gfwhmerr, gcounts, gcountserr

def fitone(datadir, filename, plot=False):
# read in an Amptek MCA data file then fit the data with a Gaussian
  t = fitmca(datadir, filename, plot=plot)
  if plot:
    plt.show()
  return t

def fitall(datadir, plotfile):
# find all the Amptek MCA data files in a directory, fit each with a Gaussian,
# then print the results
  if plotfile <> 'Screen' and plotfile <> 'None' :
    pdffile = PdfPages(plotfile+'.pdf')
    print 'Will save plots to a pdf file: ', plotfile+'.pdf'
  filelist = [file for file in os.listdir(datadir) if file.lower().endswith('.mca')]
  print 'Files in: ', datadir
  print 'File name, centroid, FWHM, counts (in peak)'
  for f in filelist:
    r = fitmca(datadir, f, plot=True) # fit Amptek MCA file with Gaussian
    if plotfile == 'Screen' :
      plt.show()
      pdb.set_trace() # type c to continue
    else:
      if plotfile <> 'None' :
        pdffile.savefig()
    rs = ' %.1f, %.1f, %i'% (r[2], r[4], r[6])
    print f+','+rs
  if plotfile <> 'Screen' and plotfile <> 'None' :
    pdffile.close()

def fitallMXS(datadir, plotfile):
# find all the Amptek MCA data files in a directory, fit each with a Gaussian,
# then print the results
  pdffile = PdfPages(plotfile+'.pdf')
  print 'Will save plots to a pdf file: ', plotfile+'.pdf'
  filelist = [file for file in os.listdir(datadir) if file.lower().endswith('.mca')]
  print 'Files in: ', datadir
  print 'File name, centroid'
  plt.clf()
  for f in filelist:
    # read in MCA file, plot, and request ROI from user
    data = readmca(f)
    counts = np.array(data[2])
    chan = arange(np.size(counts))
    plt.plot(chan[0:100],counts[0:100])
    low = input('Lower bound? ')
    high = input('Upper bound? ')+1
    # guess parameters and fit
    ampguess = max(counts[low:high])
    centguess = np.mean(chan[np.where(counts==ampguess)])
    sigguess = sqrt(sum(counts[low:high]*(chan[low:high]-centguess)**2)/ \
        sum(counts[low:high]))*.5
    fit = curve_fit(func_specshape,chan[low:high],counts[low:high],\
          p0=[ampguess,centguess,sigguess],sigma=1+sqrt(counts[low:high]+.75))
    # Calc Error
    amp = fit[0][0]
    cen = fit[0][1]
    sig = fit[0][2]
    cenerr = fit[1][1,1]
    sigerr = fit[1][2,2]
    reserr = sqrt((2.35*sig/cen**2)**2*cenerr**2+(2.35/cen)**2*sigerr**2)
    # Plot results
    plt.plot(chan,func_specshape(chan,amp,cen,sig))
    plt.text(125,max(counts)*.8,'Line Centroid: '+str(cen))
    plt.text(125,max(counts)*.7,'Line FWHM: '+str(2.35*sig))
    plt.text(125,max(counts)*.6,'Resolution: '+\
        str(2.35*sig/cen))
    plt.text(125,max(counts)*.5,'Res Uncertainty: '+str(reserr))
    plt.title(f)
    plt.xlabel('Channel')
    plt.ylabel('Counts')
    pdffile.savefig()
    plt.clf()
    print f+','+str(fit[0][1])
  pdffile.close()

def test1():
  datadir = '/work/gems/thermal/GEMS_ETU2_Thermal_Test_Data/etu2_9_20_2011/'
  print fitone(datadir, '962_590_150_Fe55.mca', plot=True)

def testall():
  datadir = '/work/gems/thermal/GEMS_ETU2_Thermal_Test_Data/etu2_9_20_2011/'
  #print fitall(datadir, 'Screen')
  fitall(datadir, 'etu2_9_20_2011')

def thermplots(runsfile, plotmca=False):
# fit a Gaussian to the Amptek MCA data files listed in the file runsfile
# then plot the results
# runsfiles lists data files with ancillary information
# first line is data directory
# next lines are: "filename", coarse gain, pressure, temperature
# ended by line with END
  # read in the run list
  rf = open(runsfile+'.txt', 'r') # open the run list file
  datadir = rf.readline()[0:-1] # first line is data directory
  mcafiles, cgain, pressure, temperature = [], [], [], []
  a = rf.readline() # read next lines have: "filename", coarse gain, pressure, temperature
  while (a[0:3] <> 'END' and a <> ''): # stop if line says END
    t = a.split(",") # split string at commas
    mcafiles.append(t[0])
    cgain.append(float(t[1]))
    pressure.append(float(t[2]))
    temperature.append(float(t[3]))
    a = rf.readline()
  rf.close()  
  # out the plot output file
  pdffile = PdfPages(runsfile+'.pdf')
  print 'Will save plots to a pdf file: ', runsfile+'.pdf'
  # analyze all the MCA files
  print 'Files in: ', datadir
  atime, cent, centerr, fwhm, fwhmerr, counts = [], [], [], [], [], []
  print 'File name, centroid, FWHM, counts (in peak)'
  for i in range(len(mcafiles)) :
    r = fitmca(datadir, mcafiles[i], plot=plotmca) # fit Amptek MCA file with Gaussian
    if plotmca :
      pdffile.savefig() # save the MCA histogram plot
    atime.append(r[0])
    cent.append(r[2])
    centerr.append(r[3])
    fwhm.append(r[4])
    fwhmerr.append(r[5])
    counts.append(r[6])
    rs = ' %.1f, %.1f, %i'% (r[2], r[4], r[6])
    print mcafiles[i]+','+rs, atime[i]
  # calculate the gain
  cent = array(cent)
  centerr = array(centerr)
  fwhm = array(fwhm)
  fwhmerr = array(fwhmerr)
  cgain = array(cgain)
  atime = array(atime)
  gain = (cent*200/cgain*0.00024134+0.0013149)/1.602E-7/5900.0*26.0
  gainerr = gain*centerr/cent
  print temperature
  # make plots
  tzero, tzerostr = 263.0, "9/20/2011"
  plt.clf()
  plt.ylabel('GEM gain')
  plt.xlabel('Time (days since 9/20/2011)')
  y = gain
  ylo = min(y)-0.2*(max(y)-min(y))
  yhi = max(y)+0.2*(max(y)-min(y))
  plt.plot(atime-263, y, 'D', markersize=20)
  plt.axis([0, 1.1*max(atime-263), ylo, yhi])
  plt.show()
  pdffile.savefig()
  pdb.set_trace() # type c to continue
  plt.clf()
  plt.ylabel('GEM gain')
  plt.xlabel('Temperature (C)')
  y = gain
  ylo = min(y)-0.2*(max(y)-min(y))
  yhi = max(y)+0.2*(max(y)-min(y))
  plt.plot(temperature, y, 'D', markersize=20)
  plt.axis([15, 35, ylo, yhi])
  plt.show()
  pdffile.savefig()
  pdb.set_trace() # type c to continue
  plt.clf()
  plt.ylabel('GEM FWHM/Centroid')
  plt.xlabel('Time (days since '+tzerostr+')')
  y = fwhm/cent
  ylo = min(y)-0.2*(max(y)-min(y))
  yhi = max(y)+0.2*(max(y)-min(y))
  plt.plot(atime-263, fwhm/cent, 'D', markersize=20)
  plt.axis([0, 1.1*max(atime-263), ylo, yhi])
  plt.show()
  pdffile.savefig()
  pdb.set_trace() # type c to continue
  plt.clf()
  plt.ylabel('GEM FWHM/Centroid')
  plt.xlabel('Temperature (C)')
  y = fwhm/cent
  ylo = min(y)-0.2*(max(y)-min(y))
  yhi = max(y)+0.2*(max(y)-min(y))
  plt.plot(temperature, fwhm/cent, 'D', markersize=20)
  plt.axis([15, 35, ylo, yhi])
  plt.show()
  pdffile.savefig()
  pdb.set_trace() # type c to continue
  pdffile.close()

from numpy import *
import Image
from matplotlib.pyplot import *
import pdb
import scipy.stats as stat
from plotting import *
import transformations
import gaussfitter
import sys,time
import tifffile as tiff

#Returns an array of pixel counts
#First index = y coordinate
#Second index = x coordinate
def loadTIFFarray(filename):
    im = Image.open(filename)
    new = reshape(array(im.getdata()),(1024,1024))
    return new

#Subtract a background level from entire image
#User must supply a background region over which
#to average
def BGsub(im,xr,yr):
    bg = mean(im[yr[0]:yr[1],xr[0]:xr[1]])
    im = im-bg
    #Correct negative pixels
    for x in range(shape(im)[0]):
        for y in range(shape(im)[1]):
            if im[x,y] < 0:
                im[x,y] = 0
    return im

#Extract x,y coordinates for each photon from image data
#In order to limit array sizes, specify x and y ranges
#based on contour plot
#xr = two element list specifying range of second index
#yr = two element list specifying range of first index
def xyintercepts(im,xr,yr):
    #Loop through image pixels and adds to x, y arrays
    for xi in arange(xr[0],xr[1]+1):
        for yi in arange(yr[0],yr[1]+1):
            #Add current x,y coordinates to arrays N times,
            #where N is the value of the current pixel
            #If arrays do not already exist:
            if im[yi,xi]>0:
                try:
                    x = concatenate((x,repeat(xi,im[yi,xi])))
                    y = concatenate((y,repeat(yi,im[yi,xi])))
                except:
                    x = repeat(xi,im[yi,xi])
                    y = repeat(yi,im[yi,xi])
    return [x,y]

#Fit line to x, y intercept data
def xyline(x,y):
    #Compute unique x, y vectors
    yu = unique(y)
    xu = []
    w = []
    for i in yu:
        ind = where(y==i)
        w.append(size(ind))
        xu.append(mean(x[where(y==i)]))
    xu = array(xu)
    w = array(w)

    #Fit line
    fit = polyfit(xu,yu,1,w=w)

    return fit

#Gaussfit
def mygaussfit(x,y):
    h = myhist(y,bins=arange(min(x),max(x)+1,1))
    gy = h[0]
    gx = h[1]
    
    ampguess = gy.max()
    centroid = mean(y)
    widthguess = sqrt(mean((y-mean(y))**2))
    fit = gaussfitter.onedgaussfit(gx,gy,err=1+sqrt(gy+.75),params=[0,ampguess,\
                centroid,widthguess],fixed=[False,False,False,False])
    finex = linspace(gx[0],gx[-1],num=100)
    plot(finex,gaussfitter.onedgaussian(finex,*fit[0]))
    plot(gx,gy,'.')

#Calculate RMS scatter about line
def xyscatter(x,y,n=0,rotate=False):

    #If rotate keyword is True, rotate line to horizontal
    ang = 0.
    if rotate==True:
        #Fit line to data
        fit = xyline(x,y)
        #Rotate data to horizontal line
        resid = 100.
        while (abs(resid) > 50e-6):
            #Rotate to horizontal
            rot = transformations.rotation_matrix\
                  (-arctan(fit[0]),[0,0,1])[0:2,0:2]
            ang = ang + -arctan(fit[0]) #Keep track of rotation
            x,y = dot(rot,[x,y])
            #Fit new line
            fit = xyline(x,y)
            resid = fit[0]
    else:
        y = x

    #Create histogram
    #default n==0 implies to just bin by integers as in the CCD
    clf()
    if n==0:
        n = arange(min(y)-.5,max(y)+.5+1)
    h = myhist(y,bins=n)
    hist(y,bins=n)

    #Fit gaussian
    gx = h[1]
    gy = h[0]

    #Subtract BG - 5 bins on edges are BG
##    bg = mean(concatenate((gy[0:5],gy[-1:-6:-1])))
##    gy = gy  - bg
    
    ampguess = gy.max()
    centroid = mean(gx[where(gy==gy.max())])
    widthguess = 1.5#sqrt(mean((y-mean(y))**2))
    fit = gaussfitter.onedgaussfit(gx,gy,err=1+sqrt(gy+.75),params=[0,ampguess,\
                centroid,widthguess],fixed=[False,False,False,False])
    finex = linspace(gx[0],gx[-1],num=100)
    plot(finex,gaussfitter.onedgaussian(finex,*fit[0]))
    plot(gx,gy,'.')

    #return rms
    return fit[0][3],fit[2][3],ang,fit[0][2]

#Analyze Cr K lines with 0 Order
def CrK():
    resim = loadTIFFarray('goni1.791_0.5graze_5min_CrKalphabeta_0.72'
                          'yaw_LinesWithZero.TIF')
    #Analyze zero order
    x,y = xyintercepts(resim,[80,105],[430,533])
    zerowidth,err,ang = xyscatter(x,y,0)
    zerox = mean(x)

    #Analyze K beta
    x,y = xyintercepts(resim,[860,894],[360,500])
    betawidth,err,ang = xyscatter(x,y,0)
    betax = mean(x)

    #Analyze K alpha
    x,y = xyintercepts(resim,[930,975],[360,500])
    alphawidth,err,ang = xyscatter(x,y,0)
    alphax = mean(x)

    #Calculate dispersion (Angstrom/mm)
    betadiff = (betax-zerox)*13./1000.
    betadisp = (12400./5946.71)/betadiff

    alphadiff = (alphax-zerox)*13./1000. #millimeters of displacement
    alphadisp = (12400./5414.72)/alphadiff #angstrom/mm

    print 'Alpha Disp: ' + str(alphadisp)
    print 'Beta Disp: ' + str(betadisp)

    return alphadisp, betadisp

#Analyze Mg K 0 Order
def Mg0(n):
    im = loadTIFFarray('goni0.791_1.5graze_300sec_0order_0.72yaw_'
                       '25.60Drill_835.7Wall_825.3Vert.TIF')
    x,y = xyintercepts(im,[400,500],[300,600])

    width = xyscatter(x,y,n)
    return width

#Analyze Mg 1st order
def Mg1(n):
    im = tiff.open('goni0.791_1.5graze_300sec_1stOrder'
                       '_0.72yaw_29.80Drilll_Vert.TIF')
    #Analyze k alpha
    x,y = xyintercepts(im,[435,475],[580,700])
    alphawidth = xyscatter(x,y,0)[0]

    #Analyze k beta
    x,y = xyintercepts(im,[410,435],[580,700])
    betawidth = xyscatter(x,y,0)[0]

    print 'Alpha width: ' + str(alphawidth*2.35*13)
    print 'Beta width: ' + str(betawidth*2.35*13)

    return alphawidth, betawidth

disp = 0.20446944362612615 #From CrK analysis, AA/mm
fwhm = 2*sqrt(2*log(2)) #Convert from sigma to FWHM
#Analyze Mg 3rd Order
def Mg3():
    #Read in all 3rd order data, concatenate, and photon count
    im = tiff.imread('3rdOrderStack30sec60x.TIF')
    im2 = tiff.imread('3rdOrderStack60sec30x.TIF')
    im3 = tiff.imread('3rdOrderStack30sec20x.TIF')

    pim1 = photoncount(im)
    pim2 = photoncount(im2)
    pim3 = photoncount(im3)

    x = concatenate((pim1[0],pim3[0]))
    y = concatenate((pim1[1],pim3[1]))

    pim = array([x,y])
    
##    imc = concatenate((im,im3))
##    pim = array(photoncount(imc))

    pdb.set_trace()

    #Select only K alpha peak
    i = logical_and(pim[0]>50,pim[0]<100)
##    i = pim[0]<80
    pim = pim[:,i]

    pdb.set_trace()

    #Fit resolution before rotation
    res = xyscatter(pim[0],pim[1],0)
    unc = res[1]/res[0] #Fractional uncertainty in sigma
    width = fwhm * res[0] * 13.e-3 #mm of width
    position = (12400./1253.6)/disp * 3
    resolution = position/width
    print 'Position: ' + str(position)
    print 'Width:' + str(width)
    print 'Counts: ' + str(size(pim[0]))
    print 'Unrotated Width: ' + str(width)
    print 'Unrotated Bin Size: ' + str(res[0])
    print 'Unrotated Position: ' + str(position)
    print 'Unrotated Resolution: ' + str(resolution)
    print 'Fractional Uncertainty: ' + str(unc)
    print 'Mean Pixel: ' + str(res[3])

##    #Perform rotation and fit resolution
##    res = xyscatter(pim[0],pim[1],rotate=True)
##    width2 = fwhm * res[0] * 13.e-3 #mm of width
##    ang = res[2] * 180/pi
##    resolution2 = position/width2
##    print 'Rotated Angle: ' + str(ang)
##    print 'Rotated Width: ' + str(width2)
##    print 'Rotated Resolution: ' + str(resolution2)
##
##    print 'Rotation Boost: ' + str((resolution2-resolution)/resolution)

#Investigate 60sec30x
def inv60sec30x():
    im = tiff.imread('3rdOrderStack30sec20x.TIF')
    bg = mean(im) + 5 * std(im)
    pdb.set_trace()

    #Loop through images, photon count, and determine mean of counts
    #within 50-100 x pixels
    m = []
    for i in range(shape(im)[0]):
##        if stat.poisson.cdf(sum(im[i]>bg),mu) > .999:
##            print 'Threw out: ' + str(i)
##            continue #Disregard image if it has too many counts
        x = []
        y = []
        charge = []
        while (im[i].max() > bg):
            #Sum charge in 3x3 array
            ind1,ind2 = where(im[i]==im[i].max())
            ind1 = ind1[0]
            ind2 = ind2[0] #Handles multiple pixels with same ADU
            charge.append(sum(im[i][ind1-1:ind1+2,ind2-1:ind2+2]))
            x.append(ind2)
            y.append(ind1)
            im[i][ind1-1:ind1+2,ind2-1:ind2+2] = 0.
        x = array(x)
        i = logical_and(x<100,x>50)
        if sum(i)>0:
            m.append(mean(x[logical_and(x<100,x>50)]))
    pdb.set_trace()

#Analyze Mg 4th Order
def Mg4():
    #Read in 4th order data and photon count
    im = tiff.imread('4thOrderStack30sec60x.TIF')
    pim = array(photoncount(im,toss=True))

    #Select only K alpha peak
    i = logical_and(pim[0]>77,pim[0]<155)
    pim = pim[:,i]

    #Fit resolution before rotation
    width = fwhm * xyscatter(pim[0],pim[1],0)[0] * 13.e-3 #mm of width
    position = (12400./1253.6)/disp * 4
    resolution = position/width
    print 'Unrotated Width: ' + str(width)
    print 'Unrotated Position: ' + str(position)
    print 'Unrotated Resolution: ' + str(resolution)

    #Perform rotation and fit resolution
    res = xyscatter(pim[0],pim[1],rotate=True)
    width2 = fwhm * res[0] * 13.e-3 #mm of width
    ang = res[2] * 180/pi
    resolution2 = position/width2
    print 'Rotated Angle: ' + str(ang)
    print 'Rotated Width: ' + str(width2)
    print 'Rotated Resolution: ' + str(resolution2)

    print 'Rotation Boost: ' + str((resolution2-resolution)/resolution)

#Photon count a TIFF image, assume low counts (~ 10)
def photoncount(img,toss=False):
    im = copy(img)
    bg = mean(im) + 5 * std(im)
    x = []
    y = []
    charge = []

    #Test for bad images
    n = []
    for i in range(shape(im)[0]):
        n.append(sum(im[i]>bg))
    mu = mean(n) #Parameter for Poisson distribution of count rate
    
    for i in range(shape(im)[0]):
        if toss==True:
            if stat.poisson.cdf(sum(im[i]>bg),mu) > .999:
                print 'Threw out: ' + str(i)
                continue #Disregard image if it has too many counts
        bg = mean(im[i]) + 5 * std(im[i])
        while (im[i].max() > bg):
            #Sum charge in 3x3 array
            ind1,ind2 = where(im[i]==im[i].max())
            ind1 = ind1[0]
            ind2 = ind2[0] #Handles multiple pixels with same ADU
            charge.append(sum(im[i][ind1-1:ind1+2,ind2-1:ind2+2]))
            x.append(ind2)
            y.append(ind1)
            im[i][ind1-1:ind1+2,ind2-1:ind2+2] = 0.

    return x,y,charge
        
#Bin size effect simulation
def bintest(n,binsize,offset,plt=False,dofit=True,bg=0.):
    #Generate random numbers
    r = stat.norm.rvs(size=n)
    #Generate uniform background - how many counts and over what interval?
    
    #Select bin edges, need at least 5 bins
    binpos = arange(0,ceil(max(r)),binsize)
    binneg = arange(0,floor(min(r)),-binsize)
    bins = concatenate((binneg[-1:0:-1],binpos))
    #Add on extra bins to sides to make sure 5 bins
    deficit = 5 - (size(bins)-1)
    if deficit > 0:
        for i in range(int(deficit)):
            bins = append(bins,bins.max()+binsize)
    
    #Add offset to put centroid at offset in a bin
    bins = bins - offset*binsize

    #Create histogram
    h = myhist(r,bins=bins)
    centroid = mean(r)

    #Fit gaussian
    if dofit==True:
        gx = h[1]
        gy = h[0]
        widthguess = sqrt(mean((r-mean(r))**2))
        ampguess = n/sqrt(2*pi)/widthguess
        err = 1+sqrt(gy+.75)
        fit = gaussfitter.onedgaussfit(gx,gy,err=err\
                            ,params=[0,ampguess,\
                    centroid,widthguess],fixed=[True,False,False,False])
        finex = linspace(gx[0],gx[-1],num=100)
        if plt==True:
            plot(finex,gaussfitter.onedgaussian(finex,*fit[0]))
            plot(gx,gy,'.')
    ##    
    ##    #Figure out fitted mean and uncertainty
    ##    print 'Mean: ' + str(fit[0][2])
    ##    print 'Std: ' + str(fit[0][3])

        #Fit results
        fitmean = fit[0][2]
        fitsigma = fit[0][3]
        
    momentsigma = std(r)

    if dofit==True:
        return (fitmean, fitsigma, centroid, momentsigma)
    else:
        return (centroid, momentsigma)

#Write function to test a shift in bin edges on the SAME data
#Generate data, set up initial bin edges using offset, generate second
#set of bin edges with shift
#Fit both bins and calculate DIFFERENCE in standard deviation and mean
def binshifttest(n,binsize,offset,shift,plt=False):
    #Generate random numbers
    r = stat.norm.rvs(size=n)
    #Select bin edges
    binpos = arange(0,ceil(max(r)),binsize)
    binneg = arange(0,floor(min(r)),-binsize)
    bins = concatenate((binneg[-1:0:-1],binpos))
    #Add offset to put centroid at offset in a bin
    bins = bins - offset*binsize

    #Calculate shifted bins
    bins2 = bins - shift*binsize

    #Create histogram
    clf()
    h = myhist(r,bins=bins)
    hist(r,bins=bins,histtype='step')
    h2 = myhist(r,bins=bins2)
    hist(r,bins=bins2,histtype='step')

    #Fit gaussian
    gx = h[1]
    gy = h[0]
    centroid = mean(r)
    widthguess = sqrt(mean((r-mean(r))**2))
    ampguess = n/sqrt(2*pi)/widthguess
    err = 1+sqrt(gy+.75)
    pdb.set_trace()
    fit = gaussfitter.onedgaussfit(gx,gy,err=err\
                        ,params=[0,ampguess,\
                centroid,widthguess],fixed=[True,False,False,False])
    pdb.set_trace()
    gx = h2[1]
    gy = h2[0]
    centroid = mean(r)
    widthguess = sqrt(mean((r-mean(r))**2))
    ampguess = n/sqrt(2*pi)/widthguess
    err = 1+sqrt(gy+.75)
    fit2 = gaussfitter.onedgaussfit(gx,gy,err=err\
                        ,params=[0,ampguess,\
                centroid,widthguess],fixed=[True,False,False,False])
    
    finex = linspace(gx[0],gx[-1],num=100)
    if plt==True:
        plot(finex,gaussfitter.onedgaussian(finex,*fit[0]))
        plot(finex,gaussfitter.onedgaussian(finex,*fit2[0]))
        plot(gx,gy,'.')

    #Compute difference in standard deviation and mean
    meandiff = fit[0][2] - fit2[0][2]
    stddiff = fit[0][3] - fit2[0][3]

    return meandiff, stddiff

#Perform bin shift test simulation, characterizing distributions
#of mean and std differences
def binshiftsim(m,n,binsize,offset,shift):
    mdiff = []
    sdiff = []
    for i in range(m):
        res = binshifttest(n,binsize,offset,shift)
        mdiff.append(res[0])
        sdiff.append(res[1])

    #Convert to numpy arrays
    return array(mdiff), array(sdiff)

def binsim(m,n,binsize,offset):
    avg = []
    std = []
    for i in range(m):
        res = bintest(n,binsize,offset)
        avg.append(res[0])
        std.append(res[2])
        print i
        sys.stdout.flush()
    return avg,std

def randomoffset(m,n,binsize,report=False,dofit=True):
    fitavg = []
    fitstd = []
    momavg = []
    momstd = []
    start = time.time()
    for i in range(m):
        res = bintest(n,binsize,random.uniform(),dofit=dofit)
        fitavg.append(res[0])
        fitstd.append(res[1])
        momavg.append(res[2])
        momstd.append(res[3])
        if report==True:
            sys.stdout.write(str(float(i)/m)+'\r')
            sys.stdout.flush()

    print 'Time elapsed: ' + str(time.time()-start)
    return array(fitavg),array(fitstd),array(momavg),array(momstd)

def scanparam(m,n,bins,offset):
    #Initialize avg and std arrays
    avg = zeros((size(bins),m))
    std = copy(avg)
    for b in range(size(bins)):
        avg[b],std[b]=binsim(m,n,bins[b],offset)
        hist(std[b],bins=10,label=str(bins[b]),histtype='step')
        print b
        sys.stdout.flush()
    legend(loc='upper left')

    return (avg,std)

#Run simulation similar to Valentine et al.
#Reproduce figures 1 and 2
#Also compute bias
def valentinesim():
    #Loop through number of counts from 10 to 100 in steps of 10
    #Also loop through relative bin size from 1 to .001 in steps of .005
    #At each step, compute centroid bias and standard deviation and
    #sigma bias and standard deviation
    #Use 200 simulations at each step to achieve sigma std of 
    #Finally, plot results as contour plots to get complete dependence

    #Scan arrays
    binscan = 1./arange(1,26)
    Nscan = arange(10,110,10)

    #Initialize result arrays
    fitmbias = zeros((size(binscan),size(Nscan)))
    mommbias = copy(fitmbias)
    fitmsigma = copy(fitmbias)
    mommsigma = copy(fitmbias)

    fitsbias = copy(fitmbias)
    momsbias = copy(fitmbias)
    fitssigma = copy(fitmbias)
    momssigma = copy(fitmbias)
    
    for i in range(size(binscan)):
        for j in range(size(Nscan)):
            fitavg,fitstd,momavg,momstd = randomoffset(200,Nscan[j],binscan[i],\
                                                       report=False)
            #Avg results
            fitmbias[i,j] = 1.-mean(fitavg)
            mommbias[i,j] = 1.-mean(momavg)
            fitmsigma[i,j] = std(fitavg)
            mommsigma[i,j] = std(momavg)
            #Std results
            fitsbias[i,j] = 1.-mean(fitstd)
            momsbias[i,j] = 1.-mean(momstd)
            fitssigma[i,j] = std(fitstd)
            momssigma[i,j] = std(momstd)
            sys.stdout.write('Bin: ' + str(binscan[i]) + '\tCounts: ' +\
                             str(Nscan[j]) + '\r')
            sys.stdout.flush()

    return fitmbias,mommbias,fitmsigma,mommsigma,fitsbias,momsbias,\
           fitssigma,momssigma

#1st order y deviation
def ypos(yaw,energy):
    alpha = arctan(8127.*(yaw*pi/180)/70.9)
    beta = arcsin(1240./energy/160./sin(.5*pi/180) - sin(alpha))
    return 8127. * sin(.5*pi/180) * cos(beta) / cos(alpha)

def mindiff(x,y):
    diff = abs(x-y)
    diff2 = diff[invert(isnan(diff))]
    return where(diff==min(diff2))[0]

def findMoments(d):
    x,y = meshgrid(arange(shape(d)[1]),arange(shape(d)[0]))
    cx = nansum(x*d)/nansum(d)
    cy = nansum(y*d)/nansum(d)
    rmsx = nansum((x-cx)**2*d)/nansum(d)
    rmsy = nansum((y-cy)**2*d)/nansum(d)
    
    return cx,cy,sqrt(rmsx),sqrt(rmsy)

#Automatically find a spot centroid after clicking two coordinates on the image
class clickCentroid:
    def __init__(self,img):
        self.img = img
        self.x0 = None
        self.y0 = None
        self.x1 = None
        self.y1 = None
        self.gx = []
        self.gy = []
        self.fig = gcf()
        self.con = self.fig.canvas.mpl_connect('button_press_event',self.clickEvent)
    def clickEvent(self,event):
        #If x0 and y0 are undefined, set them and return
        if self.x0 is None:
            #Define first point
            self.x0 = event.xdata
            self.y0 = event.ydata
            return
        #If x1 and y1 are undefined, set them and return centroid
        #Define second point
        self.x1 = event.xdata
        self.y1 = event.ydata
        #Order points properly
        x0 = min([self.x0,self.x1])
        x1 = max([self.x0,self.x1])
        y0 = min([self.y0,self.y1])
        y1 = max([self.y0,self.y1])            
        #Compute centroid between coordinates and update centroid list
        cx,cy,rmsx,rmsy = findMoments(self.img[y0:y1,x0:x1])
        print 'X: ' + str(cx+x0)
        print 'Y: ' + str(cy+y0)
        print 'RMS X: ' + str(rmsx)
        print 'RMS Y: ' + str(rmsy)
        try:
            self.gx.append(cx+x0)
            self.gy.append(cy+y0)
        except:
            self.gx = [cx+x0]
            self.gy = [cy+y0]
        #Fit a Gaussian in X and Y
        x,y = xyintercepts(self.img,[x0,x1],[y0,y1])
        dist,bins = myhist(x,bins=50.)
        fit = gaussfitter.onedgaussfit(bins,dist,\
                            params=[0.,max(dist),x0+cx,rmsx])
        print 'Gauss Params: '
        print fit[0]
        figure()
        plot(bins,dist)
        plot(bins,fit[1])
        #Reset event
        self.x0 = None
        self.x1 = None
        self.y0 = None
        self.y1 = None
            
    def close(self):
        self.fig.canvas.mpl_disconnect(self.con)

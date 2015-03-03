from numpy import *
import transformations as tr
import pdb,glob
from plotting import *

#Read in wavefront
def wavefront(f):
    #Open file
    f = open(f)
    d = f.readlines()
    d = d[-22:]

    #Set up x and y vectors
    y = d[0]
    y = y.split(',')
    y = array(y[1:-1])
    y = y.astype('float')
    x = []
    for i in arange(1,size(d)):
        s = d[i]
        s = s.split(',')
        x.append(s[0])
    x = array(x)
    x = x.astype('float')

    #Set up wavefront array
    wave = d[-21:]
    arr = zeros((size(x),size(y)))
    for i in range(size(x)):
        arr[i,:] = wave[i].split(',')[1:-1]    

    return (x,y,arr)

def randyplots():
    #Make flat vs. figure plots for Randy
    x1,y1,arr1 = wavefront('/Users/rallured/ORGS/WFS/flat.csv')
    x2,y2,arr2 = wavefront('/Users/rallured/ORGS/WFS/130313-WavefrontStep1.csv')
    clf()
    mycontour(arr1+1,x=x1,y=y1,fmt='%.2f um')
    title('Direct View of Beam')
    xlabel('X (mm)')
    ylabel('Y (mm)')
    savefig('/Users/rallured/ORGS/WFS/flat.eps')
    clf()
    mycontour(arr2+1,x=x2,y=y2,fmt='%.2f um')
    title('Figure of Si Wafer')
    xlabel('X (mm)')
    ylabel('Y (mm)')
    savefig('/Users/rallured/ORGS/WFS/figure.eps')

    #Angular step plot
    x1, y1 = readZernike('130314AngStep1z.txt')
    x2, y2 = readZernike('130314AngStep2z.txt')
    x3, y3 = readZernike('130314AngStep3z.txt')
    x4, y4 = readZernike('130314AngStep4z.txt')

    avgs = array([mean(x1),mean(x2),mean(x3),mean(x4)])
    unc = array([std(x1),std(x2),std(x3),std(x4)])
    unc = unc/10.

    clf()
    errorbar(arange(1,5),avgs,yerr=unc,fmt='.')
    title('Angular Step Demonstration')
    xlabel('Step #')
    ylabel('Arbitrary Angle (arcseconds)')
    xlim([0,5])
    ylim([6,17])
    savefig('AngSteps.eps')

#Read in WFS LabVIEW data
def readLabVIEW(f):
    #Read in file
    data = genfromtxt('WavefrontTest.txt')
    xsize = int(genfromtxt('WavefrontX.txt'))
    ysize = int(genfromtxt('WavefrontX.txt'))
    pdb.set_trace()
    #Strip off 0's
    data = data[0:xsize,0:ysize]
    return data

#Calculate surface normal
def normal(data):
    #Housekeeping
    pitch = 150 #150 micron pitch
    xslope = []
    yslope = []
    
    #Determine x and y slopes
    for i in range(size(x)):
        for j in range(size(y)-1):
            yslope.append((arr[i,j]-arr[i,j+1])/pitch)
            xslope.append((tarr[i,j]-tarr[i,j+1])/pitch)
    xslope = mean(xslope)
    yslope = mean(yslope)
    xang = arctan(xslope)
    yang = arctan(yslope)

    #Calculate surface normal
    xhat = [cos(xang),0,sin(xang)]
    yhat = [0,cos(yang),sin(yang)]
    norm = cross(xhat,yhat)
    pdb.set_trace()
    norm = norm/linalg.norm(norm)
    return norm

#Read Zernike data
def readZernike(file):
    z = genfromtxt(file)
    z = transpose(z)
    zrad = (z[2]+z[3])*250.
    xang = arctan(z[1]*2/zrad)*180/pi*60**2
    yang = arctan(z[0]*2/zrad)*180/pi*60**2
    return (xang,yang)

def hysteresis130808():
    steps = array([0,10,20,30,40,50,60])
    mx = []
    my = []
    stdx = []
    stdy = []
    for file in glob.glob('/Users/rallured/Dropbox/Hysteresis/130809Test*'):
        x,y=readZernike(file)
        mx.append(mean(x))
        my.append(mean(y))
        stdx.append(std(x)/10.)
        stdy.append(std(y)/10.)

    clf()
    errorbar(steps,mx,yerr=stdx,label='Pitch Angle',fmt='.')
    errorbar(steps,my,yerr=stdy,label='Roll Angle',fmt='.')
    xlim([-5,65])
    title('8/8/13 Step Size and Hysteresis Test')
    xlabel('Encoder Steps')
    ylabel('Arcsec')
    legend(loc='lower left')
    savefig('/Users/rallured/Dropbox/Hysteresis/130809HysteresisResults.png')

def hysteresis130813():
    steps = array([0,10,20,30,40,50,60])
    mx = []
    my = []
    stdx = []
    stdy = []
    for file in glob.glob('/Users/rallured/Dropbox/Hysteresis/RollHysteresis_081313/130813*.dat'):
        x,y=readZernike(file)
        mx.append(mean(x))
        my.append(mean(y))
        stdx.append(std(x)/10.)
        stdy.append(std(y)/10.)

    clf()
    errorbar(steps,mx-mean(mx),yerr=stdx,label='Pitch Angle',fmt='.')
    errorbar(steps,my-mean(my),yerr=stdy,label='Roll Angle',fmt='.')
    xlim([-5,65])
    title('8/13/13 Step Size and Hysteresis Test (Pitch)')
    xlabel('Encoder Steps')
    ylabel('Arcsec')
    legend(loc='upper center')
    savefig('/Users/rallured/Dropbox/Hysteresis/130813HysteresisResults1.png')

    mx = []
    my = []
    stdx = []
    stdy = []
    for file in glob.glob('/Users/rallured/Dropbox/Hysteresis/RollHysteresis_081313_2/130813*.dat'):
        x,y=readZernike(file)
        mx.append(mean(x))
        my.append(mean(y))
        stdx.append(std(x)/10.)
        stdy.append(std(y)/10.)

    clf()
    errorbar(steps,mx-mean(mx),yerr=stdx,label='Pitch Angle',fmt='.')
    errorbar(steps,my-mean(my),yerr=stdy,label='Roll Angle',fmt='.')
    xlim([-5,65])
    title('8/13/13 Step Size and Hysteresis Test (Roll)')
    xlabel('Encoder Steps')
    ylabel('Arcsec')
    legend(loc='upper center')
    savefig('/Users/rallured/Dropbox/Hysteresis/130813HysteresisResults2.png')

###Calculate surface normal
##def normal(x,y,arr):
##    #Housekeeping arrays
##    tarr = transpose(arr)
##    xstep = x[1]-x[0]
##    ystep = y[1]-y[0]
##    xslope = []
##    yslope = []
##
##    #Account for units
##    arr = arr * 1e-6
##    tarr = tarr * 1e-6
##    xstep = xstep * 1e-3
##    ystep = ystep * 1e-3
##    
##    #Determine x and y slopes
##    for i in range(size(x)):
##        for j in range(size(y)-1):
##            yslope.append((arr[i,j]-arr[i,j+1])/ystep)
##            xslope.append((tarr[i,j]-tarr[i,j+1])/xstep)'
##    xslope = mean(xslope)
##    yslope = mean(yslope)
##    xang = arctan(xslope)
##    yang = arctan(yslope)
##
##    #Calculate surface normal
##    xhat = [cos(xang),0,sin(xang)]
##    yhat = [0,cos(yang),sin(yang)]
##    norm = cross(xhat,yhat)
##    pdb.set_trace()
##    norm = norm/linalg.norm(norm)
##    return norm

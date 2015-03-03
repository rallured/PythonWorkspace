from numpy import *
import os
import pdb
from matplotlib.pyplot import *
from gaussfitter import *
from scipy.optimize import *
from plotting import *
import matplotlib

def stepsize():
    os.chdir('/Users/rallured/PythonWorkspace/Calibration/ReflectedBeam')

    #Examine hysteresis and step size for (.6,.05) pulse
    #Read in forward scan
    f = genfromtxt('121002/121002FineScan1.txt')
    b = genfromtxt('121002/121002FineScan2.txt')
    b = b[:,1:] #Subtract off null first value
    fwd = f[1]/f[0]
    background = mean(fwd[0:50])
    fwd = fwd-background
    bck = b[1]/b[0]
    background = mean(bck[0:50])
    bck = bck-background
    x = arange(size(bck))

    #Step size
    w = fwd*x
    kpeak = sum(w[72:160])/sum(fwd[72:160])
    lpeak = sum(w[188:281])/sum(fwd[188:281])
    forx = (arange(300)-kpeak)*(0.0464/2)+45.
    fstepsize = 0.861/(lpeak-kpeak)
##    fstepsize = 6.373/(lpeak-kpeak)
    print 'kpeak: ' + str(kpeak)
    print 'lpeak: ' + str(lpeak)
    print 'forward: ' + str(fstepsize)
    w = fwd*(x-kpeak)**2
    beamwidth1 = sum(w[72:160])/sum(fwd[72:160])
    w = bck*(x-lpeak)**2
    beamwidth2 = sum(w[188:281])/sum(bck[188:281])
    print 'beamwidth1: ' + str(beamwidth1)
    print 'beamwidth2: ' + str(beamwidth2)
##    forx = arange(44.105+lpeak*fstepsize,\
#                   44.105-(size(w)-lpeak-1)*fstepsize,-fstepsize)
##    fstepsize = .0464
##    forx = arange(38.163+lpeak*fstepsize,\
##                   38.163-(size(w)-lpeak-.5)*fstepsize,-fstepsize)
    pdb.set_trace()

    w = bck*x
    kpeak = mean(w[209:294])
    lpeak = mean(w[82:190])
    kpeak = sum(w[209:294])/sum(bck[209:294])
    lpeak = sum(w[82:190])/sum(bck[82:190])
    bstepsize = .861/(kpeak-lpeak)
    print 'kpeak: ' + str(kpeak)
    print 'lpeak: ' + str(lpeak)
    print 'back: ' + str(bstepsize)
    w = fwd*(x-kpeak)**2
    beamwidth1 = sum(w[209:294])/sum(fwd[209:294])
    w = bck*(x-lpeak)**2
    beamwidth2 = sum(w[82:190])/sum(bck[82:190])
    print 'beamwidth1: ' + str(beamwidth1)
    print 'beamwidth2: ' + str(beamwidth2)
    backx = arange(44.105-(size(w)-lpeak)*bstepsize,\
                   44.105+lpeak*bstepsize,bstepsize)

    #Make plots
    p = genfromtxt('PolData.txt')
    t = p[:size(p)/2]
    p = p[size(p)/2:]
    ion()
    clf()
    plot(forx,fwd)
    title('Polarizer Scan')
    xlabel('Incidence angle (Deg)')
    ylabel('Normalized Counts')
##    pdb.set_trace()
##    #savefig('/Users/rallured/PythonWorkspace/Calibration/ReflectedBeam/ForFineScan.eps')
##    clf()
##    plot(backx,bck)
##    title('12/10/02 Backward Scan')
##    xlabel('Incidence angle (Deg)')
##    ylabel('Normalized Counts')
##    #plot([71.255*bstepsize,71.255*bstepsize],[-.2,1.6],'k--')
##    plot([backx[0]+71.255*bstepsize,backx[0]+71.255*bstepsize],[-.2,1.6],'k--')
##    text(.6,1,'Hysteresis')
##    #savefig('/Users/rallured/PythonWorkspace/Calibration/ReflectedBeam/BckFineScan.eps')
    #Pol Scan and Polarization Double Plot
    fig = figure()
    ax1 = fig.add_subplot(111)
    ax1.plot(forx,fwd,'b-')
    ax1.set_xlabel('Incidence Angle (deg)')
    ax1.set_ylabel('Normalized Counts',color='b')
    for t1 in ax1.get_yticklabels():
        t1.set_color('b')
    ax2 = ax1.twinx()
    ax2.plot(t,p,'r')
    ax2.set_ylabel('Polarization (fractional)',color='r')
    for t1 in ax2.get_yticklabels():
        t1.set_color('r')
    ax1.set_title('Polarizer Response')

def kv10spec():
    os.chdir('/Users/rallured/PythonWorkspace/Calibration/VTarget/')
    hv10 = genfromtxt('111020_1680V_10G_10kV.mca',skip_header=14,skip_footer=1)
    x = arange(256) + 1

##    clf()
##    x = arange(256)+1
##    plot(x,hv10)
    
    fit = onedgaussfit(x[175:-8],hv10[175:-8],params=[0,35,205,50])
    x = x*4.950/fit[0][2] #Calibrate to energy scale

    clf()
    plot(x,hv10)
##    plot(x[175:-8],fit[1],'--')
    xlabel('Energy (keV)')
    ylabel('Counts')
    title('Spectrum of Vanadium Target')
    savefig('10kvSpec.eps')

def gaincurve():
    os.chdir('/Users/rallured/PythonWorkspace/Calibration/')
    cent = array([108.76,177.95,305.42,166.23,285.65,446])
    cent[0:3] = cent[0:3]*96.05/305.42
    volt = array([1602.,1650,1700,1750,1800,1850])
    #Convert centroids to charge
    ccr = 0.270*2E-12/41.27 #coulombs/channel ratio
    cent = cent*ccr
    #Convert centroids to gain
    ne = 5900./26.2 #Initial number of electrons
    q0 = ne * 1.60E-19 #Initial amount of charge
    #Converts charges to gain
    gain = cent/q0
    #Make gain curve plot
    clf()
    semilogy(volt/1000.,gain,'.')
    #Fit line to log plot
    fit = polyfit(log10(volt[:-1]),log10(gain[:-1]),1)
    fgain = fit[1]+fit[0]*log10(arange(1500,1900))
    plot(arange(1500,1900)/1000.,10**fgain,'--')
    xlim([1.550,1.900])
    ylim([10**4,10**6])
    title('SWPC Gain Curve')
    xlabel('Anode Voltage (kV)')
    ylabel('Detector Gain')
##    savefig('SWPCGainCurve.eps')
    return (volt,gain)

def normalization():
    os.chdir('/Users/rallured/PythonWorkspace/Calibration/ReflectedBeam/121002/')
    d = genfromtxt('121002StabilityData.txt')
    bg1 = 17.2 # pm 0.927  cts/30 sec
    bg2 = 24.9 # pm 0.70569  cts/30 sec
    m1 = d[0] - bg1
    m2 = d[1] - bg2

    clf()
    plot(m1,label='Normalization Detector')
    plot(m2,label='Sealed Detector')
    xlabel('Trial #')
    ylabel('Counts')
    legend(loc='lower left')
    title('Simultaneous Detector Data')
    savefig('Simultaneous.eps')
    clf()
    hist(m2/m1,bins=10)
    title('Counts Ratio Distribution')
    text(.825,50,'Mean: '+str(round(mean(m2/m1),2))+\
         '\nStd: '+str(round(std(m2/m1),3)),fontsize=16)
    savefig('Ratio.eps')
    
    sigc1 = std(d[0])
    sigc2 = std(d[1])
    sigb1 = sqrt(bg1)
    sigb2 = sqrt(bg2)
    sigm1 = sqrt(mean(m1))
    sigm2 = sqrt(mean(m2))
    sigs2 = sqrt((sigc2)**2-(sigm2)**2-(sigb2)**2)
    print sigs2/mean(m2)
    sigs1 = sqrt((sigc1)**2-(sigm1)**2-(sigb1)**2)
    print sigs1/mean(m1)

    return (m1,m2)

#Fit cos squared to a modulation curve
def cos2(x,const,amp,phase):
    return const + amp * cos(x-phase)**2

def fitcos2(x,y,err):
    constguess = 0.
    ampguess = np.max(y)
    phaseguess = x[where(y == np.max(y))][0]
    pdb.set_trace()
    fit = curve_fit(cos2,x,y,p0=[constguess,ampguess,\
                                      phaseguess],sigma=err)
    fdata = cos2(x,fit[0][0],fit[0][1],fit[0][2])
    chi = sum(((fdata-y)/err)**2)/(size(x)-3)
    return (fit,fdata,chi)

def modcurve():
    os.chdir('/Users/rallured/PythonWorkspace/Calibration/121031/')
    c1 = []
    c2 = []
    deg = arange(0,360,10)
    for i in deg:
        d1 = genfromtxt('121031MCA1_'+str(i)+'deg.mca',\
                        skip_header=14,skip_footer=1)
        d2 = genfromtxt('121031MCA2_'+str(i)+'deg.mca',\
                        skip_header=14,skip_footer=1)
        c1.append(sum(d1))
        c2.append(sum(d2[60:]))
    c1 = array(c1)
    c2 = array(c2)
    b2 = genfromtxt('121031MCA2_Dark.mca',skip_header=14,skip_footer=1)

    #Make bg plot
    clf()
    plot(arange(1,257),b2,label='BG Spectrum')
    title('10/31/12 Background Spectrum')
    xlabel('Channel')
    ylabel('Counts')
    plot([60,60],[0,4],'k--')
    d = genfromtxt('121031MCA2_0deg.mca',skip_header=14,skip_footer=1)
    plot(arange(1,257),d/4.,label='$0^\circ$ Spectrum')
    legend(loc='upper right')
    savefig('BGSpec.eps')
    

    b2 = sum(b2[60:])
##    b2 = 30.
    print b2
    b1 = sum(genfromtxt(\
        '/Users/rallured/PythonWorkspace/Calibration/121106MCA1Dark_160thresh.mca'\
        ,skip_header=14,skip_footer=1))
##    b2 = b2 - sqrt(b2)
    bg1err = sqrt(b1)/600.*300 * 1.55
    bg2err = sqrt(b2)*1.55
    b1 = b1/600.*300
    m1 = c1-b1
    m2 = c2-b2
    m2[0] = m2[0]-b2
    print 'ok1'
    c1err = sqrt(c1+(m1*.0198)**2+bg1err**2)
    print 'ok2'
    c2err = sqrt(c2+(m2*.0277)**2+bg2err**2)
    print 'ok3'
    r = m2/m1
    
    rerr = sqrt((c1err*m2/m1**2)**2+(c2err/m1)**2)

    clf()
    errorbar(deg,r,yerr=rerr,fmt='.')

    xlabel('Phase Angle (deg)')
    ylabel('Normalized Counts')
    title('BRP Prototype Modulation Curve')
    xlim([-10,360])

    #Average over modulation period
    c1avg = arange(18)
    c2avg = arange(18)
    c1aerr = arange(18)
    c2aerr = arange(18)
    for i in range(18):
        c1avg[i] = c1[i]+c1[i+18]
        c2avg[i] = c2[i]+c2[i+18]
        c1aerr[i] = sqrt(c1err[i]**2+c1err[i+18]**2)
        c2aerr[i] = sqrt(c2err[i]**2+c2err[i+18]**2)
    b1 = b1*2
    b2 = b2*2
    m1 = c1avg-b1
    m2 = c2avg-b2
    m2[0] = m2[0]-b2/2.
##    c1err = sqrt(c1avg+(m1*.0198)**2+(bg1err*2)**2)
##    c2err = sqrt(c2avg+(m2*.0277)**2+(bg2err*2)**2)
    r = m2/m1
    
    rerr = sqrt((c1aerr*m2/m1**2)**2+(c2aerr/m1)**2)
    f = fitcos2(arange(0,180,10)*pi/180,r,rerr)

    #Adjust for negative points
    for i in range(size(r)):
        if r[i] < 0:
            r[i] = 0

    clf()
    errorbar(arange(0,180,10),r,yerr=rerr,fmt='.')
    plot(arange(0,180,10),f[1])
    title('10/31/12 Folded Modulation Curve')
    xlabel('Polarization Angle (deg)')
    ylabel('Normalized Counts')
    xlim([-10,180])

    return (m1,m2,rerr,c1,c2)

#Plot measured and fitted reflectance curves for M1-120605
def m1_120605():
    os.chdir('/Users/rallured/PythonWorkspace/Calibration/m1-120605/')
    meas = transpose(genfromtxt('pat039587_corr2.txt'))
    fit = transpose(genfromtxt('CenterFit.txt',comments=';'))
    up = transpose(genfromtxt('pat039590_corr.abs'))
    down = transpose(genfromtxt('pat039589_corr.abs'))
    left = transpose(genfromtxt('pat039591_corr.abs'))
    right = transpose(genfromtxt('pat039592_corr.abs'))
    de = meas[0][1]-meas[0][0]
    mi = sum(meas[1]*de)
    de = up[0][1]-up[0][0]
    upi = sum(up[1]*de)
    de = down[0][1]-down[0][0]
    downi = sum(down[1]*de)
    de = left[0][1]-left[0][0]
    lefti = sum(left[1]*de)
    de = right[0][1]-right[0][0]
    righti = sum(right[1]*de)
    print mean([mi,upi,downi,lefti,righti])

    #Load in final fit data
    finalfit = np.transpose(np.genfromtxt('m1-120605FinalFit.txt',comments=';'))
    ind = np.where(np.logical_and(finalfit[0]>=490,finalfit[0]<=550))
    fit = [finalfit[0][ind],finalfit[1][ind]]
    ind = np.where(np.logical_and(finalfit[0]>=40,finalfit[0]<=1400))
    #fit = [finalfit[0][ind],finalfit[2][ind]]

    clf()
    plot(meas[0],meas[1],'.',markersize=14,label='ALS Data')
##    plot(up[0],up[1],'--',markersize=14,label='Y=+10mm')
##    plot(down[0],down[1],'--',markersize=14,label='Y=-10mm')
##    plot(left[0],left[1],'--',markersize=14,label='X=-10mm')
##    plot(right[0],right[1],'--',markersize=14,label='X=+10mm')
    plot(fit[0],fit[1],'-',label='IMD Fit')
    xlim([meas[0][0],meas[0][-1]])
##    legend(loc='upper right')
    title('M1-120605 Multilayer Response')
    xlabel('Energy (eV)')
    ylabel('Reflectance (fractional)')
    text(530,.01,'d=17.16 $\AA$',size=16)
    text(530,.009,'$\sigma$=2.59 $\AA$ rms',size=16)
    text(530,.008,'$\delta E$=1.80 eV',size=16)
##    savefig('/Users/rallured/Documents/GEMS/ReflectorPaper/M1-120605.eps')
    pdb.set_trace()

    
    clf()
    t1 = transpose(genfromtxt('pat039600.abs'))
    t2 = transpose(genfromtxt('pat039601.abs'))
    tmodel = transpose(genfromtxt('/Users/rallured/IDLWorkspace82/'+\
            'Multilayer/LLNL/m1-110525/LowETrans.txt',skip_header=24))
    plot(t1[0],t1[1])
    plot(t2[0],t2[1])
    plot(tmodel[0],tmodel[1])
    pdb.set_trace()

##    clf()
##    ang = genfromtxt('m1-120605_511eV.txt')
##    t = ang[:601]
##    rs = ang[601:601+601]
##    rp = ang[-601:]
##
##    fig = figure()
##    ax1 = fig.add_subplot(111)
##    ax1.semilogy(t,rs,'-b',label='S-Polarization')
##    ax1.semilogy(t,rp,'--b',label='P-Polarization')
##    ax1.set_xlabel('Incidence Angle (deg from normal)')
##    ax1.set_ylabel('Reflectance (fractional)',color='b')
##    ax1.legend(loc='center')
##    for t1 in ax1.get_yticklabels():
##        t1.set_color('b')
##    ax2 = ax1.twinx()
##    ax2.plot(t,(rs-rp)/(rs+rp),'r')
##    ax2.set_ylabel('Modulation Factor',color='r')
##    for t1 in ax2.get_yticklabels():
##        t1.set_color('r')
##    ax1.set_title('M1-120605 511.3 eV Response')

def moduncertainty(f):
    a = f[0][0][1]
    b = f[0][0][0]
    siga = sqrt(f[0][1][1,1])
    sigb = sqrt(f[0][1][0,0])
    sigab2 = f[0][1][1,0]
    print sigab2*2*a/(a+2*b)**4*2*b

    return (a/(a+2*b),\
            sqrt((2*b/(a+2*b)**2)**2*siga**2+(2*a/(a+2*b)**2)**2*sigb**2+\
            sigab2*2*a/(a+2*b)**4*2*b))

#Plot misalignment fit results
def misalign():
    os.chdir('/Users/rallured/PythonWorkspace/Calibration/121031/')
    c = genfromtxt('CoarseAlign6.txt')
    c = reshape(c,(31,11))
    c = transpose(c)

    clf()
    mycontour(c,x=arange(31)*.05-.75,y=arange(11)*.05)
    title('Fit Results for Misalignment Scan')
    xlabel(r'$\beta$ (deg)')
    ylabel(r'$\alpha$ (deg)')

#Double reflected search - max 0.77 Hz
def dblref():
    os.chdir('/Users/rallured/PythonWorkspace/Calibration/121023/')
    d = genfromtxt('BeamSearch2.txt')
    clf()
    plot(arange(100)*.0464/2,d[0]/d[1])
    title('10/23/12 Alignment Scan')
    xlabel('Relative Polarizer Incidence Angle (deg)')
    ylabel('Normalized Counts')

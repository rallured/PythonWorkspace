import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from os import *
import pdb

#Read in data from June 2011 measurements
chdir('/Users/rallured/IDLWorkspace82/Multilayer/LLNL/m1-110525')
osin = np.transpose(np.genfromtxt('pat037319.abs',skip_header=1))
osout = np.transpose(np.genfromtxt('pat037312.abs',skip_header=1))
correction = osin[1]/osout[1]
center = np.transpose(np.genfromtxt('pat037272.abs',skip_header=1))
up = np.transpose(np.genfromtxt('pat037276.abs',skip_header=1))
down = np.transpose(np.genfromtxt('pat037277.abs',skip_header=1))
right = np.transpose(np.genfromtxt('pat037278.abs',skip_header=1))
left = np.transpose(np.genfromtxt('pat037279.abs',skip_header=1))
rmodel = np.transpose(np.genfromtxt('RModel.txt',skip_header=20))
d450 = np.transpose(np.genfromtxt('pat037319.abs',skip_header=1))
d432 = np.transpose(np.genfromtxt('pat037323.abs',skip_header=1))
d436 = np.transpose(np.genfromtxt('pat037324.abs',skip_header=1))
d440 = np.transpose(np.genfromtxt('pat037325.abs',skip_header=1))
d444 = np.transpose(np.genfromtxt('pat037326.abs',skip_header=1))
d448 = np.transpose(np.genfromtxt('pat037327.abs',skip_header=1))
d452 = np.transpose(np.genfromtxt('pat037328.abs',skip_header=1))
d456 = np.transpose(np.genfromtxt('pat037329.abs',skip_header=1))
d460 = np.transpose(np.genfromtxt('pat037330.abs',skip_header=1))
d464 = np.transpose(np.genfromtxt('pat037331.abs',skip_header=1))
d468 = np.transpose(np.genfromtxt('pat037332.abs',skip_header=1))
tmodel = np.transpose(np.genfromtxt('TModel.txt',comments=';'))
t560 = np.transpose(np.genfromtxt('pat037333.abs',skip_header=1))
t180 = np.transpose(np.genfromtxt('pat037342.abs',skip_header=1))
t440 = np.transpose(np.genfromtxt('pat037306.abs',skip_header=1))
t760 = np.transpose(np.genfromtxt('pat037302.abs',skip_header=1))
t900 = np.transpose(np.genfromtxt('pat037303.abs',skip_header=1))
fulltrans = np.transpose(np.genfromtxt('Transmission091111.txt',skip_header=20))

#Read in final fit data
finalfit = np.transpose(np.genfromtxt('m1-110525FinalFit4.txt',comments=';'))
ind = np.where(np.logical_and(finalfit[0]>=490,finalfit[0]<=550))
rmodel = [finalfit[0][ind],finalfit[1][ind]]
ind = np.where(np.logical_and(finalfit[0]>=40,finalfit[0]<=1400))
tmodel = [finalfit[0][ind],finalfit[2][ind]]

#Read in data from December 2011 measurements
chdir('./m1-110525_Dec18_2011')
updec = np.transpose(np.genfromtxt('pat039180.abs',skip_header=1))
downdec = np.transpose(np.genfromtxt('pat039182.abs',skip_header=1))
rightdec = np.transpose(np.genfromtxt('pat039183.abs',skip_header=1))
leftdec = np.transpose(np.genfromtxt('pat039184.abs',skip_header=1))

def analyzerefdiff():
    #Calculate integrated reflectivities
    upint = np.sum(up[1]*correction*.5)
    downint = np.sum(down[1]*correction*.5)
    rightint = np.sum(right[1]*correction*.5)
    leftint = np.sum(left[1]*correction*.5)
    upintdec = np.sum(updec[1]*.5)
    downintdec = np.sum(downdec[1]*.5)
    rightintdec = np.sum(rightdec[1]*.5)
    leftintdec = np.sum(leftdec[1]*.5)
    updiff = (upint-upintdec)/upint
    print 'updiff: '+str(updiff)
    downdiff = (downint-downintdec)/downint
    print 'downdiff: '+str(downdiff)
    rightdiff = (rightint-rightintdec)/rightint
    print 'rightdiff: '+str(rightdiff)
    leftdiff = (leftint-leftintdec)/leftint
    print 'leftdiff: '+str(leftdiff)
    meanint = np.mean([upint,downint,rightint,leftint])
    meanintdec = np.mean([upintdec,downintdec,rightintdec,leftintdec])
    pdb.set_trace()
    print (meanint-meanintdec)/meanint

def makepositionplot():
    #Calculate integrated reflectivities
    centerint = np.sum(center[1]*correction*.5)
    upint = np.sum(up[1]*correction*.5)
    downint = np.sum(down[1]*correction*.5)
    rightint = np.sum(right[1]*correction*.5)
    leftint = np.sum(left[1]*correction*.5)
    meanint = np.mean([centerint,upint,downint,rightint,leftint])
    print meanint

    #Calculate center and width of energy band
    ecen = rmodel[0][np.where(rmodel[1] == np.max(rmodel[1]))]
    diff = np.abs(rmodel[1] - np.max(rmodel[1])/2)
    ewidth = 2*np.abs(ecen - rmodel[0][np.where(diff == np.min(diff))])
    print ecen
    print ewidth

    #Make plot
    plt.ion()
    plt.clf()
    plt.hold(True)
    #matplotlib.rc('text',usetex=True)
    #matplotlib.rcParams['lines.linewidth']=3
    #matplotlib.rcParams['lines.linestyle']='-.'
    #matplotlib.rcParams['font.weight']='bolder'
    plt.plot(center[0],center[1]*correction,'--',label='Center')
    plt.plot(up[0],up[1]*correction,'--',label='Y=+10mm')
    plt.plot(down[0],down[1]*correction,'--',label='Y=-10mm')
    plt.plot(left[0],left[1]*correction,'--',label='X=-10mm')
    plt.plot(right[0],right[1]*correction,'--',label='X=+10mm')
    plt.plot(rmodel[0],rmodel[1],'-',linewidth=1,label='IMD Model')
    plt.xlim([490,550])
    plt.legend()
    plt.title('M1-110525 Multilayer Response')
    plt.xlabel('Energy (eV)')
    plt.ylabel('Reflectance (fractional)')
    #plt.text(492,.015,'Center: '+str(ecen))
    #plt.text(492,.014,'FWHM: '+str(ewidth))
##    plt.savefig('PositionPlotShape.eps')

def makeangleplot():
    #Make plot
    plt.ion()
    plt.clf()
    plt.hold(True)
    matplotlib.rc('text',usetex=True)
    matplotlib.rcParams['lines.linewidth']=3
    matplotlib.rcParams['lines.linestyle']='-.'
    matplotlib.rcParams['font.weight']='bolder'
    plt.plot(d432[0],d432[1],label='theta = 43.2')
    plt.plot(d436[0],d436[1],label='theta = 43.6')
    plt.plot(d440[0],d440[1],label='theta = 44.0')
    plt.plot(d444[0],d444[1],label='theta = 44.4')
    plt.plot(d448[0],d448[1],label='theta = 44.8')
    plt.plot(d452[0],d452[1],label='theta = 45.2')
    plt.plot(d456[0],d456[1],label='theta = 45.6')
    plt.plot(d460[0],d460[1],label='theta = 46.0')
    plt.plot(d464[0],d464[1],label='theta = 46.4')
    plt.plot(d468[0],d468[1],label='theta = 46.8')
    plt.xlim([490,560])
    plt.legend()
    plt.title('Reflectance Curves at Various Incidence Angles')
    plt.xlabel('Photon Energy (eV)')
    plt.ylabel('Reflectance (fractional)')
    plt.savefig('AnglePlot.eps')

def maketransmissionplot():
    #Make plot
    plt.ion()
    plt.clf()
    plt.hold(True)
##    matplotlib.rc('text',usetex=True)
##    matplotlib.rcParams['lines.linewidth']=3
##    matplotlib.rcParams['lines.linestyle']='.'
##    matplotlib.rcParams['font.weight']='bolder'
    matplotlib.rcParams['legend.fontsize']=14
    plt.semilogy(t180[0],t180[1],'.',label='Measured Data (180-284 eV)')
    plt.plot(t440[0],t440[1],'.',label='Measured Data (440-574 eV)')
    plt.plot(t560[0][0:-2],t560[1][0:-2],'.',label='Measured Data (560-778 eV)')
    plt.plot(t760[0],t760[1],'.',label='Measured Data (760-932 eV)')
    plt.plot(t900[0][0:-1],t900[1][0:-1],'.',label='Measured Data (900-1304 eV)')
    plt.plot(tmodel[0],tmodel[1],'--',label='IMD Model')
    plt.xlim([0,1500])
    plt.ylim([10**-4,1])
    plt.legend(loc='lower right',numpoints=1)
    plt.title('M1-110525 Transmission Measurements')
    plt.xlabel('Energy (eV)')
    plt.ylabel('Transmission (fractional)')
    plt.xlim([0,1.6e3])
    plt.savefig('TransPlot.eps')
    matplotlib.rcParams['legend.fontsize']=16

def makebeforeafterplot():
    #Read in data
    chdir('/Users/ryanallured/IDLWorkspace80/Multilayer/LLNL/m1-110525')
    centerbefore = np.transpose(np.genfromtxt('pat037272.abs',skip_header=1))
    centerafter = np.transpose(np.genfromtxt('pat037319.abs',skip_header=1))
    upbefore = np.transpose(np.genfromtxt('pat037276.abs',skip_header=1))
    upafter = np.transpose(np.genfromtxt('pat037316.abs',skip_header=1))

    #Correct data for order sorter
    correction = normalize()
    centerbefore[1] = centerbefore[1]*correction
    upbefore[1] = upbefore[1]*correction
    upafter[1] = upafter[1]*correction

    #Make plot
    plt.ion()
    plt.clf()
    plt.hold(True)
    matplotlib.rc('text',usetex=True)
    matplotlib.rcParams['lines.linewidth']=1
    matplotlib.rcParams['lines.linestyle']='-'
    matplotlib.rcParams['font.weight']='bolder'
    plt.plot(centerbefore[0],centerbefore[1],label='Center, Before')
    plt.plot(centerafter[0],centerafter[1],label='Center, After')
    plt.plot(upbefore[0],upbefore[1],label='Y=+10mm, Before')
    plt.plot(upafter[0],upafter[1],label='Y=+10mm, After')
    plt.xlim([490,550])
    plt.legend()
    plt.title('Reflectance Curves Before and After Thermal Cycling')
    plt.xlabel('Photon Energy (eV)')
    plt.ylabel('Reflectance (\%)')
    plt.savefig('ThermalPlot.eps')

def makefullbrptrans():
    #Make plot
    plt.ion()
    plt.clf()
    plt.hold(True)
    matplotlib.rc('text',usetex=True)
    matplotlib.rcParams['lines.linewidth']=1
    plt.plot(fulltrans[0],100*trans[1])
    ind = np.where(fulltrans[0] == 2700)
    plt.plot(fulltrans[0][ind],100*trans[1][ind],'*')
    plt.text(3000,74,'75.2\% at 2.7 keV')
    plt.title('Predicted BRP Multilayer Reflector Transmission')
    plt.xlabel('Energy (eV)')
    plt.ylabel('Transmission (\%)')
    plt.savefig('FullBRPTrans.eps')

from plotting import *
import pdb

#Make pretty optimization plots
def optplots():
    #Make Al2O3/V plot
    chdir('/Users/rallured/IDLWorkspace82/Multilayer/')
    alopt = genfromtxt('AlOptData.txt')
    alopt = reshape(alopt,(31,21))
    alopt = transpose(alopt)
    gam = arange(.4,.61,.01)
    ener = arange(495,526)
    clf()
    mycontour(alopt*100,x=ener,y=gam,fmt='%3.1f')
    plot([506],[.52],'ro')
    plot([517],[.54],'ro')
    title(r'Al$_2$O$_3$/V MDP',verticalalignment='bottom')
    xlabel('Peak Reflectance Energy (eV)')
    ylabel('Gamma')
    savefig('AlOpt.eps')
    pdb.set_trace()

    wcopt = genfromtxt('WCOptData.txt')
    wcopt = reshape(wcopt,(21,41))
    wcopt = transpose(wcopt)
    inc = arange(45,24.5,-.5)
    ener = arange(450,555,5)
    clf()
    mycontour(wcopt*100,x=ener,y=inc)
    plot(ener[11],inc[19],'ro')
    title('WC/SiC MDP')
    xlabel('Peak Reflectance Energy (eV)')
    ylabel('Glancing Angle (deg)')
    savefig('WCOpt.eps')
    pdb.set_trace()

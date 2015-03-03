from os import *
import numpy as np
from scipy import *
from string import *
import glob
import matplotlib.pyplot as plt
import pyfits

#Reads in the spectral results from a file from each individual observation,
#these were saved in the command line using np.savetxt as one big file
def getstddata(proposalid):
    chdir('/Users/ryanallured/RXTE/MAXI_J1543-564/results/P'+proposalid)
    for path in glob.glob('96*'):
        chdir(path)
        temp = pyfits.open('top_std2.pha')
        exptime = temp[1].header['exposure']
        time = (temp[0].header['tstart']+temp[0].header['tstop'])/2
        time = 49353.0+time/86400.
        try: data
        except NameError:
            data = [np.genfromtxt('stdmodel.out',dtype='float')]
            print type(data)
            tstart = array(time)
            obsid = [path]
            ex = [exptime]
            print ex
        else:
            data.append(np.genfromtxt('stdmodel.out',dtype='float'))
            tstart = append(tstart,time)
            obsid.append(path)
            ex.append(exptime)
        chdir('..')

    return [data,tstart,obsid,ex]

#Combine the spectral results from the two proposal IDs 96371 and 96430
def combineresults():
    blah1 = getstddata('96371')
    blah2 = getstddata('96430')

    spectral1 = array(blah1[0])
    spectral2 = array(blah2[0])
    spectral = concatenate((spectral1,spectral2))

    tstart1 = array(blah1[1])
    tstart2 = array(blah2[1])
    tstart = append(tstart1,tstart2)
    
    obsid1 = array(blah1[2])
    obsid2 = array(blah2[2])
    obsid = append(obsid1,obsid2)

    spectral = spectral[argsort(tstart)]
    obsid = obsid[argsort(tstart)]
    tstart = sort(tstart)

    return [spectral,tstart,obsid]

#Reads in the spectral results that were consolidated with the previous routine into one file
def spectralresults():
    chdir('/Users/ryanallured/RXTE/MAXI_J1543-564/results/')
    spectral = np.genfromtxt('SpectralResults3.txt',names='mflux,mfluxerr1,mfluxerr2,chisq,nh,tin, \
        tinerr1,tinerr2,phoind,phoinderr1,phoinderr2,pnorm,pnormerr1,pnormerr2,dnorm, \
        dnormerr1,dnormerr2,gsigma,gsigerr1,gsigerr2,gnorm,gnormerr1,gnormerr2,unabsflux,diskflux,sflux,hflux')
    return spectral

#Reads in the timing analysis results from the power density spectra
#RMS -> integrated power in PDS from 0.1 - 10 Hz
#QPOFreq -> location of the peak in the PDS
#Qpopow -> integrated power in the QPO (quasi-periodic oscillation)
#Err -> Error in the qpo freq (not used)
def timingresults():
    chdir('/Users/ryanallured/RXTE/MAXI_J1543-564/results/')
    timing = np.genfromtxt('TimingResults2.txt',names='rms,qpofreq,qpopow,err',\
                           skip_header=1)
    return timing

#Reads in the start times for each observation
def starttimes():
    chdir('/Users/ryanallured/RXTE/MAXI_J1543-564/results/')
    tstart = np.genfromtxt('Tstart.txt')
    return tstart

#Reads in the sorted ObsIDs
def obsids():
    chdir('/Users/ryanallured/RXTE/MAXI_J1543-564/results/')
    obsid = np.genfromtxt('ObsIDs2.txt',dtype='string')
    return obsid

#Makes an ASCII file of all the locations of the good time interval files
def makegtiindex():
    chdir('/Users/ryanallured/RXTE/MAXI_J1543-564/results/P96371')

    for path in glob.glob('96371-??-??-??'):
        chdir(path)
        try: gti
        except NameError:
            gti = [getcwd() + '/xte.gti']
        else:
            gti.append(getcwd()+'/xte.gti')
        chdir('..')

    chdir('..')
    savetxt
    
    return gti

#Classifies each observation into a nominal accretion state using spectral and timing parameters
def classifystates(spectral,timing):
    disk = spectral['diskflux']/spectral['unabsflux']

    #Thermal states
    thermal = np.logical_and(disk > .75,timing['qpopow'] < .005)
    thermal = np.logical_and(thermal,timing['rms'] < .075)

    #Hard states
    hard = np.logical_and(disk < .2,spectral['phoind'] > 1.4)
    hard = np.logical_and(hard,spectral['phoind'] < 2.1)
    hard = np.logical_and(hard,timing['rms'] > .1)

    #Steep Power Law states
    spl = np.logical_and(spectral['phoind'] > 2.4,np.logical_and(timing['rms'] < .15,\
        np.logical_or(np.logical_and(disk < .8,timing['qpopow'] > .01),disk < .5)))

    #Intermediate states
    inter = np.logical_not(np.logical_or(thermal,np.logical_or(spl,hard)))

    return [thermal,hard,spl,inter]

#Plots the results and saves them to a file
def makeplots(spectral,timing,states,tstart):
    #Define secondary spectral parameters
    disk = spectral['diskflux']/spectral['unabsflux']
    hardr = spectral['hflux']/spectral['sflux']
    inters = states[3]

    #Change to plot directory
    chdir('/Users/ryanallured/RXTE/MAXI_J1543-564/results/Plots')

    #Make HID
    plt.clf()
    plt.hold(True)
    plt.plot(hardr[states[0]],spectral['mflux'][states[0]],'r.',label='Thermal',markersize=20)
    plt.plot(hardr[states[1]],spectral['mflux'][states[1]],'b.',label='Hard',markersize=20)
    plt.plot(hardr[states[2]],spectral['mflux'][states[2]],'g.',label='SPL',markersize=20)
    plt.plot(hardr[states[3]],spectral['mflux'][states[3]],'y.',label='Intermediate',markersize=20)
    plt.legend()
    plt.title('Intensity vs. Hardness')
    plt.xlabel('Hardness Ratio (6.3keV-10.5keV)/(3.8keV-7.3keV)')
    plt.ylabel('Model Flux (cts/s/cm^2)')
    plt.savefig('HID.png')

    #Make lightcurve
    plt.clf()
    merr1 = spectral['mflux']-spectral['mfluxerr1']
    merr2 = spectral['mfluxerr2']-spectral['mflux']
    plt.plot(tstart[states[0]],spectral['mflux'][states[0]],'r.',label='Thermal',markersize=20)
    plt.plot(tstart[states[1]],spectral['mflux'][states[1]],'b.',label='Hard',markersize=20)
    plt.plot(tstart[states[2]],spectral['mflux'][states[2]],'g.',label='SPL',markersize=20)
    plt.plot(tstart[states[3]],spectral['mflux'][states[3]],'y.',label='Intermediate',markersize=20)
    plt.legend()
    plt.title('Lightcurve for MAXI J1543-564')
    plt.xlabel('MJD')
    plt.ylabel('Model Flux (cts/s/cm^2')
    plt.savefig('Lightcurve.png')

    #Make power-law flux vs. disk flux
    plt.clf()
    powerlaw = spectral['unabsflux'] - spectral['diskflux']
    plt.plot(spectral['diskflux'][states[0]],powerlaw[states[0]],'.r',label='Thermal',markersize=20)
    plt.plot(spectral['diskflux'][states[1]],powerlaw[states[1]],'.b',label='Hard',markersize=20)
    plt.plot(spectral['diskflux'][states[2]],powerlaw[states[2]],'.g',label='SPL',markersize=20)
    plt.plot(spectral['diskflux'][states[3]],powerlaw[states[3]],'.y',label='Intermediate',markersize=20)
    plt.title('Power Law Flux vs. Disk Flux')
    plt.xlabel('Disk Flux (cts/s/cm^2)')
    plt.ylabel('Power Law Flux (cts/s/cm^2)')
    plt.savefig('PlawvsDflux.png')

    #Make disk fraction plot
    plt.clf()
    plt.plot(tstart[states[0]],disk[states[0]],'.r',label='Thermal',markersize=20)
    plt.plot(tstart[states[1]],disk[states[1]],'.b',label='Hard',markersize=20)
    plt.plot(tstart[states[2]],disk[states[2]],'.g',label='SPL',markersize=20)
    plt.plot(tstart[states[3]],disk[states[3]],'.y',label='Intermediate',markersize=20)
    plt.title('Disk Fraction vs. Observation Time')
    plt.xlabel('MJD')
    plt.ylabel('Disk Fraction')
    plt.savefig('Disk.png')

    plt.clf()
    #Make Photon Index vs. Hardness
    plt.subplot(3,1,1)
    plt.plot(hardr[states[0]],spectral['phoind'][states[0]],'.r',label='Thermal',markersize=20)
    plt.plot(hardr[states[1]],spectral['phoind'][states[1]],'.b',label='Hard',markersize=20)
    plt.plot(hardr[states[2]],spectral['phoind'][states[2]],'.g',label='SPL',markersize=20)
    plt.plot(hardr[states[3]],spectral['phoind'][states[3]],'.y',label='Intermediate',markersize=20)
    plt.title('Photon Index, Disk Fraction, RMS Power vs. Hardness Ratio')
    plt.ylabel('Photon Index')
    #Make Disk fraction vs. hardness
    plt.subplot(3,1,2)
    plt.plot(hardr[states[0]],disk[states[0]],'.r',label='Thermal',markersize=20)
    plt.plot(hardr[states[1]],disk[states[1]],'.b',label='Hard',markersize=20)
    plt.plot(hardr[states[2]],disk[states[2]],'.g',label='SPL',markersize=20)
    plt.plot(hardr[states[3]],disk[states[3]],'.y',label='Intermediate',markersize=20)
    plt.ylabel('Disk Fraction')
    #Make RMS power vs hardness
    plt.subplot(3,1,3)
    plt.plot(hardr[states[0]],timing['rms'][states[0]],'.r',label='Thermal',markersize=20)
    plt.plot(hardr[states[1]],timing['rms'][states[1]],'.b',label='Hard',markersize=20)
    plt.plot(hardr[states[2]],timing['rms'][states[2]],'.g',label='SPL',markersize=20)
    plt.plot(hardr[states[3]],timing['rms'][states[3]],'.y',label='Intermediate',markersize=20)
    plt.ylabel('RMS Power')
    plt.xlabel('Hardness Ratio (6.3keV-10.5keV)/(3.8keV-7.3keV)')
    plt.savefig('SpectralvsHard.png')

    #Make Photon Index vs. Tstart
    plt.clf()
    plt.plot(tstart[states[0]],spectral['phoind'][states[0]],'.r',label='Thermal',markersize=20)
    plt.plot(tstart[states[1]],spectral['phoind'][states[1]],'.b',label='Hard',markersize=20)
    plt.plot(tstart[states[2]],spectral['phoind'][states[2]],'.g',label='SPL',markersize=20)
    plt.plot(tstart[states[3]],spectral['phoind'][states[3]],'.y',label='Intermediate',markersize=20)
    plt.title('Photon Index vs. Time')
    plt.ylabel('Photon Index')
    plt.xlabel('MJD')
    plt.savefig('PhoIndvsTime.png')

    #Make RMS power vs Tstart
    plt.clf()
    plt.plot(tstart[states[0]],timing['rms'][states[0]],'.r',label='Thermal',markersize=20)
    plt.plot(tstart[states[1]],timing['rms'][states[1]],'.b',label='Hard',markersize=20)
    plt.plot(tstart[states[2]],timing['rms'][states[2]],'.g',label='SPL',markersize=20)
    plt.plot(tstart[states[3]],timing['rms'][states[3]],'.y',label='Intermediate',markersize=20)
    plt.title('RMS Power vs. Time')
    plt.ylabel('RMS Power')
    plt.xlabel('MJD')
    plt.savefig('RMSvsTime.png')

    #Make RMS power vs Tstart with errors
    plt.clf()
    stddev = sqrt(timing['err'])
    plt.errorbar(tstart,timing['rms'],stddev,stddev)
    plt.title('RMS Power vs. Time')
    plt.ylabel('RMS Power')
    plt.xlabel('MJD')
    plt.savefig('RMSErr.png')

    #Make QPO power vs Tstart
    plt.clf()
    plt.plot(tstart[states[0]],timing['qpopow'][states[0]],'.r',label='Thermal',markersize=20)
    plt.plot(tstart[states[1]],timing['qpopow'][states[1]],'.b',label='Hard',markersize=20)
    plt.plot(tstart[states[2]],timing['qpopow'][states[2]],'.g',label='SPL',markersize=20)
    plt.plot(tstart[states[3]],timing['qpopow'][states[3]],'.y',label='Intermediate',markersize=20)
    plt.title('QPO Power vs. Time')
    plt.ylabel('QPO Power')
    plt.xlabel('MJD')
    plt.savefig('QPOvsTime.png')

    #Make Chisq vs Tstart
    plt.clf()
    plt.plot(tstart,spectral['chisq'])
    plt.title('Chi Squared vs. Time')
    plt.ylabel('Chi Squared')
    plt.ylabel('MJD')
    plt.savefig('Chisq.png')

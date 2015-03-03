from numpy import *
from matplotlib.pyplot import *
import pdb

#Read in Zygo ASCII file
def readzygo(filename):
    #Open file
    f = open(filename,'r')

    #Read third line to get intensity shape
    for i in range(3):
        l = f.readline()
    l = l.split(' ')
    iwidth = int(l[2])
    iheight = int(l[3])

    #Read fourth line to get phase shape
    l = f.readline()
    l = l.split(' ')
    pwidth = int(l[2])
    pheight = int(l[3])

    #Read eighth line to get scale factors
    for i in range(4):
        l = f.readline()
    l = l.split(' ')
    scale = float(l[1])
    wave = float(l[2])
    o = float(l[4])
    latscale = float(l[6])

    #Read eleventh line to get phase resolution
    f.readline()
    f.readline()
    l = f.readline()
    l = l.split(' ')
    phaseres = l[0]
    if phaseres is 0:
        phaseres = 4096
    else:
        phaseres = 32768

    #Read through to first '#' to signify intensity
    while (l[0]!='#'):
        l = f.readline()

    #Read intensity array
    #If no intensity, l will be '#' below
    l = f.readline()
    while (l[0]!='#'):
        #Convert to array of floats
        l = array(l.split(' '))
        l = l[:-1].astype('float')
        #Merge into intensity array
        try:
            intensity = concatenate((intensity,l))
        except:
            intensity = l
        #Read next line
        l = f.readline()

    #Reshape into proper array
    try:
        intensity = reshape(intensity,(iheight,iwidth))
    except:
        intensity = NaN

    #Read phase array
    l = f.readline()
    while (l!=''):
        #Convert to array of floats
        l = array(l.split(' '))
        l = l[:-1].astype('float')
        #Merge into intensity array
        try:
            phase = concatenate((phase,l))
        except:
            phase = l
        #Read next line
        l = f.readline()

    phase = reshape(phase,(pheight,pwidth))
    phase[where(phase==phase.max())] = nan
    phase = phase*scale*o*wave/phaseres
    f.close()
    print wave, scale, o, phaseres

    return intensity, phase, latscale

#Convert Zygo ASCII to easily readable ASCII format
def convertzygo(filename):
    #read in zygo data
    intensity,phase,latscale = readzygo(filename)

    savetxt(filename.split('.')[0]+'.txt',phase,header='Lat scale: '+\
            str(latscale)+'\n'+'Units: meters')

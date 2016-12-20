import numpy as np
import matplotlib.pyplot as plt
import reflectivity as ref
from scipy.optimize import curve_fit
import pdb

#topThick = np.linspace(10.,100.,21)
#botThick = np.linspace(100.,700.,101)

def computeFieldModulation(filename,topThick,botThick):
    """
    Load in reflectivity array as computed in IMD
    For each top coating thickness, compute modulation
    and therefore intensity ratio from front and back
    reflections.
    topThick should be vector specifying top thickness
    botThick should be vector specifying bottom (substrate) thickness
    """
    #Load in data
    d = np.genfromtxt(filename)
    d = d.reshape((len(topThick),len(botThick)))

##    #Compute max and min of modulation curve
##    #Fit cos**2 distribution using curve_fit module
##    cos2 = lambda botThick,m,amplitude,phase,freq: m+amplitude*\
##           np.cos(2*np.pi*freq*botThick+phase)**2
##    fit = [curve_fit(cos2,botThick,di,\
##                    p0=[np.min(di),np.max(di)-np.min(di),\
##                        0.,1./(632.8/1.5)]) for di in d]
##
##    #Test out validity of fit
##    plt.plot(botThick,d[0])
##    diff = [d[i]-cos2(botThick,*fit[i][0]) for i in range(11)]
##    pdb.set_trace()

    #Use min/max and compute intensity ratios of front and back contributions
    eRatio = []
    for i in range(len(topThick)):
        eSum = np.sqrt(max(d[i]))
        eDiff = np.sqrt(min(d[i]))
        eFront = (eSum+eDiff)/2.
        eBack = (eSum-eDiff)/2.
        eRatio.append((eBack/eFront)**2)

    return eRatio

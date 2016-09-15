import numpy as np
import matplotlib.pyplot as plt
import utilities.imaging.man as man
import utilities.imaging.analysis as anal
import utilities.imaging.fitting as fit

def repeatability(d1,d2):
    """
    
    """
    #Get zoomed in view
    d1sub = man.stripnans(d1)
    d2sub = man.stripnans(d2)
    #Perform misalignment fits
    res1 = fit.fitCylMisalign(d1sub)
    res2 = fit.fitCylMisalign(d2sub)
    #Remove misalignment
    d1sub = d1sub - res1[0]
    d2sub = d2sub - res2[0]
    #Set original images based on filtered subapertures
    d1[~np.isnan(d1)] = d1sub[~np.isnan(d1sub)]
    d2[~np.isnan(d2)] = d2sub[~np.isnan(d2sub)]
    return d1,d2

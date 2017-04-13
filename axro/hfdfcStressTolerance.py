import numpy as np
import matplotlib.pyplot as plt
import axro.solver as slv
import utilities.imaging.man as man
import legendremod as leg

uniIF = '/home/rallured/Dropbox/AXRO/InfluenceFunctions/CoatingStress/2D/170216_BlockUniformIF.fits'

def sagExplanation(minmode='slope',pzt_min_strain=0.,pzt_max_strain=100.,\
                   d=None,azweight=.015):
    """
    Create uniform sag distortion
    Run 5mm uniform stress IFs on quadratic sag distortion map
    Plot correction, residual, and stress map
    """
    slv.setupAlexey(uniIF,.5,pzt_min_strain=pzt_min_strain,\
                    pzt_max_strain=pzt_max_strain,\
                    minmode=minmode,azim_slope_weight=azweight)
    slv.preMath()

    if d is None:
        x,y = man.autoGrid(np.zeros((200,200)))
        d = leg.singleorder(x,y,0,2)
    else:
        d = man.newGridSize(d,(200,200))
    
    res = slv.pyExecute(d)
    res2 = slv.pyExecute(-d)

    fig = plt.figure()
    fig.add_subplot(2,2,1)
    p = plt.imshow(res[0])
    plt.colorbar()
    p.axes.set_xticks([])
    p.axes.set_yticks([])
    plt.title('Closest Match')
    fig.add_subplot(2,2,2)
    p = plt.imshow(man.remove2DLeg(res[0],xo=10))
    plt.colorbar()
    p.axes.set_xticks([])
    p.axes.set_yticks([])
    plt.title('Flattened Match')
    fig.add_subplot(2,2,3)
    p = plt.imshow(man.remove2DLeg(res[1],xo=10))
    plt.colorbar()
    p.axes.set_xticks([])
    p.axes.set_yticks([])
    plt.title('Flattened Residual')
    fig.add_subplot(2,2,4)
    p = plt.imshow(res[2].reshape((21,21)))
    plt.colorbar()
    p.axes.set_xticks([])
    p.axes.set_yticks([])
    plt.title('Stress Map')

    return res,res2

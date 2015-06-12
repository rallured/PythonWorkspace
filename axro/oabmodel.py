import numpy as np
import pyfits
from scipy.interpolate import griddata
import pdb
import traces.conicsolve as con
import traces.PyTrace as PT

#Global directory variables for problem
datadir = '/home/rallured/data/solve_pzt/'

def convertOABInfluence(filename,Nx,Ny,method='cubic'):
    """Read in Vanessa's CSV file for 1 m RoC AXRO mirror
    Convert the azimuthal cartesian coordinate into an arc length
    Regrid the perturbed coordinates onto a regular arc and axial
    coordinate grid
    """
    d = np.transpose(np.genfromtxt(filename,skip_header=1,delimiter=','))
##    x0,y0,z0 = pyfits.getdata('150519StartingNodes.fits')
##    x = x0 + d[5]
##    y = y0 + d[6]
##    z = z0 + d[7]
    x = d[2] + d[5]
    y = d[3] + d[6]
    z = d[4] + d[7]
    pdb.set_trace()

    theta = np.arctan2(y,-z)
    rPert = 1. - np.sqrt(x**2+z**2) #Positive result is toward OA
    arc = theta * 1.

    #Now have x for cylindrical axis, arc for azimuthal axis,
    #and rPert for the radial perturbation
    #Regrid to a regular grid and you're done
    gy = np.linspace(x.min(),x.max(),Nx+2)
    gx = np.linspace(arc.min(),arc.max(),Ny+2)
    gx,gy = np.meshgrid(gx,gy)
    d = np.transpose(griddata((x,arc),rPert,(gy,gx),method=method))
    d[np.isnan(d)] = 0.

    return d[1:-1,1:-1]

def oabDistortion(amp,freq,phase,filename):
    """Similar to the original createDistortion, this introduces
    a sinusoidal ripple to the 1 meter RoC geometry.
    """
    #Create distortion array
    y,x = np.mgrid[0:200,0:200]
    x,y = x*.5,y*.5
    d = amp*np.sin(2*np.pi*freq*y+phase)

    #Save as fits file
    pyfits.writeto(datadir+'distortions/'+filename,d,clobber=True)

    return d

def wolterDistortion(filename,Nx,Ny):
    """This creates the distortion map of a 1 meter RoC
    Wolter I as compared to a 1 meter cylinder
    Should be matched to height/theta grid
    Of influence functions
    """
    #Determine Wolter prescription by finding
    #intersection radius closest to 1 meter
    #at 1125 mm above focal point
    rint = np.linspace(990.,1000.,1000)
    r0 = con.primrad(10125.,rint,10000.)
    rselect = rint[np.argmin(abs(r0-1000.))]
    a,p,d,e = con.woltparam(rselect,10000.)
    
    #Need to get rays onto Wolter I surface
    tspan = .04952303577 #maximum theta from IFs
    zspan = 50.*np.cos(2*a)
    t,z = np.meshgrid(np.linspace(-tspan,tspan,200),\
                      np.linspace(10125.-zspan,10125.+zspan,200))
    #Set up source rays along optical axis
    PT.z = z.flatten()
    PT.l = np.cos(t.flatten())
    PT.m = np.sin(t.flatten())
    PT.n = np.zeros(np.size(PT.z))
    PT.x = PT.l*1000.
    PT.y = PT.m*1000.
    PT.ux = np.zeros(np.size(PT.z))
    PT.uy = np.zeros(np.size(PT.z))
    PT.uz = np.zeros(np.size(PT.z))
    #Trace to Wolter I surface and go to cylinder
    PT.wolterprimary(rselect,10000.)
    PT.transform(con.primrad(10125.,rselect,10000.),0,10125.,0,0,0)
    PT.transform(0,0,0,0,a,0)
    PT.transform(-1000.,0,0,0,0,0)
    pdb.set_trace()

    #Compute radial perturbations, reshape, and return distortion
    rpert = 1000.-np.sqrt(PT.x**2+PT.y**2)
    return rpert.reshape((200,200)),z.reshape((200,200)),t.reshape((200,200))

def analyzeOABCorrection():
    """Compute the correction to the OAB cylinder and
    analyze the resulting performance. Should compute
    a resulting RMS diameter assuming equal and uncorrelated
    performance with the secondary.
    In reality, the performance *will* likely be correlated
    due to sag addition.
    """
    

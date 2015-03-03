from numpy import *
from matplotlib.pyplot import *
import os
from scipy.optimize import *
from gaussfitter import *

os.chdir('/Users/ryanallured/Documents/Research/GEMS/Fluorescence/')
cnts = genfromtxt('120410FluorescenceTest4500s.mca',\
                  skip_header=12,skip_footer=1)
chan = arange(1,513)

#Initial fit to get Al K energy
fit = onedgaussfit(chan[42:63],cnts[42:63],err=1+sqrt(cnts[42:63]+.75),\
                   params=[0,58,49,50],fixed=[True,False,False,False])
al = fit[0][2]
alsig = fit[0][3]
sigconst = alsig/sqrt(al)
conv = 103.2/(5.9-3.2)
esc_mnka = 103.2#(5.9-3.2)*conv
esc_mnkasig = sigconst*sqrt(esc_mnka)
esc_mnkb = (6.49-3.2)*conv
esc_mnkbsig = sigconst*sqrt(esc_mnkb)
mnka = 204.
mnkasig = sigconst*sqrt(mnka)
mnkb = 204*6.49/5.9
mnkbsig = sigconst*sqrt(mnkb)

def specshape(x,a0,a1,a2,a3,a4,w0,w1,w2,w3,w4,e0,e1,e2,e3,e4):
    return onedgaussian(x,0,a0,e0,w0)+onedgaussian(x,0,a1,e1,w1)+\
           onedgaussian(x,0,a2,e2,w2)+\
           onedgaussian(x,0,a3,e3,w3)+onedgaussian(x,0,a4,e4,w4)

guess = [37,155,16,1250,125,15,20,20,30,30,98,190,200,355,375]
guess = real(guess)
guess = guess.astype('float')

fit2 = curve_fit(specshape,chan[88:],cnts[88:],\
                p0=guess,sigma=1+sqrt(cnts[88:]+.75))

clf()
plot(chan,cnts)
plot(chan[88:],specshape(chan[88:],*fit2[0]))

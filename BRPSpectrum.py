from numpy import *
from matplotlib.pyplot import *
from gaussfitter import *
from os import *

chdir('/Users/ryanallured/Documents/Research/GEMS/AnalysisDocuments/PredictedSpectrum')


vnoise = 50.
vsat = 4000
field = .2
voltages = arange(0,5001,1)

c = 25.720503401310726
brpres = c/sqrt(515.)
alkres = c/sqrt(1500.)
brpeff = sqrt(brpres**2+field**2)
aleff = sqrt(alkres**2+field**2)

#Minimum Gain
vbrp = vnoise/(1-brpeff/2)
gmin = vbrp/.83E-3/(515./26.2)
val = vbrp*(1500./515.)
brpgauss = onedgaussian(voltages,0,1,vbrp,brpeff*vbrp/2.35)
algauss = onedgaussian(voltages,0,1,val,aleff*val/2.35)
ion()
clf()
plot(voltages,brpgauss,label='515 eV')
plot(voltages,algauss,label='1500 eV')
plot(voltages,brpgauss+algauss)
plot([vnoise,vnoise],[0,2],'k--')
plot([vsat,vsat],[0,2],'k--')
ylim([0,1.2])
legend()
title('Predicted Spectrum: G='+str(gmin))
xlabel('Voltage (V)')
ylabel('Dimensionless Flux')
savefig('MinimumGain2.png')

#Maximum Gain
val = vsat/(1+aleff/2)
gmax = val/.83E-3/(1500/26.2)
vbrp = val*(515./1500.)
brpgauss = onedgaussian(voltages,0,1,vbrp,brpeff*vbrp/2.35)
algauss = onedgaussian(voltages,0,1,val,aleff*val/2.35)
clf()
plot(voltages,brpgauss,label='515 eV')
plot(voltages,algauss,label='1500 eV')
plot(voltages,brpgauss+algauss)
plot([vnoise,vnoise],[0,2],'k--')
plot([vsat,vsat],[0,2],'k--')
ylim([0,1.2])
legend()
title('Predicted Spectrum: G='+str(gmax))
xlabel('Voltage (V)')
ylabel('Dimensionless Flux')
savefig('MaximumGain2.png')

#Optimal Gain
gopt = gmin+(gmax-gmin)/2.
gopt2 = (vsat+vnoise)/(((1-.5*brpeff)*515.+(1+.5*aleff)*1500.)*(.83E-3/26.2))
vbrp = .83E-3*(515./26.2)*gopt2
val = (1500./515.)*vbrp
brpgauss = onedgaussian(voltages,0,1,vbrp,brpeff*vbrp/2.35)
algauss = onedgaussian(voltages,0,1,val,aleff*val/2.35)
clf()
plot(voltages,brpgauss,label='515 eV')
plot(voltages,algauss,label='1500 eV')
plot(voltages,brpgauss+algauss)
plot([vnoise,vnoise],[0,2],'k--')
plot([vsat,vsat],[0,2],'k--')
ylim([0,1.2])
legend()
title('Predicted Spectrum: G='+str(gopt2))
xlabel('Voltage (V)')
ylabel('Dimensionless Flux')
savefig('OptGain2.png')

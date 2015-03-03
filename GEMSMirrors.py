#from numpy import *
#from matplotlib.pyplot import *
import os

ion()
os.chdir('/Users/ryanallured/IDLWorkspace80/GEMSMirrors/')
data2 = genfromtxt('AuReflectance.txt',skip_header=15)
data2 = transpose(data2)
data10 = genfromtxt('AuReflectance10keV.txt',skip_header=15)
data10 = transpose(data10)
angle = 90-data2[0]
inducedpol10 = (data10[2]-data10[3])/(data10[2]+data10[3])
inducedpol10_2 = (data10[2]**2-data10[3]**2)/(data10[2]**2+data10[3]**2)
inducedpol2 = (data2[2]-data2[3])/(data2[2]+data2[3])
inducedpol2_2 = (data2[2]**2-data2[3]**2)/(data2[2]**2+data2[3]**2)

clf()
semilogy(angle,data2[2],label='2 keV, S Polarization')
plot(angle,data2[3],label='2 keV, P Polarization')
plot(angle,data10[2],label='10 keV, S Polarization')
plot(angle,data10[3],label='10 keV, P Polarization')
legend(loc='upper right')
title('Reflectance from Au')
ylabel('Reflectance (fractional')
xlabel('Incidence Angle (0=glancing)')
savefig('RefAu.eps')

clf()
loglog(angle,inducedpol2,label='2 keV, 1 Reflection')
plot(angle,inducedpol2_2,label='2 keV, 2 Reflections')
plot(angle,inducedpol10,label='10 keV, 1 Reflection')
plot(angle,inducedpol10_2,label='10 keV, 2 Reflections')
legend(loc='lower right')
title('Induced Polarization Due to Au Reflection')
ylabel('Induced Polarization (fractional')
xlabel('Incidence Angle (0=glancing)')
savefig('InducedPolarization.eps')

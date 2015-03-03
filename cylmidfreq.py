from numpy import *
from matplotlib.pyplot import *
import PyTrace as PT
import pdb
import reconstruct
from conicsolve import primrad
import time

#Test out Wolter-I using L-L surfaces
def testwolter():
    tstart = time.time()
    #Set up beam
    PT.subannulus(primrad(8425.,220.,8400.),primrad(8525.,220.,8400.),pi/6,1e6)
    PT.transform(0,0,8800.,pi,0,0)

    #Trace down to primary
    PT.primaryLL(220.,8400.,8525.,8425.,pi/6,[-0.01,.0001],[0,0],[4,0])
    PT.reflect()
    PT.secondaryLL(220.,8400.,8375.,8275.,pi/6,0.001,0,0)
    PT.reflect()
    PT.flat()
    print time.time()-tstart
    

#Trace source to test optic
def tracesource(num,cylNullX=0.,cylNullY=0.,cylNullZ=0.,\
                cylNullRX=0.,cylNullRY=0.,cylNullRZ=0.):
    #Set up source
    PT.circularbeam(12.5,num)

    #Go to cylindrical field lens
    PT.transform(-cylNullX,-cylNullY,-cylNullZ,-cylNullRX,-cylNullRY,-cylNullRZ)
    PT.LJ1653L2()
    PT.transform(cylNullX,cylNullY,cylNullZ,cylNullRX,cylNullRY,cylNullRZ)

    #Go to test optic
    PT.transform(0,0,420.,0,0,0)
    PT.flat()

    return 0

#Trace ideal cylindrical test optic
def idealcylinder(testX=0.,testY=0.,testZ=0.,testRX=0.,testRY=0.,testRZ=0.):
    PT.transform(-testX,-testY,-testZ,-testRX,-testRY,-testRZ)
    PT.transform(0,0,-220.,0,0,0)
    PT.cyl(220.)
    PT.reflect()
    PT.transform(testX,testY,testZ,testRX,testRY,testRZ)
    PT.transform(0,0,0,pi,0,0) #Get rays going in +z direction

    return 0

#Trace ideal Wolter primary
def idealprimary(testX=0.,testY=0.,testZ=0.,testRX=0.,testRY=0.,testRZ=0.):
    PT.transform(-testX,-testY,-testZ,-testRX,-testRY,-testRZ)
    PT.wolterprimtan(220.,8400.)
    PT.reflect()
    PT.transform(testX,testY,testZ,testRX,testRY,testRZ)
    PT.transform(0,0,0,pi,0,0) #Get rays going in +z direction
    PT.transform(0,0,220,0,0,0) #Reference frame at line focus

#Trace rippled Wolter primary
def rippleprimary(amp,freq,testX=0.,testY=0.,testZ=0.,testRX=0.,testRY=0.,testRZ=0.):
    PT.transform(-testX,-testY,-testZ,-testRX,-testRY,-testRZ)
    PT.woltersinetan(220.,8400.,amp,freq)
    PT.reflect()
    PT.transform(testX,testY,testZ,testRX,testRY,testRZ)
    PT.transform(0,0,0,pi,0,0) #Get rays going in +z direction
    PT.transform(0,0,220,0,0,0) #Reference frame at line focus

#Trace using a L-L deformed Wolter primary
def primaryLL(coeff,axial,az,testX=0.,testY=0.,testZ=0.,testRX=0.,testRY=0.,testRZ=0.):
    PT.transform(-testX,-testY,-testZ,-testRX,-testRY,-testRZ)
    PT.primaryLLtan(220.,8400.,8525.,8425.,pi/6.,coeff,axial,az)
    PT.reflect()
    PT.transform(testX,testY,testZ,testRX,testRY,testRZ)
    PT.transform(0,0,0,pi,0,0) #Get rays going in +z direction
    PT.transform(0,0,220,0,0,0) #Reference frame at line focus    

#Trace from test optic - plane is assumed to be at point
#tangent to cylindrical optic, optical axis going in +z direction
def tracefromtest(cylFieldX=0.,cylFieldY=0.,cylFieldZ=0.,\
                  cylFieldRX=0.,cylFieldRY=0.,cylFieldRZ=0.,\
                  cylNullX=0.,cylNullY=0.,cylNullZ=0.,\
                  cylNullRX=0.,cylNullRY=0.,cylNullRZ=0.):
    PT.transform(0,0,200.,0,0,0)#Go back to cylindrical null lens
    PT.transform(cylNullX,cylNullY,cylNullZ,cylNullRX,cylNullRY,cylNullRZ)
    PT.flat()
    PT.LJ1653L2(reverse=True)
    PT.transform(-cylNullX,-cylNullY,-cylNullZ,-cylNullRX,-cylNullRY,-cylNullRZ)
    PT.transform(0,0,100.,0,0,0) #Go to spherical collimator lens
    PT.flat()
    PT.AC254_400_A(reverse=True)
    PT.transform(0,0,331.29,0,0,0) #Go to cylindrical field lens
    PT.flat()
    PT.transform(cylFieldX,cylFieldY,cylFieldZ,cylFieldRX,cylFieldRY,cylFieldRZ)
    PT.LJ1629L2()
    PT.transform(-cylFieldX,-cylFieldY,-cylFieldZ,\
                 -cylFieldRX,-cylFieldRY,-cylFieldRZ)
    PT.transform(0,0,314.250,0,0,0) #Go to field lens
    PT.flat()
    PT.AC508_200_A()
    PT.transform(0,0,161.293,0,0,0) #Go to WFS plane
    PT.flat()
    
#Compute wavefront of perfect primary with respect
#to perfect cylindrical reference - should see the cone
def computeidealinf(num):
    #Trace reference wavefront
    tracesource(num)
    idealcylinder()
    tracefromtest()

    #Subtract all rays outside of WFS
    ind = where(logical_and(abs(PT.y) < 7.3,abs(PT.x) < 7.3))
    PT.vignette(ind=ind)

    #Create reference WFS data
    xang, yang, phase = reconstruct.southwellbin(PT.x,PT.y,PT.l,PT.m,\
                                                 .114,130,130)
    
    #Trace conical mirror
    tracesource(num)
    primaryLL(.01,2.,0.)
    tracefromtest()

    #Subtract all rays outside of WFS
    ind = where(logical_and(abs(PT.y) < 7.3,abs(PT.x) < 7.3))
    PT.vignette(ind=ind)

    #Create conical mirror WFS data
    xang2, yang2, phase2 = reconstruct.southwellbin(PT.x,PT.y,PT.l,PT.m,\
                                                    .114,130,130)

    #Construct phase and angle arrays for wavefront reconstruction
    influence = PT.referencedWavefront(xang,yang,phase,xang2,yang2,phase2)

    return influence

#Compute wavefront of rippled primary with respect
#to perfect cylindrical reference - should see the cone
#Handles cylindrical lens and test optic misalignments
def computerippleinf(num,amp,freq,testX=0.,testY=0.,testZ=0.,\
                    testRX=0.,testRY=0.,testRZ=0.,\
                    cylFieldX=0.,cylFieldY=0.,cylFieldZ=0.,\
                    cylFieldRX=0.,cylFieldRY=0.,cylFieldRZ=0.,\
                    cylNullX=0.,cylNullY=0.,cylNullZ=0.,\
                    cylNullRX=0.,cylNullRY=0.,cylNullRZ=0.):
    #Trace reference wavefront
    tracesource(num,cylNullX=cylNullX,cylNullY=cylNullY,cylNullZ=cylNullZ,\
                cylNullRX=cylNullRX,cylNullRY=cylNullRY,cylNullRZ=cylNullRZ)
    idealcylinder(testX=testX,testY=testY,testZ=testZ,\
                testRX=testRX,testRY=testRY,testRZ=testRZ)
    tracefromtest(cylNullX=cylNullX,cylNullY=cylNullY,cylNullZ=cylNullZ,\
            cylNullRX=cylNullRX,cylNullRY=cylNullRY,cylNullRZ=cylNullRZ,\
            cylFieldX=cylFieldX,cylFieldY=cylFieldY,cylFieldZ=cylFieldZ,\
            cylFieldRX=cylFieldRX,cylFieldRY=cylFieldRY,cylFieldRZ=cylFieldRZ)

    #Subtract all rays outside of WFS
    ind = where(logical_and(abs(PT.y) < 7.3,abs(PT.x) < 7.3))
    PT.vignette(ind=ind)

    #Create reference WFS data
    xang, yang, phase = reconstruct.southwellbin(PT.x,PT.y,PT.l,PT.m,\
                                                 .114,130,130)
    
    #Trace conical mirror
    tracesource(num,cylNullX=cylNullX,cylNullY=cylNullY,cylNullZ=cylNullZ,\
                cylNullRX=cylNullRX,cylNullRY=cylNullRY,cylNullRZ=cylNullRZ)
    rippleprimary(amp,freq,testX=testX,testY=testY,testZ=testZ,\
                testRX=testRX,testRY=testRY,testRZ=testRZ)
    tracefromtest(cylNullX=cylNullX,cylNullY=cylNullY,cylNullZ=cylNullZ,\
            cylNullRX=cylNullRX,cylNullRY=cylNullRY,cylNullRZ=cylNullRZ,\
            cylFieldX=cylFieldX,cylFieldY=cylFieldY,cylFieldZ=cylFieldZ,\
            cylFieldRX=cylFieldRX,cylFieldRY=cylFieldRY,cylFieldRZ=cylFieldRZ)

    #Subtract all rays outside of WFS
    ind = where(logical_and(abs(PT.y) < 7.3,abs(PT.x) < 7.3))
    PT.vignette(ind=ind)

    #Create conical mirror WFS data
    xang2, yang2, phase2 = reconstruct.southwellbin(PT.x,PT.y,PT.l,PT.m,\
                                                    .114,130,130)

    #Construct phase and angle arrays for wavefront reconstruction
    influence = PT.referencedWavefront(xang,yang,phase,xang2,yang2,phase2)

    return influence

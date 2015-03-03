import PyTrace as PT

#Test out Zernike trace
#create list of photons within radius 1 and send in -z direction
#Only put defocus term in, so correct behavior of rays should be obvious

n = 10**4
rho = random.rand(n)*100.
theta = random.rand(n)*2*pi
PT.x = rho*cos(theta)
PT.y = rho*sin(theta)
PT.z = repeat(100.,n)
PT.l = repeat(0.,n)
PT.m = repeat(0.,n)
PT.n = repeat(-1.,n)
PT.ux = copy(l)
PT.uy = copy(l)
PT.uz = copy(l)

pos = [x,y,z]

coeff = array([0.,0.,0.,.1])
rorder,aorder = zernikemod.zmodes(4)

#Put rays in -z position, pointing in positive z axis (toward lens)
transform(0.,0.,0.,pi,0,0)
#Trace to xy plane for interaction with lens
flat()
#Call collimator
lens(2704.01,779.37,10.92,64.01*2,1.64363)
pdb.set_trace()
transform(0,0,-.1,0,0,0) #Gap in achromat
flat()
lens(780.87,-1131.72,15.42,64.01*2,1.51501)
transform(0,0,-1892.019,0,0,0) #Go to focal plane
flat()
pdb.set_trace()

###Call field lens (222 mm focal length)
##lens(205.86,-205.86,5.05,50.,1.51501)
##pdb.set_trace()
###Trace up to focal length
##transform(0.,0.,-195.628,0,0,0)
##pdb.set_trace()
###Trace to focal plane
##flat()
##pdb.set_trace()


#Trace to sphere centered 100 mm above ray starting point
#with RoC 1000 mm
##
##tran.reflect(l,m,n,ux,uy,uz)
##
##tran.transform(x,y,z,l,m,n,ux,uy,uz,0.,0.,250.,0.,0.,0.)
##tran.flat(x,y,z,l,m,n,ux,uy,uz)
##
##pdb.set_trace()

###Trace to 45 degree mirror, reflect, and trace out
##transform(0.,0.,0.,45.*pi/180,0.,0.)
##flat()
##pdb.set_trace()
##reflect()
##pdb.set_trace()
##transform(0.,0.,0.,45*pi/180,0.,0.)
##pdb.set_trace()
##transform(0.,0.,100.,0.,0.,0.)
##pdb.set_trace()
##flat()
##pdb.set_trace()

#Trace to 10 degree surface, refract with N-BK7, verify correct output angle
##tstart = time.time()
##transform(0.,0.,0.,10.*pi/180,0.,0.)
##print time.time()-tstart
##flat()
##tstart = time.time()
##refract(1.,1.51501)

from numpy import *
from matplotlib.pyplot import *

#Test out simulating lens focusing a perfect plane wave
wave = 633.e-6 #633 nm
k = 2*pi/wave

#Define Gaussian beam, make sure phase array is larger than beam
#This prevents aliasing
meshhalfwidth = 40.
N = 2000
x = linspace(-meshhalfwidth,meshhalfwidth,N)
y = copy(x)

xg,yg = meshgrid(x,y)
beamwidth = 2.
##beam = exp(-(xg**2+yg**2)/(2*pi*beamwidth))

#Cylindrical beam
ind = where(sqrt(xg**2+yg**2)<beamwidth)
beam = zeros(shape(xg))
beam[ind] = 1.

#Multiply by quadratic phase factor for focusing lens
f = 500.
beam = beam * exp(-1.j*k/(2*f)*(xg**2+yg**2))

#Propagate to focus, multiply by quadratic phase factor
#and take FFT
#What does frequency mean? There will be scaling to do
#We're at the focus, so we know roughly what the spot size should be
d = 50.
beam = beam * exp(1.j*k/(2*d)*(xg**2+yg**2)) #Basically cancels out the lens
beamfoc = fft.fft2(beam)
freq = fft.fftshift(fft.fftfreq(size(x),x[1]-x[0])*wave*d)

beammag = fft.fftshift(real(beamfoc*conjugate(beamfoc)))

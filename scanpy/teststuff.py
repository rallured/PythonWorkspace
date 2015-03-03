#!/usr/bin/python
from scipy impot *
import curvefitting

#Stuff to generate a test function.

#Generate a gaussian curve and overlay random noise - all in bin form.

#Generate a gaussian in millivolts.
signal = normal(loc = 1500, scale = 50, size = (500)).astype('int')

#Generate Poisson noise
noise = zeros(256)

#Model every bin as its own Poisson process.
for i in range(len(noise)):
    noise[i] += poisson(lam = 12)


sig_hist = histogram(signal, bins = 256, range = (0, 5000))[0]

#Total signal + poisson bin noise.
total = noise + sig_hist

#Range of x values - just the bin #'s
x = arange(256).astype('int')

#Fit to a gaussian with a constant offset (this is the nterms = 4 part)
fit_out = curvefitting.gaussguess(x, total, nterms = 4)

#Plot actual data
plot(x, total)

#Overlay fit
plot(x, curvefitting.gaussian((fit_out[0], fit_out[1], fit_out[2][0], fit_out[3]), x))

#These lines added for SVN commit test.

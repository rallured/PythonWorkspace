from numpy import *
import matplotlib.pyplot as plt
import pyfits
import os
import pdb
import sys
import scipy.interpolate as inter
import matplotlib as mpl
import scipy.ndimage.interpolation as interpol

from ProgressBar import ProgressBar

############################################################################################
# Initialize
##############################################
# For the progress bars.
custom_options = {
    'end':100,
    'width':100
}

scan_length = 1000
N = 10**2
############################################################################################  
# Setting up the correct test distributions.
##############################################
def make_distribution(location,sigma,N_counts):
    test_dist = random.normal(loc = location, scale = sigma, size = random.poisson(lam = N_counts))
    plt.ioff()
    your_mom = plt.hist(test_dist,bins=50)
    counts = your_mom[0]
    x = asarray([mean(your_mom[1][k:k+2]) for k in range(len(counts))])    
    ind = where(counts > 0)
    x_scrub,counts_scrub = x[ind],counts[ind]
    return x_scrub,counts_scrub

def make_all_counts_arrays(x_scrub,LSF_scrub):
    LSF_all_counts = zeros(sum(LSF_scrub)) + 1
    for i in range(len(x_scrub)-1):
	counts = LSF_scrub[i]
	bin_size = 1.0/counts
	if i == 0:
	    x_all_counts = i + arange(0,1,bin_size) + 0.5*bin_size
	else:
	    x_all_counts= hstack((x_all_counts,i+arange(0,1,bin_size) + 0.5*bin_size))
    return x_all_counts,LSF_all_counts

def measure_HEW_paper_method(x,counts):
    pdb.set_trace()
    N = sum(counts)
    x_all_counts,LSF_all_counts = make_all_counts_arrays(x,counts)
    cdf = cumsum(LSF_all_counts)
    percentage = cdf/N
    sigma_percentage = percentage*sqrt(1.0/cdf + 1.0/N)
    
    hew = x_all_counts[abs(percentage - 0.75).argmin()] - x_all_counts[abs(percentage - 0.25).argmin()]
    hew_max = x_all_counts[abs(percentage - sigma_percentage - 0.75).argmin()] - x_all_counts[abs(percentage + sigma_percentage - 0.25).argmin()]
    hew_min =  x_all_counts[abs(percentage + sigma_percentage - 0.75).argmin()] - x_all_counts[abs(percentage - sigma_percentage - 0.25).argmin()]
    return hew,hew_max,hew_min,x_all_counts,cdf

# Measure once for comparison.
x,counts = make_distribution(1000.0,5.0,N)
nom_HEW,nom_HEW_max,nom_HEW_min,nom_x,nom_cdf = measure_HEW_paper_method(x,counts)

# Set up Monte Carlo method.
scan_hews = zeros(scan_length)
probar = ProgressBar(**custom_options)
progress_array = asarray(range(0,scan_length,scan_length/100))

for j in range(scan_length):
    if any(j == progress_array):
	print '\b'*250 + 'Statistics Simulation Progress: ', probar+1,
	sys.stdout.flush()
    
    sim_x,sim_counts = make_distribution(1000,5.0,N)
    scan_hews[j] = measure_HEW_paper_method(sim_x,sim_counts)[0]

# Bootstrap method HEW error. Paper method HEW error is nom_HEW_max - nom_HEW or nom_HEW - nom_HEW_min (roughly, but not quite, symmetric).
MC_bootstrap_hew_error = std(scan_hews)



#Quick function defs to take and plot data on the windows lappy using Daly's 
#MCA

import daqserial
from scipy import *
from matplotlib import pyplot

root = 'C:\\Transferred Users\\ZPrieskorn\\My Documents\\SeanMattingly\\Python\\AutoScan\\HandData\\'

def get_hist(ac_time = 120):
    a,b = daqserial.acquire(120)
    absc = arange(256)
    ord_a = histogram(a, bins = 256, range = (0, 255))
    ord_b = histogram(b, bins = 256, range = (0, 255))
    return (absc, ord_a, ord_b)

def get_save_hist(a_fname, b_fname, ac_time = 120):
	a,b = daqserial.acquire(ac_time)
	absc = arange(256)
	ord_a = histogram(a, bins = 256, range = (0, 255))[0]
	ord_b = histogram(b, bins = 256, range = (0, 255))[0]
	savetxt(root + a_fname, ord_a)
	savetxt(root + b_fname, ord_b)
	return (absc, ord_a, ord_b)


#balh balh blah for SVN commit test
#Second blahblah for SVN Commit test from Windows.

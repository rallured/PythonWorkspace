import numpy as np
import pdb

def createDist(N):
    """Return an array of positions with number of counts
    varying via Poisson statistics"""
    return np.random.normal(loc=1000.,scale=5.,size=np.random.poisson(lam=N))

def computeHEWsCasey(N):
    """Formulate CDF """
    #Create distribution
    x = createDist(N)
    #Form CDF
    N = np.size(x)
    cdf = np.cumsum(np.repeat(1./N,N))
    #Reformulate x positions into absolute deviations from the mean
    x = np.abs(x - np.mean(x))
    x.sort()
    #At this point you can plot(x,cdf) to see that it is indeed the CDF
    #This isn't necessary for calculating the HEWs, but it
    #illustrates the idea.
    #plot(x,cdf)
    
    #Pick off HEW by choosing the median position...half the counts are
    #within this radius, half are out
    hew = np.median(x)*2. #Could have used x[N/2]
    #Now choose upper and lower bounds by choosing the position where
    #(N+sqrt(N))/2 and (N-sqrt(N))/2 counts are within radius
    hewmax = x[(N+np.sqrt(N))/2]*2.
    hewmin = x[(N-np.sqrt(N))/2]*2.
    #We now have the hew and the bounds from Casey's method
    return hew,hewmax-hew,hew-hewmin

def bruteForce_MC(N,M):
    """Brute force HEW calculation Monte Carlo. Simulation runs
    M times. Assume counts vary by Poisson statistics with mean N
    counts. For each run, create a distribution and return its
    HEW (median deviation from mean). Then compute the standard
    deviation of the HEWs. This standard deviation should be
    equal to Casey's hewmax-hew or hew-hewmin
    """
    hewlist = np.zeros(M)
    for i in range(M):
        x = createDist(N)
        x = np.abs(x-np.mean(x))
        x.sort()
        hewlist[i] = np.median(x)*2.
    return np.mean(hewlist), np.std(hewlist)

#The bounds returned by computeHEWsCasey should be approximately
#equal to the bound returned by the brute force MC, but they are not
print computeHEWsCasey(100000)
print bruteForce_MC(100000,1000)

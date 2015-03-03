import testflex95
import PytranTrace
import time
from numpy import *

tstart =time.time()
inp = random.randn(1000000)
print testflex95.singleloop(inp)
##PytranTrace.huygenstest(.2,arange(1000.),arange(1000.))
print time.time()-tstart

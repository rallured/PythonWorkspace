import userKov.pyGeneralRoutines.generalRoutines
import dislin
import os
from math import *
import pylab
from numpy import *
import Ccoating

def plotta():
	dislin.disini()
	dislin.scrmod ('revers')
	dislin.pagera ()
	dislin.color('white')
	dislin.graf   (1., 20., 5., 5., 0., 1., 0., 0.1)
	n=len(xTeo)-1
	dislin.color  ('red')
	#dislin.plot (xTeo,yTeo)
	dislin.curve  (xTeo, yTeo,n)
	dislin.color  ('green')
	dislin.curve  (xTeo, yFresnel,n)
	dislin.disfin()

def plotta2(x,y):
    dislin.plot(x,y)
    dislin.disfin()

'''
if __name__=="old__main__":
	refTeoFile="IrImd.dat"
	xTeo=load(refTeoFile)[:,0]
	xTeo=[x/1000 for x in xTeo]
	yTeo=load(refTeoFile)[:,1]
	#print xTeo,yTeo
	index=loadInd()
	yFresnel=fresnel(1,index,xTeo)
	k=open("kovIr.dat","w")
	#dislin.plot(xTeo,yFresnel)
	#dislin.plot(xTeo,yTeo)
	dislin.disfin()
	#mette output su file
	for i,x in enumerate(xTeo):
		k.write("%s\t%s\n"%(x*1000,yFresnel[i]))
	k.close()
'''

if __name__=="__main__":
	refTeoFile="IrEner.txt"
	xTeo=pylab.load(refTeoFile)[:,0]
	#xTeo=[x/1000 for x in xTeo]
	#xTeo=range(90)
	yTeo=pylab.load(refTeoFile)[:,1]
	iridioTot=Ccoating.Ccoating() #inizializza con valore di default (iridio)
	
	yFresnel=iridioTot.reflex(xTeo,0.0030)
	k=open("kovIr.dat","w")
	#dislin.plot(xTeo,yFresnel)
	#dislin.plot(xTeo,yTeo)
	#dislin.disfin()
	for i,x in enumerate(xTeo):
		k.write("%s\t%s\n"%(x,yFresnel[i]))
	k.close()







    
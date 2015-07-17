import os
from math import *
from numpy import *
import pylab 
from numpy import interp


def snellCos(n1,n2,angle):
    '''restituisce il coseno dell'angolo di rifrazione, dati indici di rif e angolo
    di incidenza'''
    return sqrt(n2**2-sin(angle)**2)[0]/n2

def critAng(n1,n2):
    '''restituisce l'angolo critico calcolato dalle parti reali,
    dati gli indici di rifrazione, se negativo, dalla parte di n2'''
    if n1>n2:
        ca=asin(n2/n1)
    else:
        ca=-asin(n1/n2)
    return ca



class Ccoating(object):
    '''e' un coating (per ora solo monostrato)
    properties: refIndex,coatingFile,ener;
    methods:reflex(self,angle,ener=0)'''

    def readFile(coatingFile):
        ''' read a file with format:
        Thickness(A)	Material	Roughness(A)
        and set the related properties.'''
        logging.debug('reading coating file %s...'%coatingFile)
        try:
            coating=asciitable.read(coatingFile)
            self.dspacing=float(coating["Thickness(A)"])
            self.roughness=float(coating["Roughness(A)"])
            self.materials=coating["Material"]
        except:
            logging.critical('Problem in reading file %s'%coatingFile)
        
    def __init__(self,coatingFile="",ener=0,angle=0):
        #self.refIndex=(0.99996558,-2.3678801e-06)
        if ener: self.ener=ener
        if angle: self.angle=angle
        self.file=coatingFile
        if coatingFile: readFile(coatingFile)
        self.__changed=True

    def __getInd(self):
        if self.__changed :     #ricarica solo se cambiato file o energie
            n=pylab.load(self.matFile)
            n[:,1]=1-n[:,1]
            n[:,2]=-n[:,2]
            n1=interp(n[:,1],n[:,0],self.__ener)
            n2=interp(n[:,2],n[:,0],self.__ener)
            self.__index=[(nn1,nn2) for nn1,nn2 in zip(n1,n2)]
            self.__changed=False
        return self.__index
    def __setInd(self, value):
        	 self.__index= value
    def __delInd(self): pass
    refIndex = property(__getInd, __setInd, __delInd, "interpolated refraction index")

    def __getMatFile(self): return self.__mf
    def __setMatFile(self, value):
        self.__mf = value
        self.__changed=True
    def __delMatFile(self): pass
    matFile = property(__getMatFile, __setMatFile, __delMatFile, "material refraction index file") 

    def __getMatDict(self): return self.__mf
    def __setMatDict(self, value):
        self.__mf = value
        self.__changed=True
    def __delMatDict(self): pass
    matDict = property(__getMatDict, __setMatDict, __delMatDict, "material refraction index file") 

    def __getEn(self): return self.__ener
    def __setEn(self, value):
        self.__ener = array(value)
        self.__changed=True
    def __delEn(self): pass
    ener = property(__getEn, __setEn, __delEn, "energie(s) list for calculation.") 

    def __getAn(self): return self.__angle
    def __setAn(self, value):
        self.__ener = array(value)
        self.__changed=True
    def __delAn(self): pass
    ener = property(__getAn, __setAn, __delAn, "Angle(s) list for reflex calculation.") 
    
    def fresnel(self,angle):
        ''' restituisce i valori di riflettivita' per i valori di indice <n> e
        lista angoli incidenza <angle>, indice scalare complesso
        (tupla o lista di tuple) e angolo come lista o singolo valore'''
        '''proprieta': indici, reflex'''
  
        ref=[]
        n1=complex(1,0)
        n=self.refIndex
        if ener==0:
            logging.error('Energy not defined.')
            return 0
        if  iterable(angle) and iterable(n):    #tutte e due liste
            #da implementare restituzione di matrice
            logging.warning('Iterable angle AND index, option not implemented.')
            pass
        elif iterable(angle):   #angoli lista, energia scalare
            logging.warning('Iterable angle, option not implemented.')
            pass
        elif iterable(n):
            #angle is a scalar
            complex*16 re2(pmaxEn), rs2(pmaxEn), ro2(pMaxEn)
            ro2(1:nener)=(1-delox(1:nener)-ci*betox(1:nener))**2 !(refraction index "odd")**2
            re2(1:nener)=(1-delex(1:nener)-ci*betex(1:nener))**2 !(refraction index "even")**2
            rs2(1:nener)=(1-delsx(1:nener)-ci*betsx(1:nener))**2 !(refraction index substrate)**2
            reflexCore(dSpacing,n_layers,reflex,nener,angle,self.roughness)
            return ref
        else:
            #single angle and single energy. Return a scalar
            logging.warning('scalar angle and energy, option not implemented.')
            pass

    def reflex(self,angle=0,ener=0):
        self.ener=ener
        self.angle=angle    
        return self.fresnel(angle) #se cambiato ricalcola gli indici per self.__changed
        
    

if __name__=="__main__":
	refTeoFile="IrEner.txt"
	xTeo=pylab.load(refTeoFile)[:,0]
	yTeo=pylab.load(refTeoFile)[:,1]
	iridioTot=Ccoating() #inizializza con valore di default (iridio)
	
	yFresnel=iridioTot.reflex(xTeo,0.0030)
	k=open("kovIr.dat","w")
	for i,x in enumerate(xTeo):
		k.write("%s\t%s\n"%(x,yFresnel[i]))
	k.close()        

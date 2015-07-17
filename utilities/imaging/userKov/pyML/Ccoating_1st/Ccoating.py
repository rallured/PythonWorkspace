import os
from numpy import *
import pylab 
from numpy import interp

def Cexp(n):        #n=a+ib
    l=exp(n.real)   #fattore moltiplicativo e**a
    i=n.imag        #b
    ei=complex(cos(i),sin(i))            #e**ib
    return l*ei

def Csin(ang):
    j=complex(0,1)
    return (Cexp(complex(ang)*j)-Cexp(-complex(ang)*j))/(2*j) 

def Ccos(ang):
    j=complex(0,1)
    return (Cexp(complex(ang)*j)+Cexp(-complex(ang)*j))/2

def Cacos(val):
    alfa=val.real
    beta=val.imag
    gamma=sqrt(alfa**2+beta**2)
    b=log(gamma+sqrt(gamma**2+1))
    a=acos(alfa/sinh(b))
    return complex(a,b)

def Ctan(ang):
    return Csin(ang)/Ccos(ang)

def pro(a):
    return ["%s: la radice e' %s, infatti al quadrato da'-> %s"% (a,i,i**2)
            for i in Csqrt(a)]

def Csqrt(n):
    beta=n.imag
    alfa=n.real
    if (beta==0 and alfa >= 0):return complex (sqrt(alfa)),-complex (sqrt(alfa))
    b1=sqrt((-alfa+(sqrt(alfa**2+beta**2)))/2)
    #b2=sqrt((-alfa-(sqrt(alfa**2+beta**2)))/2)
    a1=beta/(b1*2)
    #a2=beta/(b2*2)
    return complex(a1,b1),complex(-a1,-b1) #,complex(a2,b2)

def snellCos(n1,n2,angle):
    '''restituisce il coseno dell'angolo di rifrazione, dati indici di rif e angolo
    di incidenza'''
    return Csqrt(n2**2-Csin(angle)**2)[0]/n2

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
    properties: refIndex,matFile,ener;
    methods:reflex(self,angle,ener=0)'''
    def __init__(self,matFile="iridio.dat",ener=-1):
        #self.refIndex=(0.99996558,-2.3678801e-06)
        self.ener=ener
        self.matFile=matFile
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

    def __getEn(self): return self.__ener
    def __setEn(self, value):
        self.__ener = array(value)
        self.__changed=True
    def __delEn(self): pass
    ener = property(__getEn, __setEn, __delEn, "energie(s) list for calculation") 
    
    def fresnel(self,angle):
        ''' restituisce i valori di riflettivita' per i valori di indice <n> e
        lista angoli incidenza <angle>, indice scalare complesso
        (tupla o lista di tuple) e angolo come lista o singolo valore'''
        '''proprieta': indici, reflex'''
        ref=[]
        n1=complex(1,0)
        n=self.refIndex
        if  iterable(angle) and iterable(n):    #tutte e due liste
            #da implementare restituzione di matrice
            
            pass
        elif iterable(angle):   #angoli lista, energia scalare
            for ang1 in angle:
                rad=Csqrt(n2**2-Ccos(ang1)**2)[1]
                c=Csin(ang1)
                refPar=(c-rad)/(c+rad)
                refPer=(c*n2**2-rad)/(c*n2**2+rad)
                ref.append((abs(refPar)**2+abs(refPer)**2)/2)
            return ref
        elif iterable(n):
            ang1=angle
            for i,j in n:
                n2=complex(i,j)
                rad=Csqrt(n2**2-Ccos(ang1)**2)[1]
                c=Csin(ang1)
                refPar=(c-rad)/(c+rad)
                refPer=(c*n2**2-rad)/(c*n2**2+rad)
                ref.append((abs(refPar)**2+abs(refPer)**2)/2)
            return ref
        else:
            rad=Csqrt(n2**2-Ccos(ang1)**2)[1]
            c=Csin(ang1)
            refPar=(c-rad)/(c+rad)
            refPer=(c*n2**2-rad)/(c*n2**2+rad)
            ref=(abs(refPar)**2+abs(refPer)**2)/2
            return ref

    def reflex(self,angle,ener=0):
        if ener:
            self.ener=ener
            self.__changed=True
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

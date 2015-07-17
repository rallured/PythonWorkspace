import os
from pylab import *

'''coating a doppio strato, fatto grezzamente raddoppiando tutto.'''
'''
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
'''
def snellCos(n1,n2,angle):
    '''restituisce il coseno dell'angolo di rifrazione, dati indici di rif e angolo
    di incidenza'''
    return sqrt(n2**2-sin(angle)**2)[0]/n2

def critAng(n1,n2):
    '''restituisce l'angolo critico, dati gli indici di rifrazione, se negativo, dalla parte di n2'''
    if n1>n2:
        ca=asin(n2/n1)
    else:
        ca=-asin(n1/n2)
    return ca

class layer(object):
    '''refraction index for a list of energies and a given material'''
    def __init__(self,mat="iridio.dat",ener=-1):
        #self.refIndex=(0.99996558,-2.3678801e-06)
        self.ener=ener
        self.mat1=mat
        self.__changed=True
        
    def __getInd1(self):
        if self.__changed :     #ricarica solo se cambiato file o energie
            n=load(self.mat)
            n[:,1]=1-n[:,1]
            n[:,2]=-n[:,2]
            n1=interp(self.__ener,n[:,1],n[:,0])
            n2=interp(self.__ener,n[:,2],n[:,0])
            self.__index=[(nn1,nn2) for nn1,nn2 in zip(n1,n2)]
            self.__changed=False
        return self.__index
    def __setInd1(self, value):
        	 self.__index= value
        	 self.__changed=True
    def __delInd1(self): pass
    refIndex1 = property(__getInd, __setInd, __delInd, "interpolated refraction index")

    def __getMat(self): return self.__mf
    def __setMat(self, value):
        self.__mf = value
        self.__changed=True
    def __delMat(self): pass
    matFile = property(__getMat, __setMat, __delMat, "%s refraction index file"%(self.__mf))
    

class Ccoating(object):
    '''e' un coating bistrato '''
    def __init__(self,mat1="iridio.dat",mat2="carbon.dat",ener=-1):
        self.ener=ener
        self.layer1=layer(mat1,ener)
        self.layer2=layer(mat2,ener)
        self.__changed=True

    def __getMat1(self): return self.layer1.mat1
    def __setMat1(self, value):
        self.layer1.mat = value
        self.layer1.changed=True
    def __delMat1(self): pass
    matFile = property(__getMat1, __setMat1, __delMat1, "layer1 %s refraction index file"%(self.layer1.mat1))

    def __getMat2(self): return self.layer1.mat1
    def __setMat2(self, value):
        self.layer2.mat = value
        self.layer2.changed=True
    def __delMat2(self): pass
    matFile = property(__getMat2, __setMat2, __delMat2, "layer2 %s refraction index file"%(self.layer1.mat1)) 

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
        n=self.refIndex1
        nn=self.refIndex2
        ang1=angle
        for c,i,j in zip(n):
            n2=complex(i,j)
            rad=sqrt(n2**2-cos(ang1)**2)[1]
            c=sin(ang1)
            refPar=(c-rad)/(c+rad)
            refPer=(c*n2**2-rad)/(c*n2**2+rad)
            ref.append((abs(refPar)**2+abs(refPer)**2)/2)
        return ref

    def reflex(self,angle,ener=0):
        if ener:
            self.ener=ener
            self.__changed=True
        return self.fresnel(angle) #se cambiato ricalcola gli indici per self.__changed
        
def test():
    a=Ccoating('Ir.dat','a-C.dat',arange(15))
    print
    a.reflex(0.001)

if __name__=="__main__":
    test()
                
        

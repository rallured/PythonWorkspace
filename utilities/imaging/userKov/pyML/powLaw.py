import os
from numpy import *
    

def pl4(pars,n_bil):
    a,b,c,gamma=pars
    t=[]
    for i in range (1,n_bil+1):
        d=a/(b+i)**c
        t.append([d*gamma,d*(1-gamma)])
    return array(t)

def fromd1dnc((d1,dN,c),n_bil):
    '''return a,b,c parameters of pl4, as calculated from first and last tichness
    and c parameter'''
    b=(n_bil-1)/((d1/dN)**(1/c))-1
    a=dN*(b+n_bil)**c
    pars=a,b,c
    return pars,n_bil


def tFromd1dn((d1,dN,c),n_bil,gamma):
    '''return thickness from pl4 applied to the parameters extracted from
    d1,dn,c '''
    a,b,c=fromd1dnc((d1,dN,c),n_bil)[0]
    return pl4((a,b,c,gamma), n_bil)


def dist(x,y):
    return(x-y)**2
    
def findc(target,toll,c=[0.1,4.],law=powLaw):
    '''find c in crange=[cmin,cmax], using bisection method with toll
    tollerance in tichkness units, for the target value [index,value]  '''
    d1,dN,n_bil,gamma=271.,33.,290,1.
    tgtInd,tgtVal=target
    def sp(c): return tFromd1dn((d1,dN,c),n_bil,gamma)[:,0]
    cmin,cmax=c
    newc=(cmin+cmax)/2

    if dist(tgtVal,sp(cmin)[tgtInd])>dist(tgtVal,sp(cmax)[tgtInd]): #pone quello che sara' il limite del nuovo intervallo
        newLim=cmax
        trash=cmin
    else:
        newLim=cmin
        trash=cmax
    print tgtVal,newLim,sp(newLim)[tgtInd],c,dist(tgtVal,sp(cmin)[tgtInd]),\
          dist(tgtVal,sp(newc)[tgtInd]),dist(tgtVal,sp(cmax)[tgtInd])    
    if dist(tgtVal,sp(newc)[tgtInd])> dist(tgtVal,sp(trash)[tgtInd]):
        print "fuori dai margini: ",c
        return None
    else:
        c=[newc,newLim]
        c.sort()
        if dist(tgtVal,sp(newc)[tgtInd])<toll:
            cres=fromd1dnc((d1,dN,newc),n_bil)[0][2]
        else:
            cres=findc(target,toll,c)
    return cres

if __name__=="__main__":
    c=findc([144.,38],10e-7)
    print c   

# integrare nel seguito facendo classe...
    
testSi=["indici","del","Si"]
testW=["indici","del","W"]

# schema proprieta' get e set ##############
class C(object):
    def __init__(self): self.__x = None
    def getx(self): return self.__x
    def setx(self, value): self.__x = value
    def delx(self): del self.__x
    x = property(getx, setx, delx, "I'm the 'x' property.")
# ###########################################    
    
class multilayer(object):
    def reflexCalc(self,Angle,EnRange=(0,100),npoints=200):
        self.__reflex=array([[i,1] for i in range (1,100)])
        
    
    def __init__(self,n_bil,par,testInd=testW,law=pl4):
        self.n_bil=n_bil
        self.__par=par
        self.testInd=testInd
        self.law=law
        self.d=self.law(self.__par,self.n_bil)        

    def __setpar(self,val):
        self.__par=val
        self.d=self.law(val,self.n_bil)   
    def __getpar(self):return self.__par

    def __getref(self):
        self.reflexCalc()
        return self.__reflex
    
    par = property(__getpar, __setpar, "Multilayer parameters")
    reflex = property(__getref, "Multilayer reflectivity")    

    
    

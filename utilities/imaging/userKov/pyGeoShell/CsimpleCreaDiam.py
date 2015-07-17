import namelist_class
import os
import math
from userKov.pyML import Ccoating
from pylab import *

class shellGeo:
    '''classe che regola il dimensionamento delle shell, a partire
    dai parametri. mettere eccezioni per file mancanti'''
    XMM=(0.0016,-0.02)
    SX=(0.0006,-0.04)


        
    def aeff(self,ener,outFile=None):
        '''restituisci l'area efficace totale sommata sulle shell
        alle energie della lista ener, se si passa un nome di file scrive
        i valori'''
        aTot=zeros(len(ener))
        i=0
        for ang,shellD in zip(self.angle,self.d):
            i=i+1
            print i
            dR=math.tan(ang)*self.SH*10
            acoll=math.pi*dR*shellD/100           #math.pi*dR*(shellD+dR)/100 c'e fattore 2 per diam?
            print shellD,acoll,ang
            aeff=acoll*array(self.coating.reflex(ang,ener))**2
            open("shell_%03i.txt"%i,"w").write("\n".join(["%s\t%s" %(e,a) for e,a in zip(ener, aeff)]))
            aTot=aTot+aeff
        #if outFile: self.aWrite(outFile)
        return aTot


    def diam(self):
        '''tutte le misure in cm,
        risultati convertiti in mm all'uscita, restituisce lista di diametri e
        lista di spessori'''
        print "stai usando il metodo diam in via di estinzione, modifica il"
        print "programma chiamante per usare proprieta' diam e spes con metodo"
        print "calcGeo"
        ind= range(self.nshell)
        m = self.thickVal[0]
        q = self.thickVal[1]
        ds=[]
        spes=[]
        d=self.Dmax
        for i in ind:
            spes.append(m*d*10+q) #la lista degli spessori e' in mm
            #print (d*10*m + q),spes[-1]
            d=d-2*(spes[-1]/10)
            ang= math.atan2(d,(2*self.FL*100))/4
            dmed=d-2*self.SH*math.tan(ang)
            ds.append(dmed*10)
            d=dmed-2*self.SH*math.sin(math.radians(self.ff))
        return ds, spes

    def diam2 (self):
        #tutte le misure in mm,imita fortran
        self.angle=[]
        ind= range(self.nshell)
        m = self.thickVal[0]
        q = self.thickVal[1]
        ds=[]
        spes=[]
        r=self.Dmax/2
        SHmm=self.SH*10
        FLmm=self.FL*1000
        alkw=SHmm*math.sin(math.radians(self.ff))
        for i in ind:
            s=2*r*m + q
            if s<self.thickVal[2]:s=self.thickVal[2]
            spes.append(s)
            r=r-(spes[-1])
            th= math.atan2(r,FLmm)/4
            self.angle.append(th)
            a=2*math.tan(th)*SHmm
            r0=(math.sqrt(a**2+4*r**2)-a)/2
            while True:
                th= math.atan2(r0,FLmm)/4
                a=2*math.tan(th)*SHmm
                r1=(math.sqrt(a**2+4*r**2)-a)/2
                #print i,abs(r1-r0)
                if (abs(r1-r0)<=0.0001):break
                r0=r1
            rmed=r1
            ds.append(rmed*2)
            r=rmed-alkw
        return ds, spes

    def __weights(self):
        diam=self.d
        spes=self.spes
        ang=self.angle
        W=[]
        for i in range(len(diam)):
            Dmed=diam[i]
            Dmax=Dmed+2*math.tan(ang[i])*self.SH*10
            Dmin=Dmed-2*math.tan(3*ang[i])*self.SH*10
            p=((Dmax+Dmin)/2+Dmed)*math.pi*spes[i]*self.SH*self.wd/100000
            #print Dmax,Dmed,Dmin,spes[i]
            W.append(p)
        self.weights=W        

    def calcGeo(self):
        '''calcola i parametri geometrici della configurazione:
        diametri, spessori, angoli, pesi '''

        #tutte le misure in mm,imita fortran, da diam2
        self.angle=[]
        ind= range(self.nshell)
        m = self.thickVal[0]
        q = self.thickVal[1]
        ds=[]
        spes=[]
        r=self.Dmax/2
        SHmm=self.SH*10
        FLmm=self.FL*1000
        alkw=SHmm*math.sin(math.radians(self.ff))
        for i in ind:
            s=2*r*m + q
            if s<self.thickVal[2]:s=self.thickVal[2]
            spes.append(s)
            r=r-(spes[-1])
            th= math.atan2(r,FLmm)/4
            #self.angle.append(th)
            a=2*math.tan(th)*SHmm
            r0=(math.sqrt(a**2+4*r**2)-a)/2
            while True:
                th= math.atan2(r0,FLmm)/4
                a=2*math.tan(th)*SHmm
                r1=(math.sqrt(a**2+4*r**2)-a)/2
                if (abs(r1-r0)<=0.0001):break
                r0=r1
            rmed=r1
            th= math.atan2(rmed,FLmm)/4
            self.angle.append(th)
            ds.append(rmed*2)
            r=rmed-alkw
        self.d,self.spes= ds, spes
        self.__weights()

    def a1Distr (self,n,oaArcMin):
        '''return as a function the distribution of incidence angles on the first surface
        of the n-th shell'''
        alfa=self.angle[n]
        theta=(math.pi/180)*oaArcMin/60
        f=cos(alfa)*sin(theta)
        g=cos(theta)*sin(alfa)
        def fun (x): return arccos(f*cos(x)+g)
        print "test: passata def funzinoe f(x)=acos(f*cos(x)+g) con f,g = %s%s"%(f,g)
        print "funz(0)=",fun(0)
        return fun
        
        

    def __init__(self,focal,shellHeight,Dmax,fillingFactor,nshell,wallDensity,thickVal=XMM):
        self.__geoChanged=1
        self.FL=focal
        self.SH=shellHeight
        self.Dmax=Dmax
        self.ff=fillingFactor
        self.focal=focal
        self.nshell=nshell
        self.wd=wallDensity
        self.thickVal=thickVal
        self.coating=Ccoating.Ccoating("iridio.dat")
        self.calcGeo()

    def mean_ff(self):
        '''restituisce il filling factor (inteso come rapporto tra area proiettata e
        area occupata) reale e quello calcolato con
        formula approssimata, ottenuta dagli integrali'''
        ff=self.ff
        F=self.FL
        H=self.SH*10
        Rmax=self.Dmax/2
        Rmin=self.d[-1]/2
        aGeo=math.pi*(Rmax**2-Rmin**2)
        thF4=4*math.tan(ff)*F
        print '----------'

        #calcola area raccolta totale 
        Ash=0
        rIn=[]
        for shd,ang in zip(self.d, self.angle):
            rIn.append(shd/2+H*math.tan(ang))   #crea lista dei diametri all'ingresso
            Ash=Ash+math.pi*((rIn[-1])**2-(shd/2)**2)    
        ffRatio=Ash/aGeo
        print 'da calcolo: ', ffRatio

        #calcola con sommatoria
        aShSum=0
        rIn.append(self.d[-1]/2-math.tan(ff)*H)
        for i in range(self.nshell):
            R1=rIn[i+1]/2
            R2=rIn[i]/2
            aShSum=aShSum+math.pi*(R2**2-R1**2)*R1/(R1+thF4)
            #print, apar
        print 'da sommatoria: ', aShSum/aGeo

        #ora lo calcola con la formula
        teo=2*math.pi*((Rmax-Rmin)+(thF4**2)*math.log((thF4+Rmax)/(thF4+Rmin)))
        print 'da formula: ', teo/aGeo
        
class CSimpleCreaDiam(namelist_class.Namelist,shellGeo):
    '''classe per creazione di diametri con parametri letti da namelist,
    manca qualche perfezionamento per gestire tutti i casi della namelist,
    come i vari flag o il diametro max non fornito direttamente, accetta sia
    namelist di oa che di oa2'''
    
    def __init__(self,nl):
        namelist_class.Namelist.__init__(self,nl)
        self.dirorigin=os.path.join(self["WORKDIR"].strip("\"' "),self["DIRORIGIN"].strip("\"' "))
        self.dirrisult=os.path.join(self["WORKDIR"].strip("\"' "),self["DIRRISULT"].strip("\"' "))
        try:   #prova a usare i parametri per OA2
            shellGeo.__init__(self,self["F_LENGTHdaImp_m"],
                              self["F_HEIGHTdaImp_cm"],self["massimo"],self["FieldOfViewDeg"],
                              self["nshDaImp"],self["WallDensity"],(self["ThickM"],self["ThickQ"],self["minThickmm"]))
        except: #senno' inizializza come OA
            q=self["ThickExtSh"]/self["massimo"]
            shellGeo.__init__(self,self["F_LENGTHdaImp_m"],
                          self["F_HEIGHTdaImp_cm"],self["massimo"],self["FieldOfViewDeg"],
                          self["nshDaImp"],self["WallDensity"],(0,q,0))

   # def __getattr__(self,name):
    #    print "sto recuperando il valore di %s"%(name)
     #   return shellGeo.__getattr__(name)


def plottaConfronto(nl,comblist):
    '''presa da versione di test, ma non testata in questo contesto. crea file di plot
    per confronto con vecchi risultati'''
    p=open("confrDiam.plt","w")
    for comb in comblist:
        ndf="F"+str(comb[0]).split(".")[0]+"_ff"+str(comb[1]).split(".")[1].zfill(3)+"_D"+str(comb[2])+".txt"
        nl2.FL,nl2.ff,nl2.Dmax=comb
        di,sp=nl2.diam2()
        #per plottare confronto con vecchi diametri
        p.write("plot '"+ndf+"' u 0:1 w lp")
        fcnf1="mlff"+str(comb[1]).split(".")[1].zfill(3)+str(comb[0]).split(".")[0]+"mmax"+str(comb[2])+".txt"
        fcnf2="m2ff"+str(comb[1]).split(".")[1].zfill(3)+str(comb[0]).split(".")[0]+"mmax"+str(comb[2])+".txt"
        if os.path.exists(os.path.join("af_files",fcnf1)):
           p.write(",\\\n'"+os.path.join("af_files",fcnf1)+"' u 0:1 w lp")
        if os.path.exists(os.path.join("af_files",fcnf2)):
           p.write(",\\\n'"+os.path.join("af_files",fcnf2)+"' u 0:1 w lp")
        '''if os.path.exists(os.path.join("frlike",ndf)):
           p.write(",\\\n'"+os.path.join("frlike",ndf)+"' u 0:1 w lp")
        fcnf3="ml_ff"+str(comb[1]).split(".")[1].zfill(3)+"_"+str(comb[0]).split(".")[0]+"m_max"+str(comb[2])+".txt"
        if os.path.exists(os.path.join("extractedDiam",fcnf3)):
           p.write(",\\\n'"+os.path.join("extractedDiam",fcnf3)+"' u 0:1 w lp")
        fcnf4="ir_ff"+str(comb[1]).split(".")[1].zfill(3)+"_"+str(comb[0]).split(".")[0]+"m_max"+str(comb[2])+".txt"
        if os.path.exists(os.path.join("extractedDiam",fcnf4)):
           p.write(",\\\n'"+os.path.join("extractedDiam",fcnf4)+"' u 0:1 w lp")'''
        p.write("\n")
        p.write("pause -1\n\n")

        p.write("set term png\n")
        p.write("set out '"+ndf+".png'\n")
        p.write("replot\n")
        p.write("set term win\n set out\n\n")
    p.close()                  

if "__name__"=="x__main__":
    nl2=CSimpleCreaDiam("imp_OffAxis.txt")

    focal=[i*2.5 for i in range(8,12)]
    ff=[0.00,0.07,0.15]
    D=[60,70]
    comblist=[[f1,f2,f3] for f1 in focal for f2 in ff for f3 in D]

    #comblist=[[20.0,0.0,60]]
    for comb in comblist:
        ndf="F"+str(comb[0]).split(".")[0]+"_ff"+str(comb[1]).split(".")[1].zfill(3)+"_D"+str(comb[2])+".txt"
        nl2.FL,nl2.ff,nl2.Dmax=comb
        di,sp=nl2.diam()
        f= open(ndf,"w")
        for d in di[-1:0:-1]:
            f.write(str(d)+"\n")
        print ndf,nl2.FL,nl2.SH,nl2.Dmax,nl2.ff,nl2.nshell,nl2.wd,(nl2.thickVal[0],nl2.thickVal[1])
        print "\n"
        f.close()
    
    #plottaConfronto(nl2,comblist)
       
if __name__=="__main__":
    a=CSimpleCreaDiam("imp_OffAxis.txt")
    #print a.nonesiste
    a.mean_ff()
    
    

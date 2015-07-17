import namelist_class
import os
import math

'''2012/04/02 modified for web interface'''
    #list of notable thickness M and Q parameters
    XMM=(0.0016,-0.02,0)
    SX=(0.0006,-0.04,0)
    HXMT=(0.0006,0,0.2)
    HXMT_light=(0.0006,0,0.1)

class shellGeo(object):
    '''classe che regola il dimensionamento delle shell, a partire
    dai parametri. '''

    def __init__(self,focal,shellHeight,Dmax,fillingFactor,nshell,wallDensity,thickVal=XMM,profileType='ph'):
        '''thickval is a vector containing [thickM,thickQ,minThick]'''
        self.FL=focal
        self.SH=shellHeight
        self.Dmax=Dmax
        self.ff=fillingFactor
        self.nshell=nshell
        self.wd=wallDensity
        self.thickVal=thickVal
        self.profileType=profileType

    def read_pars_from_namelist(self,filename):
        '''Read all the geometrical parameters from a namelist in file FILENAME.'''
        shellGeo.__init__(self,self["focal_length"],self["shell_length"],self["maxDiam"],self["angular_shell_separation_deg"],
        self["nShells"],self["WallDensity"],(self["ThickM"],self["ThickQ"],self["minThick"]))

        l=open(filename,'r').readlines()
        l=[ll.strip() for ll in l if ll.count('=') == 1 ]  #consider only valid line (one and only one = sign in the line).
        dict={}
        for ll in l:
            s=ll.split('=')
            dict[s[0]]=s[1]
        return dict
        

    def diam (self):
        #tutte le misure in cm,
        #risultati convertiti in mm all'uscita
        ind= range(self.nshell)
        m = self.thickVal[0]
        q = self.thickVal[1]
        ds=[]
        spes=[]
        d=self.Dmax
        for i in ind:
            s=d/self.Dmax
            if s<self.thickVal[2]:s=self.thickVal[2]
            spes.append(s) #la lista degli spessori e' in mm
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

    def weights(self,diam,spes):
        ang=[math.atan(d/(2000*self.FL))/4 for d in diam]
        W=[]
        for i in range(len(diam)):
            Dmed=diam[i]
            Dmax=Dmed+2*math.tan(ang[i])*self.SH*10
            Dmin=Dmed-2*math.tan(3*ang[i])*self.SH*10
            p=((Dmax+Dmin)/2+Dmed)*math.pi*spes[i]*self.SH*self.wd/100000
            #print Dmax,Dmed,Dmin,spes[i]
            W.append(p)
        return W

    def fill (self,Dmin):
        '''returns the number of shell to fill the diamters larger than dmin, with
        the current configuration parameters by self'''
        store =self.nshell
        self.nshell=50
        n=self.nshell
        while n==self.nshell:
            self.nshell=self.nshell+1
            n=[j>=Dmin for j in self.diam2()[0]].count(True)            
        res=n
        self.nshell=store
        return res


       
class CSimpleCreaDiam(namelist_class.Namelist,shellGeo):
    '''classe per creazione di diametri con parametri letti da namelist,si inizializza con il nome del
    file delle impostazioni. manca qualche perfezionamento
    per gestire tutti i casi della namelist, come i vari flag o il diametro max non fornito direttamente'''
    ##attenzione - le variabili di shellGeo e Namelist non sono sincronizzate, quando cambio i valori di una
    ##non si rispecchiano sull'altra.
    
    def __init__(self,nl):
        namelist_class.Namelist.__init__(self,nl)
        #self.dirorigin=os.path.join(self["WORKDIR"].strip("\"' "),self["DIRORIGIN"].strip("\"' "))
        #self.dirrisult=os.path.join(self["WORKDIR"].strip("\"' "),self["DIRRISULT"].strip("\"' "))
        shellGeo.__init__(self,self["focal_length"],
                          self["shell_length"],self["maxDiam"],self["angular_shell_separation_deg"],
                          self["nShells"],self["WallDensity"],(self["ThickM"],self["ThickQ"],self["minThick"]))


def printPeso(spesLaw=CSimpleCreaDiam.HXMT,conf=None,NMaxSh=50):
    '''questa invece lo calcola leggendo una configurazione e variandone una serie di
    parametri passati in conf (lista di liste) con legge di spessori e il massimo numero
    di shell. restituisce stringa con riassunto'''

    if conf==None: conf=[[2.45,30,300,30],[2.65,30,300,30],[2.85,30,300,30],
          [2.45,30,330,30],[2.65,30,330,30],[2.85,30,330,30],
          [2.45,30,350,30],[2.65,30,350,30],[2.85,30,350,30]]    
    
    def setShell(a,v):
            F,H,Dmax,nsh=v
            a.FL=F
            a.SH=H
            a.Dmax=Dmax
            a.nshell=nsh

    def findNshell(a,Dmin,Nmax):
        '''calcolo di quante shell ci vogliono per coprire un
        range di diametri'''
        if Nmax!=-1:
            n=Nmax
            ntry=Nmax
            while (n==ntry and ntry<=Nmax):
                ntry=ntry+1
                a.nshell=ntry
                d,s=a.diam2()        
                Dmin=90
                n = [j>=Dmin for j in d].count(True)
            return min(n,Nmax)
        else: return a.nshell
        
    def printDrange(a):
        rv=[]
        d,s=a.diam2()
        w=a.weights(d,s)
        rv.append( "F=%s, H=%s, Dmax=%s, nshell=%s: "%(a.FL,a.SH,a.Dmax,a.nshell))
        rv.append( "diametri da %s a %s"%(d[-1],d[0]))
        rv.append( "spessori da %s a %s"%(s[-1],s[0]))
        rv.append( "Pesi: da %s a %s"%(w[-1],w[0]))
        rv.append( "Totale: %s\n"%reduce(lambda x, y: x+y, w))
        return "\n".join(rv)

    file=r"F:\vince\risultati\HXMT\dati\hxmt_9\hxmtStart_focals\f245\imp_offAxis.txt"
    a=CSimpleCreaDiam(file)

    Dmin=90
    Nmax=NMaxSh
    strv=[]
    for v in conf:
        setShell(a,v)
        a.thickVal=(spesLaw)
        ns=findNshell(a,Dmin,Nmax)
        strv.append( "Dmin=%s, ci vogliono %s shell,"%(Dmin,ns))
        strv.append( "infatti:")
        a.nshell=ns-1
        strv.append(printDrange(a))
    return "\n".join(strv)


def NshellPerD (listaFile):
    '''analizza una serie di configurazioni a partire dai file delle impostazioni contenuti nelle
    cartelle elencate nel file listaFile.
    ne modifica il diametro massimo e calcola il numero di shell per lo stesso diametro minimo.
    genera i nuovi file di impostazioni per il raytracing. serve per raggruppare piu' configurazioni
    in un singolo raytracing per avere plot incrementale'''

    outDir="inc_conf"
    confs=open(listaFile,"r").readlines()
    if not os.path.isdir(outDir): os.mkdir(outDir)
    count=1
    for c in confs:
        confName,Dmin=c.split()[0:2]
        wolter=CSimpleCreaDiamOA(os.path.join(confName,"imp_OffAxis.txt"))
        wolter.Dmax=700
        wolter["massimo"]=700
        wolter["nshDaImp"]=wolter.fill(float(Dmin))
        wolter.write(os.path.join(outDir,confName+"_imp.txt"))
        open(os.path.join(outDir,confName+"_imp.txt"),"a").write(str(wolter["nshDaImp"]))
        print "scritti i risultati per la directory %s,\nanalizzati %s cartelle" %(confName,count)
        print "servono %s shell"%wolter["nshDaImp"]
        count=count+1


        
if __name__=="__main__":
    os.chdir(r"E:\work\WTDf")
    a=CSimpleCreaDiamOA("imp_OffAxis.dat")
    
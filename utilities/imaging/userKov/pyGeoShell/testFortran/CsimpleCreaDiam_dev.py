import namelist_class
import os
import math

class shellGeo:
    '''classe che regola il dimensionamento delle shell, a partire
    dai parametri. '''
    XMM=(0.0016,-0.02)
    SX=(0.0006,-0.04)
    def __init__(self,focal,shellHeight,Dmax,fillingFactor,nshell,wallDensity,thickVal=XMM):
        self.FL=focal
        self.SH=shellHeight
        self.Dmax=Dmax
        self.ff=fillingFactor
        self.nshell=nshell
        self.wd=wallDensity
        self.thickVal=thickVal

    def diam (self):
        #tutte le misure in cm,
        #risultati convertiti in mm all'uscita
        ind= range(self.nshell)
        m = self.thickVal[0]
        q = self.thickVal[1]
        ds=[]
        spes=[]
        d=self.Dmax
        if d==60:ee=0.3
        if d==70:ee=0.36
        for i in ind:
            s=ee*d/self.Dmax
            spes.append(s) #la lista degli spessori e' in mm
            #print (d*10*m + q),spes[-1]
            d=d-2*(spes[-1]/10)
            ang= math.atan2(d,(2*self.FL*100))/4
            dmed=d-2*self.SH*math.tan(ang)
            ds.append(dmed*10)
            d=dmed-2*self.SH*math.sin(math.radians(self.ff))
        return ds, spes    
        
class CSimpleCreaDiam(namelist_class.Namelist,shellGeo):
    '''classe per creazione di diametri con parametri letti da namelist, manca qualche perfezionamento
    per gestire tutti i casi della namelist, come i vari flag o il diametro max non fornito direttamente'''
    
    def __init__(self,nl):
        namelist_class.Namelist.__init__(self,nl)
        self.dirorigin=os.path.join(self["WORKDIR"].strip("\"' "),self["DIRORIGIN"].strip("\"' "))
        self.dirrisult=os.path.join(self["WORKDIR"].strip("\"' "),self["DIRRISULT"].strip("\"' "))
        shellGeo.__init__(self,self["F_LENGTHdaImp_m"],
                          self["F_HEIGHTdaImp_cm"],self["massimo"],self["FieldOfViewDeg"],
                          self["nshDaImp"],self["WallDensity"],(self["ThickM"],self["ThickQ"]))
        
       
nl2=CSimpleCreaDiam("imp_OffAxis.txt")

focal=[i*2.5 for i in range(8,12)]
ff=[0.00,0.07,0.15]
D=[60,70]
comblist=[[f1,f2,f3] for f1 in focal for f2 in ff for f3 in D]

#comblist=[[20.0,0.0,60]]
p=open("confrDiam.plt","w")
#p2=open("confrDiam2.plt","w")
#p2.write("plot ")
for comb in comblist:
    ndf="F"+str(comb[0]).split(".")[0]+"_ff"+str(comb[1]).split(".")[1].zfill(3)+"_D"+str(comb[2])+".txt"
    nl2.FL,nl2.ff,nl2.Dmax=comb
    di,sp=nl2.diam2()
    f= open(ndf,"w")
    for d in di[-1:0:-1]:
        f.write(str(d)+"\n")
    print ndf,nl2.FL,nl2.SH,nl2.Dmax,nl2.ff,nl2.nshell,nl2.wd,(nl2.thickVal[0],nl2.thickVal[1])
    print "\n"
    f.close()
    #per plottare confronto con vecchi diametri
    p.write("plot '"+ndf+"' u 0:1 w lp")
    fcnf1="mlff"+str(comb[1]).split(".")[1].zfill(3)+str(comb[0]).split(".")[0]+"mmax"+str(comb[2])+".txt"
    fcnf2="m2ff"+str(comb[1]).split(".")[1].zfill(3)+str(comb[0]).split(".")[0]+"mmax"+str(comb[2])+".txt"
    if os.path.exists(os.path.join("af_files",fcnf1)):
       p.write(",\\\n'"+os.path.join("af_files",fcnf1)+"' u 0:1 w lp")
    if os.path.exists(os.path.join("af_files",fcnf2)):
       p.write(",\\\n'"+os.path.join("af_files",fcnf2)+"' u 0:1 w lp")
    if os.path.exists(os.path.join("frlike",ndf)):
       p.write(",\\\n'"+os.path.join("frlike",ndf)+"' u 0:1 w lp")
    fcnf3="ml_ff"+str(comb[1]).split(".")[1].zfill(3)+"_"+str(comb[0]).split(".")[0]+"m_max"+str(comb[2])+".txt"
    if os.path.exists(os.path.join("extractedDiam",fcnf3)):
       p.write(",\\\n'"+os.path.join("extractedDiam",fcnf3)+"' u 0:1 w lp")
    fcnf4="ir_ff"+str(comb[1]).split(".")[1].zfill(3)+"_"+str(comb[0]).split(".")[0]+"m_max"+str(comb[2])+".txt"
    if os.path.exists(os.path.join("extractedDiam",fcnf4)):
       p.write(",\\\n'"+os.path.join("extractedDiam",fcnf4)+"' u 0:1 w lp")
    p.write("\n")
    p.write("pause -1\n\n")

    p.write("set term png\n")
    p.write("set out '"+ndf+".png'\n")
    p.write("replot\n")
    p.write("set term win\n set out\n\n")

    #p2.write("'"+ndf+"' u 0:1 w lp,\\\n")    

p.close()    
#p2.close()        
    


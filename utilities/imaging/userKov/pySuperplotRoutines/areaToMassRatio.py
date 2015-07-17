from numpy import *
from userKov.pyGeneralRoutines.generalRoutines import loadCol
import os


def areaToMassRatio(folder,enTarget,npunti,coatInd=None):

    if coatInd == None: coatInd=0
    coatingPars=[[1,3],['M2','Pt']] #la prima lista sono gli indici delle colonne da leggere,
                                    #la seconda le stringhe da attaccare al nome del file.
    coatCol=coatingPars[0][coatInd]
    coatStr=coatingPars[1][coatInd]
    
    #recupera i pesi
    ff=open(os.path.join(folder,"pesi.txt"),"r")
    pes=map (float,ff.readlines()[4:])
    pes.sort()
    ff.close()
    nshells=len(pes)

    #recuperai diametri
    diam=loadCol(os.path.join(folder,'diamplot.txt'),1,skip=1,nlines=nshells)
    diam.sort()
    
    #recupera le aree
    ener=loadCol(os.path.join(folder,'001','aeff0001.dat'),1,nlines=npunti)
    aree=[loadCol(os.path.join(folder,'%03i'%(i),'aeff%04i.dat'%(i)),coatCol+1,nlines=npunti) for i in range(1,nshells+1)]
    aree=array([interp(enTarget,ener,a) for a in aree])
    ratio=array([aa/pes for aa in aree.transpose()])
    for i,en in enumerate(enTarget):
        outstr=["Nshell\tdiam\tAreaToMass\tArea\tMass\ten=%s"%en]
        for ns,(d,r) in enumerate(zip(diam,ratio[i,:])):
            outstr.append("%s\t%s\t%s\t%s\t%s"%(nshells-ns,d,r,aree[ns,i],pes[ns]))
        open(os.path.join(folder,"%03i_%sarea.txt"%(en,coatStr)),"w").write("\n".join(outstr))

if __name__=="__main__":
##    #directory con i risultati di superplot, mettere dentro anche file dei pesi
##    folder=r"F:\Nuova\src_plt\hexitSat2009\fullsizeD450_F10_1plt\upperlimit"
##    enTarget=[10.,20.,30.,50.,60.,70.]
##    npunti=800  #per evitare errori con stringa alla fine
##    areaToMassRatio(folder,enTarget,npunti,0)
    
    #directory con i risultati di superplot, mettere dentro anche file dei pesi
    folder=r"E:\work\usb_F\Nuova\src_plt\NHXM_descoping\modulated_PtC50_1plt\upperlimit"
    enTarget=[10.,20.,30.,50.,60.,70.]
    npunti=1000  #per evitare errori con stringa alla fine
    areaToMassRatio(folder,enTarget,npunti,0)    
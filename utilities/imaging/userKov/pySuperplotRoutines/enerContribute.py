from numpy import *
from userKov.pyGeneralRoutines.generalRoutines import loadCol


def enerContribute(folder,entarget,outdir="",coating=0,weightDir=""):
    ''' given a <folder> with the results of superplot (to read the file diamplot.txt)
    writes in a file with suffix _contrib the single shell contributes for the energy
    entarget. coating: 0=multilayer, 1=single layer
    '''
    
    if coating!=0: coating=1 #per selezionare la giusta colonna per il tipo di coating
    coatingStrings=["ML","sl"]
    coatingCol=[0,2] #colonne corrispondenti alla coatingString
    diam=loadCol(folder+'diamplot.txt',0,skip=1)
    acoll=loadCol(folder+'diamplot.txt',1,skip=1)
    angle=loadCol(folder+'diamplot.txt',2,skip=1)
    if weightDir <> "":
        peso=loadCol(weightDir+'pesi.txt',0,skip=4)
    else:
        peso=zeros(len(diam))
    if (diam[0]<diam[-1] and peso[0]>peso[-1]) or (diam[0]<diam[-1] and peso[0]>peso[-1]) : peso=peso[::-1]
    aeff30=[]
    for i in range(0,len(angle)):
            #read the right column for mono/multilayer
        aeff=loadCol(folder+'%003i\\aeff%0004i.dat'%(i+1,i+1),1+coatingCol[coating],nlines=-5) 
        ener=loadCol(folder+'%003i\\aeff%0004i.dat'%(i+1,i+1),0,nlines=-5)
        aeff30.append(interp([entarget],ener,aeff))
    if outdir=="":outdir=folder
    tot=sum(aeff30)
    s=["nshell\tdiam(mm)\tangle(rad)\tacoll(cm^2)\taeff(cm^2)at %s\taeff contrib.\tweight(kg)"%entarget]
    for i in range(0,len(angle)):s.append("%s\t%s\t%s\t%s\t%s\t%s\t%s"%(i+1,diam[i],angle[i],acoll[i],float(aeff30[i]),float(aeff30[i])/tot,peso[i]))
    f=open(outdir+"%s_contrib_%3.3i.dat"%(coatingStrings[coating],entarget),'w').write("\n".join(s))
    

if __name__=="__main__":
##    folder1=r'D:\work\Nuova\src_plt\hexitSat2009\F10D350ff010_thsx_1plt\upperlimit\\'
##    #r'D:\work\Nuova\src_plt\hexitSat2009\F10D350ff005_100sh_thsx_1plt\upperlimit\\'    
##    entarget=[1.,30.,60.,70.]
    folder1=r'E:\work\usb_F\Nuova\src_plt\NHXM_descoping\modulated_mlphB_1plt\upperlimit\\'
    entarget=[1.,20.,35.,60.,70.]
    for en in entarget: enerContribute(folder1,en,coating=1,weightDir=folder1)
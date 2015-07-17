from userKov.pyGeneralRoutines.generalRoutines import loadCol
from userKov.pyOffAxis.extractInfo import loadEn
import os

spiderFac=0.9
folder=r"H:\vince\prog\f_2005\ufficiali\prog_superplotGen1\simx2008_structuresPtC\PtC_4sum\ml115_200_oc\somma"
listaFile="lista.txt"
os.chdir(folder)

for f in open(listaFile,"r").read().split('\n'):
    if f.strip() != '':
        f=f.strip()
        aree=loadCol(f,2)
        ener=loadEn(f)
        name=os.path.split(f)[0]
        out=open(name.split(os.path.sep)[0]+"_spider.txt","w")
        a=aree[:]
        for k in range(len(a)):
                a[k]=[area*spiderFac for area in a[k]]
        for aa in a:
            for t in zip(ener,aa):out.write(str(t[0])+"\t"+str(t[1])+"\n")
            out.write("\n"*2)        
        out.close()        
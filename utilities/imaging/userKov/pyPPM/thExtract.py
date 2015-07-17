from userKov.pyGeneralRoutines.generalRoutines import *
import sys
from numpy import *
import os

if __name__=="__main__":
    if (len(sys.argv)==1):
        filename= raw_input('give a filename')
    else:
        filename=sys.argv[1]
        
    stTh=loadCol(filename,0)[0]
    stTh=stTh[1:]  #remove substrate
    if (len(stTh)!=int(len(stTh)/2)*2):
        stTh.append(0)
    stTh=array(stTh)
    #crea dati su due colonne e li scrive su file
    stTh=reshape(stTh,(-1,2))
    outfile="_2col".join(os.path.splitext(filename))
    of=open (outfile,'w')
    for l in stTh:
        of.write("%s\t%s\n"%tuple(l))
        #print l
    of.close()
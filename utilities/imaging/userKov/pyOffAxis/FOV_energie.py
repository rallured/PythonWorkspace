# -*- coding: latin-1 -*-
'''calcola fov, grasp e grasp/platescale per i file (risultanti dal programma fortran offAxis)
elencati nei file di testo in lista fl. i risultati sono nei file di testo con 4 colonne:
energia, fov diameter(arcmin), grasp (FOV**2/Aeff**2),grasp / plateScale(Rspot(=15")/Focal).
attenzione alla definizione della griglia di punti per l'output in entarget.
I valori sono calcolati fino ad un energia massima per la quale si ha Aeff(Emax)=cutoffEn*Aeff[E0],
i valori di taglio vengono elencati in cutoffEnergy.txt.
Il programma crea un file per plottare le singole FOV con peso e area a 30 kev, per plottare 
gli altri risultati usare le routines in plotRoutines.py

Manca un fattore pigreco/4 nel calcolo del grasp'''

import os
from userKov.pyOffAxis import extractInfo
from userKov.pyGeneralRoutines.generalRoutines import loadCol
from numpy import interp
#from plotRoutines import *

def calcolaFOV_Grasp(resdir,cutoffRatio=0.05,fileLista='lista.txt',Rspot=15,focalm=10.):
    '''calcola i risultati di FOV e GRASP per le cartelle (risultati di offaxis)
    elencate nel file fl. i risultati vanno in resdir
    cutoffRatio rapporto tra aeffMax e aeff(E) al di sotto del quale si ritiene
    non affidabile il calcolo del FOV'''
    #------------------------------------
    #parametri
    npoints=100
    filelista='lista.txt'
    os.chdir(resdir)
    entarget= [(80.)/npoints*(i+1) for i in range(npoints)]
    entarget2=30.
    Wfactor=1.3
    SpiderFactor=0.9
    taroccoFactor=1  #fattore di moltiplicazione per il FOV per taroccare i dati
    #----------------------------------------------
    #inizializzazione
    if not os.path.exists(resdir): os.mkdir(resdir)
    plt=open(os.path.join(resdir,"FOVvsEnergy.plt"),"w")
    cut=open(os.path.join(resdir,"cutoffEnergy.txt"),"w")
    tab=open(os.path.join(resdir,"tabella.txt"),"w")
    tabella=[]
    f=open(os.path.join(resdir,"note.txt"),"w")
    f.write(__doc__)
    f.close()
    print __doc__
    #----------------------------------------------
    
    listaDir=open(filelista,"r")
    a=listaDir.read().split()
    for config in a:
        print config
        ener=extractInfo.loadEn(os.path.join(config,"aree.txt"))   
        aree=loadCol(os.path.join(config,"aree.txt"),1)[0] #recupera vettore delle aree in asse
        m=max(aree)
        cutoffEn=ener[[i>m*cutoffRatio for i in aree].count(True)-1]
        cut.write(config+"\t"+str(cutoffEn)+"\n")
        entargetCut=entarget[0:[i<=cutoffEn for i in entarget].count(True)-1]
        FOV_vec_diam=[2*taroccoFactor*extractInfo.offAxisInfo(config,e)[2] for e in entargetCut]
        Weight,ang,FOV30,onAx,fovAr=extractInfo.offAxisInfo(config,entarget2)
        tabella.append([config,Weight*Wfactor,FOV30*2*taroccoFactor,onAx*SpiderFactor,ang])
            
        outf=open(os.path.join(resdir,config+"_grasp.txt"),"w")
        grasp=[(FOV_vec_diam[i]**2)*interp(entargetCut,ener,aree)[i] for i in range(len(FOV_vec_diam))]
        grasp2=[float(focalm)*g/Rspot for g in grasp]
        for i in range (len(entargetCut)):
            outf.write(str(entarget[i])+"\t"+str(FOV_vec_diam[i])+"\t"+str(grasp[i])+"\t"+str(grasp2[i])+"\n")
        outf.close()
        
        plt.write("set title ' %s, ang. shell sep.=%s, FL= %s m, Dmax= %s cm' \n"  %tuple(extractInfo.extractTitle(config)))
        plt.write("set term win\n")
        plt.write("set out\n")
        plt.write("set key box\n")
        plt.write("set grid\n")
        plt.write("set xrange [0:80]\n")
        plt.write("set yrange [0:18]\n")
        plt.write("unset label\n")
        plt.write("set xlabel 'Energy (keV)'\n")
        plt.write("set ylabel 'Interpolated FOV diameter (arcmin)'\n")
        label='set label "WEIGHT (+30%%): %4.2f Kg\\n'%(Weight*Wfactor)
        label=label+ 'FOV@%4.2f keV: %4.2f arcmin\\n' %(entarget2,FOV30*2*taroccoFactor)
        label=label+ 'On-axis Aeff@%4.2f keV: %4.2f cm^2" at graph 0.02,0.15\n' %(entarget2,onAx*SpiderFactor)
        plt.write(label)        
        plt.write("plot '" + config + "_grasp.txt' u 1:2 w l title '"+config+"'\n")
        
        plt.write("pause -1\n")
        plt.write("set term png\n")
        plt.write("set output 'FOVvsEn"+config+".png'\n")
        plt.write("replot\n")
        plt.write("set term win\n")
        plt.write("set out\n")
        plt.write("\n")
    listaDir.close()

    plt.close()
    cut.write("Aeff cutoff Ratio: "+str(cutoffRatio))
    cut.close()

    tab.write("Coating\tangular shell sep (°) \tFocal lenght(m)\tMax diameter(cm)\t")
    tab.write("Weight(Kg)with structure\tFOV diam(arcmin)\tArea on Axis(cm^2)x%s\tSampling Angles(arcmin)\n"%SpiderFactor)
    for el in tabella:
        for e in extractInfo.extractTitle(el[0]): tab.write (e+"\t")
        for dato in el[1:-1]: tab.write(str(dato)+"\t")
        tab.write(str(el[-1])+"\n")
    tab.close() 
 

if __name__=="__main__":
    #resdir="FOV_energy"
    
    #calcolaFOV_Grasp(resdir,0.05,"80kev.txt")
    
    #plottaGrasp(resdir,(0,2,3))
    #plottaFOV(resdir,(0,2,3))
    #plottaGrasp(resdir,(0,1))
    
    resdir=r'E:\work\workOA\traie8\NHXMphB_mlArea'
    calcolaFOV_Grasp(resdir,0.05,"FOV.txt")
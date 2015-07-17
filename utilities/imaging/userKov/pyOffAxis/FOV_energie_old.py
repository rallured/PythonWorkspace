# -*- coding: latin-1 -*-

import generaPlot
import os

def interpola (xsource,ysource,xtarget):
    '''dati vettori x e y e vettore xtarget delle ascisse per cui si vogliono
    valori della y interpolata, restuisce il vettore delle y interpolate'''
    indici=[[i<target for i in xsource].count(True)-1 for target in xtarget]
    y=[]
    for i,e in enumerate(xtarget):
        j=indici[i]
        y.append(ysource[j]+(e-xsource[j])*(ysource[j+1]-ysource[j])/(xsource[j+1]-xsource[j]))
    return y

if __name__=="x__main__":
    
    resdir="fov_test"
    cutoffRatio=0.05 #rapporto tra aeffMax e aeff(E) al di sotto del quale
                    #si ritiene non affidabile il calcolo del FOV
    entarget= [0.1+((60.-0.1)/99*i) for i in range(100)]
    entarget2=30.
    if not os.path.exists(resdir): os.mkdir(resdir)
    plt=open(os.path.join(resdir,"FOVvsEnergy.plt"),"w")
    cut=open(os.path.join(resdir,"cutoffEnergy.txt"),"w")
    fl=['IrTutti.txt','tutti.txt']
    #fl=['testIr.txt'] 
    tabella=[]

    for filelista in fl:
        listaDir=open(filelista,"r")
        a=listaDir.read().split()
        for config in a:
            print config
            ener=generaPlot.loadEn(os.path.join(config,"aree.txt"))   
            aree=generaPlot.loadCol(os.path.join(config,"aree.txt"),2)[0]
            Area30=aree[[i<entarget2 for i in ener].count(True)-1]
            m=max(aree)
            cutoffEn=ener[[i>m*cutoffRatio for i in aree].count(True)-1]
            cut.write(config+"\t"+str(cutoffEn)+"\n")
            entargetCut=entarget[1:[i<cutoffEn for i in entarget].count(True)-1]
            FOV_vec=[generaPlot.offAxisInfo(config,e)[2] for e in entargetCut]
            Weight,ang,FOV30,fovAr=generaPlot.offAxisInfo(config,entarget2)
            tabella.append([config,Weight,FOV30,fovAr,ang])    
                
            outf=open(os.path.join(resdir,config+".txt"),"w")
            graspf=open(os.path.join(resdir,config+"_grasp.txt"),"w")
            grasp=[(FOV_vec[i]**2)*interpola(ener,aree,entargetCut)[i] for i in range(len(FOV_vec))]
            for i in range (len(entargetCut)):
                outf.write(str(entarget[i])+"\t"+str(FOV_vec[i])+"\n")
                graspf.write(str(entarget[i])+"\t"+str(grasp[i])+"\n")
            outf.close()
            graspf.close()
            
            plt.write("set title ' %s, ang. shell sep.=%s, FL= %s m, Dmax= %s cm' \n"  %tuple(generaPlot.extractTitle(config)))
            plt.write("set term win\n")
            plt.write("set out\n")
            plt.write("set key box\n")
            plt.write("set grid\n")
            plt.write("unset label\n")
            plt.write("set xlabel 'Energy (keV)'\n")
            plt.write("set ylabel 'Interpolated FOV diameter (arcmin)'\n")
            label='set label "WEIGHT (+30%%): %4.2f Kg\\n'%(Weight*1.3,)
            label=label+ 'FOV@%4.2f keV: %4.2f arcmin\\n' %(entarget2,FOV30*2)
            label=label+ 'On-axis Aeff@%4.2f keV: %4.2f cm^2" at graph 0.02,0.15\n' %(entarget2,Area30)
            plt.write(label)        
            plt.write("plot '" + config + ".txt' u 1:($2*2) w l title '"+config+"'\n")
            
            plt.write("pause -1\n")
            plt.write("set term png\n")
            plt.write("set output '"+config+".png'\n")
            plt.write("replot\n")
            plt.write("set term win\n")
            plt.write("set out\n")
            plt.write("\n")
        listaDir.close()

    plt.close()
    cut.write("Aeff cutoff Ratio: "+str(cutoffRatio))
    cut.close()

    tab=open(os.path.join(resdir,"tabella.txt"),"w")
    tab.write("Coating\tangular shell sep (°) \tFocal lenght(m)\tMax diameter(cm)\t")
    tab.write("Weight(Kg w/o struct)\tFOV(arcmin radius)\tArea(cm^2 @FOVangle)\tSampling Angles(arcmin)\n")
    for el in tabella:
        for e in generaPlot.extractTitle(el[0]): tab.write (e+"\t")
        for dato in el[1:-1]: tab.write(str(dato)+"\t")
        tab.write(str(el[-1])+"\n")
    tab.close()
    plotta(resdir)  
            
def plottaInsieme(p,resdir,ylabel="Interpolated FOV diameter (arcmin)",dataFileFrmt="%s_ff%03i_%im_max%s.txt",
                   outFrmt="FOV_%s_%i.png",term="png",using="u 1:($2*2)",xrange="[0:60]",yrange="[0:18]",
                   comb=[["ML","Ir"],[20,22.5,25,27.5,30],[0.0,0.07,0.15],[60,70]]):
    '''comb in formato [[coatings],[focali],[ff],[diam]]'''
    
    titFrmt="ff=%s°, Dmax=%scm"
    p.write("unset label \n")
    p.write("set xlabel 'Energy (keV)'\n")
    coatdic={"ML":"W/Si multilayer", "Ir":"Ir monolayer"}
    #print comb
    for coat in comb[0]:
        for focal in comb[1]:
            print coat,focal,ylabel
            combs= [(ff,d) for ff in comb[2] for d in comb[3]]
            combs=[cc for cc in combs if os.path.exists(os.path.join(resdir,dataFileFrmt%(coat,cc[0]*100,focal,cc[1])))]
            print len(combs)
            if len(combs)==0: continue
            p.write ("set yrange %s\n"%yrange)
            p.write ("set yrange %s\n"%xrange)
            p.write("set title '%s, %sm Focal length'\n" %(coatdic[coat],focal))
            p.write("set ylabel '%s'\n" %(ylabel)) 
            p.write("plot ")
            for cc in combs[0:-1]:
                p.write("'%s' %s title '%s' w l,\\\n"%(dataFileFrmt%(coat,cc[0]*100,focal,cc[1]),using,titFrmt%cc)) 
            cc=combs[-1]
            p.write("'%s' %s title '%s' w l \n\n"%(dataFileFrmt%(coat,cc[0]*100,focal,cc[1]),using,titFrmt%cc)) 
            p.write("pause -1\n")
            p.write("set term %s\n"%(term))
            p.write("set output '%s'\n"%(outFrmt%(coat,focal))) 
            p.write("replot\n")
            p.write("set term win\n")
            p.write("set out\n")
            p.write("\n\n")
            p.flush()


def ciclazza (possibilita, indici, base=None):
    '''data una lista possibilita di liste di elementi tra cui scegliere, restituisce
    tutte le possibili scelte ciclando sugli indici nella tupla indici.
    per gli elementi non ciclati usa i valori in base. Senza base usa 0, mentre
    passando possibilita come base li lascia immutati'''
    
    combs=[]
    indici=tuple(indici)
    if base==None:
        b=[[0]*(len(possibilita))]
    else:
        b=[base]
    for i in indici:
        for bb in b:                    #print "ricavo nuovi valori da bb ",bb
            for x in possibilita[i]:
                bbb=bb[:]
                bbb[i]=x                #print "aggiungo ",bbb," a combs: ",combs
                combs.append(bbb)
            b=combs[:]
        combs=[]
    return b


    
def groupPlot(resdir,pltFilePrefix="",ylabel="Interpolated FOV diameter (arcmin)",
              dataFileFrmt="%s_ff%03i_%im_max%s.txt",
              outPrefix=("FOV_",".png"),term="png",using="u 1:($2*2)",xrange="[0:60]",yrange="[0:18]",
              comb=[["ML","Ir"],[20,22.5,25,27.5,30],[0.0,0.07,0.15],[60,70]],
              combLbl=["Coat","Focal","FF","Dmax"],indice=(0,1),):
    '''comb in formato [[coatings],[focali],[ff],[diam]], indice 0, 1 o 2, a seconda di come raggruppare'''

    pltFileName=pltFilePrefix+"by"+"".join([combLbl[i] for i in indice])+".plt"
    p=open(os.path.join(resdir,pltFileName),"w")   
    indice=tuple(indice)
    coatdic={"ML":"W/Si multilayer", "Ir":"Ir monolayer"}
    cycVars=[i for i in range(4) if i not in indice]
    lblFrmt=", ".join([["%s","FL=%sm","ff=%s°","Dmax=%scm"][i] for i in cycVars])
    titleFrmt=", ".join([["%s","FL=%sm","ff=%s°","Dmax=%scm"][i] for i in indice])
    outFrmt=outPrefix[0]+"_".join([["%s","%sm","ff%s","D%s"][i] for i in indice])+outPrefix[1]
    p.write("unset label \n")
    p.write("set xlabel 'Energy (keV)'\n")
        
    extCyc=ciclazza(comb,indice,comb)   #crea conbinazioni esterne 
    i=0
    for cyc in extCyc:
        i=i+1
        combs= ciclazza(cyc,cycVars,cyc)    #crea conbinazioni interne
        combs= [cc for cc in combs if os.path.exists(os.path.join(resdir,dataFileFrmt%(cc[0],cc[2]*100,cc[1],cc[3])))]
        if len(combs)==0: continue
        p.write ("set yrange %s\n"%yrange)
        p.write ("set xrange %s\n"%xrange)
        p.write("set title '%sm Focal length'\n" %(titleFrmt%(cyc[indice[0]],cyc[indice[1]])))
        p.write("set ylabel '%s'\n" %(ylabel)) 
        p.write("plot ")
        for cc in combs[0:-1]:
            p.write("'%s' %s title '%s' w l,\\\n"%(dataFileFrmt%(cc[0],cc[2]*100,cc[1],cc[3]),
                                                   using,lblFrmt%(cc[cycVars[0]],cc[cycVars[1]]))) 
        cc=combs[-1]
        p.write("'%s' %s title '%s' w l \n\n"%(dataFileFrmt%(cc[0],cc[2]*100,cc[1],cc[3]),
                                               using,lblFrmt%(cc[cycVars[0]],cc[cycVars[1]])))
        p.write("pause -1\n")
        p.write("set term %s\n"%(term))
        p.write("set output '%s'\n"%(outFrmt%(cyc[indice[0]],cyc[indice[1]]))) 
        p.write("replot\n")
        p.write("set term win\n")
        p.write("set out\n")
        p.write("\n\n")
        p.flush()
    p.close() 

            
def plotta(resdir):
    
    p=open(os.path.join(resdir,"byFocal.plt"),"w")
    plottaInsieme(p,resdir)
    plottaInsieme(p,resdir,ylabel="Grasp (Aeff x FOV^2)",dataFileFrmt="%s_ff%03i_%im_max%s_grasp.txt",
                   outFrmt="Grasp_%s_%i.png",using="u 1:2",yrange="[0:120000]")
    p.close()
    print "fine"

def testa(resdir):
    #plotta(resdir)
    '''resdir="fov_test2"
    groupPlot(resdir,filePrefix="grasp_")'''
    groupPlot(resdir,pltFilePrefix="Grasp_",ylabel="Grasp (Aeff x FOV^2)",dataFileFrmt="%s_ff%03i_%im_max%s_grasp.txt",
                   outPrefix=("Grasp_",".png"),using="u 1:2",yrange="[0:120000]",indice=(0,2))  

    '''possibilita=[[1,2,3],[4,5],[6,7,8],[9,10]]
    indici=(0,)
    print ciclazza(possibilita,indici)'''

if __name__=="__main__":
    testa(resdir)
# -*- coding: latin-1 -*-
import userKov
from userKov.pyGeneralRoutines import  generalRoutines
from userKov.pyOffAxis import FOV_energie
from extractInfo import Wrecover
import os
    
def groupPlot(resdir,indice=(0,1),xrange="[*:*]",yrange="[*:*]",
              comb=[["ML","Ir"],[20,22.5,25,27.5,30],[0.0,0.07,0.15],[60,70]],
              pltFilePrefix="",outPrefix="FOV_",term="",combLbl=["Coat","Focal","FF","Dmax"],
              dataFileFrmt="%s_ff%03i_%im_max%s.txt",using="u 1:2",
              xlabel="Energy (keV)",ylabel="Interpolated FOV diameter (arcmin)"):
    '''plotta i file, raggruppandoli per caratteristiche comuni, secondo impostazioni
    resdir: directory per i risultati,pltFilePrefix: prefisso da attaccare al file di plot
    xrange,yrange: range del plot
    comb in formato [[coatings],[focali],[ff],[diam]],
    indice 0, 1 o 2, a seconda di quali elementi tenere raggruppati
    outPrefix: prefisso per eventuale figura di output
    term: terminale per figura di output ("png"|"ps"|"bmp"|"")
    combLbl: nome dei campi corrispondenti agli elementi di comb. usato per generare il nome del file di plot
    dataFileFrmt: formato del nome del file da cui leggere i dati
    using: colonne da plottare
    xlabel,ylabel: etichette degli assi'''

    termDic={"png":("png",".png"),"ps":("postscript color",".eps"),"bmp":("bitmap",".bmp"),"":("",""),"win":("","")}    
    termExt=termDic[term]
    pltFileName=pltFilePrefix+"by"+"".join([combLbl[i] for i in indice])+".plt"
    p=open(os.path.join(resdir,pltFileName),"w")   
    indice=tuple(indice)
    coatDic={"M2":"Pt/C multilayer","ML":"W/Si multilayer", "Ir":"Ir monolayer"}
    cycVars=[i for i in range(4) if i not in indice]
    lblFrmt=", ".join([["%s","FL=%sm","ff=%s°","Dmax=%scm"][i] for i in cycVars])
    titleFrmt=", ".join([["%s","FL=%sm","ff=%s°","Dmax=%scm"][i] for i in indice])
    outFrmt=outPrefix+"_".join([["%s","%sm","ff%s","D%s"][i] for i in indice])+termExt[1]
        
    p.write("unset label \n")
    p.write("set xlabel 'Energy (keV)'\n")
        
    extCyc=generalRoutines.ciclazza(comb,indice,comb)   #crea conbinazioni esterne
    i=0
    for cyc in extCyc:
        i=i+1
        combs= generalRoutines.ciclazza(cyc,cycVars,cyc)    #crea conbinazioni interne
        combs= [cc for cc in combs if os.path.exists(os.path.join(resdir,dataFileFrmt%(cc[0],cc[2]*100,cc[1],cc[3])))]
        if len(combs)==0: continue
        d=range(len(comb))
        for i,dd in enumerate(cyc):
            if i in cycVars:
                d[i]=cyc[i][0]
            else:
                d[i]=cyc[i]
        d=[coatDic[d[0]],d[1],d[2],d[3]]
        titVars=tuple([d[i] for i in indice])
        p.write ("set yrange %s\n"%yrange)
        p.write ("set xrange %s\n"%xrange)
        p.write("set title '%s'\n" %(titleFrmt%(titVars)))
        p.write("set ylabel '%s'\n" %(ylabel)) 
        p.write("plot ")
        for cc in combs[0:-1]:
            d=[coatDic[cc[0]],cc[1],cc[2],cc[3]]
            lblVars=tuple([d[i] for i in cycVars])
            p.write("'%s' %s title '%s' w l lw 2,\\\n"%(dataFileFrmt%(cc[0],cc[2]*100,cc[1],cc[3]),
                                                       using,lblFrmt%(lblVars))) 
        cc=combs[-1]
        d=[coatDic[cc[0]],cc[1],cc[2],cc[3]]
        lblVars=tuple([d[i] for i in cycVars])
        p.write("'%s' %s title '%s' w l lw 2 \n\n"%(dataFileFrmt%(cc[0],cc[2]*100,cc[1],cc[3]),
                                               using,lblFrmt%(lblVars)))
        p.write("pause -1\n")
        if term[0]<>"":
            p.write("set term %s\n"%(termExt[0]))
            p.write("set output '%s'\n"%(outFrmt%(tuple([cyc[i] for i in indice])))) 
            p.write("replot\n")
            p.write("set term win\n")
            p.write("set out\n")
        p.write("\n\n")
        p.flush()
    p.close() 


def cmpPlot(resdir,groupFile="gruppi.txt",overwrite=0,
            xrange="[*:*]",yrange="[*:*]",
            pltFilePrefix="",outPrefix="FOV_",term="",
            dataFileFrmt="%s.txt",using="u 1:2",
            title="",
            xlabel="Energy (keV)",ylabel="Interpolated FOV diameter (arcmin)"):
    '''plotta risultati elencati in groupFile per gruppi (separati riga bianca). se in numero maggiore
    di nmax(=6) suddivide in gruppi di peso crescente con differenze di peso max secondo quanto impostato
    in DwMax nella routine dividi.
    resdir: directory per i risultati,pltFilePrefix: prefisso da attaccare al file di plot
    xrange,yrange: range del plot
    pltFilePrefix: prefisso per il file plt generato
    outPrefix: prefisso per eventuale figura di output
    term: terminale per figura di output ("png"|"ps"|"bmp"|"")
    dataFileFrmt: formato del nome del file da cui leggere i dati (N.B.:non funziona, perche' il file no
    e' stato ancora creato quando viene usato, funziona solo se si rilancia una seconda volta, cosi' e' inutile!
    using: colonne da plottare
    title: titolo del grafico
    xlabel,ylabel: etichette degli assi'''
    
    def findFirstFree(basename,i=1):
        while os.path.exist(basename%(i)):
            i=i+1
        return i

    def dividi(l,nmax,DwMax=45):
        if len(l)<=nmax: return [l]
        w=[Wrecover(s) for s in l]
        dw=dict([(i,v) for i,v in enumerate(w)])
        clusters=generalRoutines.QTclust(dw,DwMax)
        lnew=[]
        for c in clusters:
            lnew.append([l[i] for i in c])
        for i,ll in enumerate(lnew[:]):
            wd=dict([(v,w[l.index(v)]) for v in ll])
            lnew[i]=generalRoutines.sortByVal(wd)[0]
        return lnew

    def creaGruppi(fileName,nmax=5):
        f=open(fileName,"r")
        a=[i.split() for i in f.read().split("\n\n")]
        f.close()
        c=[]
        for i in range(len(a)): c.extend(dividi(a[i],nmax))
        cdic=dict([(i,Wrecover(cc[0])) for i,cc in enumerate(c)])
        csort=[c[i] for i in generalRoutines.sortByVal(cdic)[0]]
        return csort
    
    termDic={"png":("png",".png"),"ps":("postscript color",".eps"),"bmp":("bitmap",".bmp"),"":("","")}    
    termExt=termDic[term]
    daPlottare= creaGruppi (groupFile)
    pltName=pltFilePrefix+os.path.splitext(groupFile)[0]+".plt"
    p=open(os.path.join(resdir,pltName),"w")
    p.write("unset label \n")
    p.write("set xlabel 'Energy (keV)'\n")
    p.write ("set yrange %s\n"%yrange)
    p.write ("set xrange %s\n"%xrange)
    p.write ("set ylabel '%s'\n" %(ylabel))
    p.write ("set title '%s'\n" %(title))
    
    outBaseName=outPrefix+os.path.splitext(groupFile)[0]+"%03i"+termExt[1]
    i=1      
    for gr in daPlottare:
        listafiles=[dataFileFrmt%ii for ii in gr if os.path.exists(os.path.join(resdir,dataFileFrmt%ii))]
        if (len(listafiles)==0):continue
        if not overwrite and os.path.exists(pltName):
            i=findFirstFree(os.path.join(resdir,outBaseName),i)
            return
        outName=outBaseName%(i)        

        p.write ("plot ")
        for f in listafiles[0:-1]:
            t=f+" W=%.2f kg"%Wrecover(os.path.splitext(f)[0])
            graspFile=list(os.path.splitext(f))
            graspFile.insert(1,"_grasp")                              
            p.write("'%s' %s title '%s' w l lw 2,\\\n"%("".join(graspFile),using,t))
        f=listafiles[-1]
        t=f+" W=%.2f kg"%Wrecover(os.path.splitext(f)[0])
        graspFile=list(os.path.splitext(f))
        graspFile.insert(1,"_grasp")             
        p.write("'%s' %s title '%s' w l lw 2\n\n"%("".join(graspFile),using,t))
        p.write("pause -1\n")
        if term!="":
            p.write("set term %s\n"%(termExt[0]))
            p.write("set output '%s'\n"%(outName)) 
            p.write("replot\n")
            p.write("set term win\n")
            p.write("set out\n")
        p.write("\n\n")
        p.flush()
        i=i+1
    p.close()


def plottaGrasp(resdir,indice=(0,1),xrange="[0:60]",yrange="[0:450000]",
                comb=[["ML","Ir"],[20,22.5,25,27.5,30],[0.0,0.07,0.15],[60,70]],
                pltFilePrefix="Grasp_",outPrefix="Grasp_",term="png"):
    groupPlot(resdir=resdir,pltFilePrefix=pltFilePrefix,xrange=xrange,yrange=yrange,
              comb=comb,indice=indice,outPrefix=outPrefix,term=term,  
              dataFileFrmt="%s_ff%03i_%im_max%s.txt",using="u 1:3",
              ylabel="Grasp (Aeff x FOV^2) [cm^2 deg^2]")

def plottaGraspPS(resdir,indice=(0,1),xrange="[0:60]",yrange="[0:500000]",
                comb=[["ML","Ir"],[20,22.5,25,27.5,30],[0.0,0.07,0.15],[60,70]],
                pltFilePrefix="GraspsuPlate_",outPrefix="GraspSc_",term="png"):
    groupPlot(resdir=resdir,pltFilePrefix=pltFilePrefix,xrange=xrange,yrange=yrange,
              comb=comb,indice=indice,outPrefix=outPrefix,term=term,  
              dataFileFrmt="%s_ff%03i_%im_max%s.txt",using="u 1:4",
              ylabel="Grasp/PlateScale [(Aeff x Focal x FOV^2)/Rspot]")  

def plottaFOV (resdir,indice=(0,1),xrange="[0:60]",yrange="[0:18]",
              comb=[["ML","Ir"],[20,22.5,25,27.5,30],[0.0,0.07,0.15],[60,70]],
              pltFilePrefix="FOV_",outPrefix="FOV_",term="png"):
    groupPlot(resdir=resdir,indice=indice,pltFilePrefix=pltFilePrefix,xrange=xrange,yrange=yrange,
              comb=comb,outPrefix=outPrefix,term=term,using="u 1:2")  

def plotTemplate(resdir,indice=(0,1),xrange="[*:*]",yrange="[*:*]",
              comb=[["ML","Ir"],[20,22.5,25,27.5,30],[0.0,0.07,0.15],[60,70]],
              pltFilePrefix="",outPrefix="FOV_",term="",  
              combLbl=["Coat","Focal","FF","Dmax"],
              dataFileFrmt="%s_ff%03i_%im_max%s.txt",using="u 1:2",
              xlabel="Energy (keV)",ylabel="Interpolated FOV diameter (arcmin)"):
    
    '''le variabili sono ordinate in modo da far venire dopo quelle fisse per i plot
    di un certo tipo (per es. etichette degli assi o formato dei file da leggere).
    quelle che possono essere modificate (ad esempio i range vengono prima).
    In questo modo e' piu' facile creare routines per plottare file di un certo tipo,
    eliminando le variabili non desiderate che possono generare errori dall'intestazione della routine
    e lasciandole impostate nella chiamata di groupPlot'''

    groupPlot(resdir=resdir,indice=indice,xrange=xrange,yrange=yrange,comb=comb,
              pltFilePrefix=pltFilePrefix,outPrefix=outPrefix,term=term,  
              combLbl=combLbl,dataFileFrmt=dataFileFrmt,using=using,
              xlabel=xlabel,ylabel=ylabel)

if __name__=="__main__":
    os.chdir(r"C:\work2\traie6\simx2008_ufficiale")
    resdir="fov_Energy"
    if not os.path.exists(resdir): os.mkdir(resdir)
    #calcolaFOV_Grasp(resdir)
    plottaGrasp(resdir,(0,2,3))
    plottaFOV(resdir,(0,2,3))
    plottaGraspPS(resdir,(0,2,3))
    
    plottaFOV(resdir,(0,1))
    plottaGrasp(resdir,(0,1))
    plottaGraspPS(resdir,(0,1))
    cmpPlot(resdir,overwrite=1,xrange="[0:60]",yrange="[0:18]",pltFilePrefix="ps",
            outPrefix=os.path.join("lyx","img","FOV_"),
            term="ps",title="Field of View Diameter")
    cmpPlot(resdir,overwrite=1,xrange="[0:60]",yrange="[0:450000]",pltFilePrefix="psGrasp_",
            outPrefix=os.path.join("lyx","img","Grasp_"),term="ps",
            using="u 1:3",ylabel="Grasp (Aeff x FOV^2)",title="GRASP")
    cmpPlot(resdir,overwrite=1,xrange="[20:60]",yrange="[0:140000]",pltFilePrefix="psGrasp2060_",
            outPrefix=os.path.join("lyx","img","Grasp2060_"),term="ps",
            using="u 1:3",ylabel="Grasp (Aeff x FOV^2)",title="GRASP")
    cmpPlot(resdir,overwrite=1,xrange="[0:60]",yrange="[0:500000]",pltFilePrefix="psGraspPs_",
            outPrefix=os.path.join("lyx","img","GraspPs_"),term="ps",
            using="u 1:4",ylabel="Grasp/PlateScale [(Aeff x Focal x FOV^2)/Rspot]",
            title="GRASP / PlateScale")
    cmpPlot(resdir,overwrite=1,xrange="[20:60]",yrange="[0:500000]",pltFilePrefix="psGraspPs2060_",
            outPrefix=os.path.join("lyx","img","GraspPs2060_"),term="ps",
            using="u 1:4",ylabel="Grasp/PlateScale [(Aeff x Focal x FOV^2)/Rspot]",
            title="GRASP / PlateScale")
    #plottaFOV(resdir,(0,1))
    '''cmpPlot(resdir,overwrite=1,pltFilePrefix="test_",xrange="[0:60]",yrange="[0:18]")'''
    '''cmpPlot(resdir,overwrite=1,xrange="[0:60]",yrange="[0:18]",pltFilePrefix="best_ps",
            outPrefix=os.path.join("lyx","img","best_FOV_"),
            term="ps",title="Field of View Diameter")
    cmpPlot(resdir,overwrite=1,xrange="[0:60]",yrange="[0:450000]",pltFilePrefix="best_psGrasp_",
            outPrefix=os.path.join("lyx","img","best_Grasp_"),term="ps",
            using="u 1:3",ylabel="Grasp (Aeff x FOV^2)",title="GRASP")
    cmpPlot(resdir,overwrite=1,xrange="[20:60]",yrange="[0:140000]",pltFilePrefix="best_psGrasp2060_",
            outPrefix=os.path.join("lyx","img","best_Grasp2060_"),term="ps",
            using="u 1:3",ylabel="Grasp (Aeff x FOV^2)",title="GRASP")
    cmpPlot(resdir,overwrite=1,xrange="[0:60]",yrange="[0:500000]",pltFilePrefix="best_psGraspPs_",
            outPrefix=os.path.join("lyx","img","best_GraspPs_"),term="ps",
            using="u 1:4",ylabel="Grasp/PlateScale [(Aeff x Focal x FOV^2)/Rspot]",
            title="GRASP / PlateScale")
    cmpPlot(resdir,overwrite=1,xrange="[20:60]",yrange="[0:500000]",pltFilePrefix="best_psGraspPs2060_",
            outPrefix=os.path.join("lyx","img","best_GraspPs2060_"),term="ps",
            using="u 1:4",ylabel="Grasp / PlateScale [(Aeff x Focal x FOV^2)/Rspot]",
            title="GRASP/PlateScale")'''

    cmpPlot(resdir,overwrite=1,pltFilePrefix="FOV_",xrange="[0:60]",yrange="[0:18]",
            term="png",title="Field of View Diameter")
    cmpPlot(resdir,overwrite=1,xrange="[0:60]",yrange="[0:18]",pltFilePrefix="ps",
            outPrefix=os.path.join("lyx","img","FOV_"),
            term="ps",title="Field of View Diameter")
    

if __name__=="x__main__":
    os.chdir(r"D:\2\erosita")
    resdir="grasp"
    cmpPlot(resdir,overwrite=1,groupFile="scelti.txt",xrange="[0:60]",yrange="[0:18]",pltFilePrefix="FOVpng_",
            outPrefix="FOV_",
            term="png",title="Field of View Diameter")
    cmpPlot(resdir,overwrite=1,groupFile="scelti.txt",xrange="[0:60]",yrange="[0:450000]",pltFilePrefix="Grasppng_",
        outPrefix="Grasp_",term="png",
        using="u 1:3",ylabel="Grasp (Aeff x FOV^2)",title="GRASP")
    cmpPlot(resdir,overwrite=1,groupFile="scelti.txt",xrange="[0:60]",yrange="[0:500000]",
        pltFilePrefix="GraspPspng_",
        outPrefix="GraspPs_",term="png",
        using="u 1:4",ylabel="Grasp/PlateScale [(Aeff x Focal x FOV^2)/Rspot]",
        title="GRASP / PlateScale")

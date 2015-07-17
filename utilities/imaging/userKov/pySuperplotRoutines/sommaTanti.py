
import os
import string

'''prende i nomi di directory contenuti in scelti.txt e per ogni directory effettua la
somma usando il programma sommaSelectSet(sss.exe), con le impostazioni contenute in schema_imp_sommaSelSet.txt.
di cui cambia dirorigin e dirrisult, in base alla configurazione estratta dal nome della cartella
crea un file di plot che carica tutti i file di plot delle somme che vengono generati nelle varie
directory. Per ogni somma genera anche un plot di confronto con il risultato del ray-tracing.
i risultati del ray tracing vengono presi in modo poco standard, probabilmente questo pezzo e' da modificare'''

def sss(conf,p):
    ''' sistema il file delle impostazioni e dei parametri per la configurazione
    estratta dalla stringa conf '''

    def ifTr(conf):
        if conf[0:2]=="Ir":
            return "ml"+conf[2:]
        else:
            return conf
        
    s=conf.split("_")
    coating,ff,focal,dmax=s
    focal=focal[0:2]
    dmax=dmax[-3:]
    #print coating,ff,focal,dmax
    replace={}
    replace["dirrisult"]="'"+str(conf)+"_sum'"
    replace["dirorigin"]="'"+str(ifTr(conf))+"'"
    toRead={}
    toRead["dirorigin"]=''
    toRead["workdir"]=''

    #sistema il file imp_sommaSelSet.txt   
    orig=open("schema_imp_sommaSelSet.txt","r")
    if os.path.exists("imp_sommaSelSet.txt"):
        os.remove("imp_sommaSelSet.txt")
    dest=open("imp_sommaSelSet.txt","w")
    for linea in orig.readlines():
        possKey=string.replace(string.lower(linea.split("=")[0])," ","")
        if replace.has_key(possKey):
            dest.write(linea.split("=")[0]+" = "+replace[possKey]+"\n")
        else:
            dest.write(linea)
        if toRead.has_key(possKey):
            toRead[possKey]=string.replace(string.replace(string.lower(linea.split("=")[1])," ",""),"\n","")
            toRead[possKey]=string.replace(toRead[possKey],"'","")
    if(string.lower(coating)=="ml" or string.lower(coating)=="m2"):
        dest.write("1 100 0 1  shell              ")
    elif(string.lower(coating)=="ir"):
        dest.write("1 100 1 1  shell              ")
    else:
        print "valore di coating ",string.lower(coating)," non riconosciuto"
        Beep(440,1000)
        return
    orig.close()
    dest.close()
    #scrive il file di plot
    p.write("cd '"+toRead["workdir"]+ "\\" + string.replace(replace["dirrisult"],"'","")+"'\n")
    p.write("load 'sommaPNG.plt'\n")
    p.write("cd '..\\..'\n")
    
    p.write("unset logscale y\n")

def plotCmp(conf,ofile,pathList,workdir):

    def file2cmp(conf,pathList):
        '''prova una serie di directory/nomefile (in pathList), che potrebbero contenere risultati
        precedenti per un confronto/controllo, se non trova i file non li fa plottare.'''
        pp=[]            
        for p in pathList:
            #print p[0]+conf+p[1],os.path.exists(p[0]+conf+p[1])
            if (os.path.exists(p[0]+conf+p[1])):pp.append(p[0]+conf+p[1])
        if len (pp)==0:
            return ""
        else:
            return pp

    cmpFile=file2cmp(conf,pathList)
    if cmpFile<>"":
        ofile.write("set term win\n")
        ofile.write("set out\n\n")
        ofile.write("plot '"+workdir+"\\"+str(conf)+"_sum")
        ofile.write("\\sommaTOT.txt' u 1:2 w l")
        for plf in cmpFile:
            ofile.write(",\\\n'"+plf)
            ofile.write("' u 1:2 i 0 w l")
        ofile.write("\nset term png\n")
        ofile.write("set out '"+conf+".png'\n")
        ofile.write("replot\n")
        ofile.write("set term win\n")
        ofile.write("set out\n\n")

    
dirOffAxis="G:\vince\prog\offaxis\OA\simX2006"
f=open("scelti.txt","r")
p=open("plottaSum.plt","w")
pc=open("plottaCmp.plt","w")

for conf in f.readlines():
    conf=conf[:-1]
    sss(conf,p)
    os.system("sss.exe")
    #creaFile di plot per confronto con vecchi risultati
    pathList=[("simX2006onAxis_old\\","_sum\\sommaTOT.txt"),
                                ("G:\\vince\\prog\\offaxis\\OA\\simX2006\\","\\aree.txt")]
    plotCmp(conf,pc,pathList,"simX2006onAxis")
    
f.close()
p.close()
pc.close()
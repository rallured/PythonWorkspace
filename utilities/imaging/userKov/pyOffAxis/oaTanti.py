
import os
import string

'''ricavato da simulaTanti per superplot.
per la serie di directory <dir> contenute nel file scelti.txt lancia il
programma per la simulazione delle aree efficaci (superplot) -> compilato in su(_ml | _rt).exe
per multilayer o riflessione totale
le impostazioni sono prese da schema_datiplot.txt, in cui vengono rimpiazzati
focale, dirrisult e diamfile, ricavati da <dir>.

Il file dei diametri viene generato in af_files con il nome della directory, i valori
per i diametri vengono estratti dal file transryatrace.dat contenuto in ogni sottodirectory
<dir> dei risultati raytracing (contenuti nella directory <dirOffAxis> impostata in questo prog).

la directory dei dati di partenza e' impostata in schema_datiplot.txt nella variabile dirorigin
in questa directory viene generato e poi cancellato per ogni <dir> il file parametri.txt.
'''


def creaimpo(conf,doa):
    ''' sistema i file delle impostazioni e dei parametri per la configurazione
    estratta dalla stringa conf, doa e' la directory che contiene i risultati per
    le diverse configurazioni'''
    
    s=conf.split("_")
    coating,ff,focal,dmax=s
    focal=focal[0:2]
    dmax=dmax[-3:]
    #print coating,ff,focal,dmax
    replace={}
    replace["dirrisult"]="'"+str(conf)+"'"
    replace["f_length"]=float(focal)*100
    replace["diamfile"]="'"+str(string.replace(conf,"_",""))+".txt'"
    if coating=="m2":
        replace["mat1odd"]="platinum.dat"
        replace["mat2even"]="carbon.dat"
        replace["matest1st"]="carbon.dat"
        replace["matest2nd"]="platinum.dat"
    toRead={}
    toRead["dirorigin"]=''
    toRead["workdir"]=''
    toRead["ngruppi"]=''
    
    #sistema il file 3_datiplot    
    orig=open("schema_datiplot.txt","r")
    if os.path.exists("3_datiplot.txt"):
        os.remove("3_datiplot.txt")
    dest=open("3_datiplot.txt","w")
    for linea in orig.readlines():
        #print "*",linea,string.lower(linea.split("=")[0]),replace.has_key(string.lower(linea.split("=")[0])),"*"
        possKey=string.replace(string.lower(linea.split("=")[0])," ","")
        if replace.has_key(possKey):
            dest.write(linea.split("=")[0]+" = "+replace[possKey]+"\n")
            #print string.lower(linea.split("=")[0])
        else:
            dest.write(linea)
        if toRead.has_key(possKey):
            toRead[possKey]=string.replace(string.replace(string.lower(linea.split("=")[1])," ",""),"\n","")
            toRead[possKey]=string.replace(toRead[possKey],"'","")
    orig.close()
    dest.close()

    #sistema il file parametri.txt

    diror=os.path.join(toRead["workdir"],toRead["dirorigin"])
    dirlist=open(os.path.join(diror,"dirlist.dat"),"r")
    for dirfom in dirlist.readlines():
        if dirfom <>"\n":
            print dirfom
            dirfom = string.replace(dirfom,"\n","")
            #print dirfom
            #os.rename(os.path.join(diror,dirfom,"parametri.txt"),os.path.join(diror,dirfom,"old_parametri.txt"))
            orig=open(os.path.join("schema_parametri.txt"),"r")
            if os.path.exists(os.path.join(diror,dirfom,"parametri.txt")):
                os.remove(os.path.join(diror,dirfom,"parametri.txt"))
            dest=open(os.path.join(diror,dirfom,"parametri.txt"),"w")
            for e in orig.readlines():
                possKey=string.replace(string.lower(e.split("=")[0])," ","")
                #print possKey
                if replace.has_key(possKey):
                    dest.write(e.split("=")[0]+" = "+str(replace[possKey])+"\n")
                    #print string.lower(linea.split("=")[0])
                    #print e.split("=")[0]+" = "+replace[possKey]+"\n"
                else:
                    dest.write(e)
                    #dest.write(linea)
            orig.close()
            dest.close()

            #crea il file dei diametri in af_files
            if coating=="m2":
                confor=conf[0:1]+"l"+conf[2:]
            else:
                confor=conf
            origin=open(os.path.join(doa,confor,"transRayTrace.dat"),"r")
            l=origin.readlines()
            for i,v in enumerate(l):
                if string.find(v,"radius at the intermediate pupil [mm]")<>-1:
                    break
            diam=l[i+100:i:-1]
            diam =[float(d)*2 for d in diam]
            df=open(os.path.join(toRead["workdir"],"af_files",string.replace(replace["diamfile"],"'","")),"w")
            #print df
            for d in diam:
                df.write(str(d)+"\n")
                #print d
            df.close()
            origin.close()   
    dirlist.close()
    
if __name__=="main":    
    dirOffAxis="G:\\vince\\prog\\offaxis\\OA\\simX2006"
    f=open("scelti.txt","r")    
    for conf in f.readlines():
        conf=conf[:-1]
        creaimpo(conf,dirOffAxis)
        os.system("su.exe")
        #crea plot
    f.close()    
def load (file,n,nSkip=0,index=-1):
        ''' carica da un file l'n-sima (in base 1) colonna, restituisce lista
        di liste di valori letti (una per ogni blocco di dati)
        nSkip indica quante righe saltare, se =-1 usa prima riga come intestazione delle colonne
        se =-2, salta dati non numerici,
        se =0 da' errore in caso di dati non numerici, se carattere lo usa come
        carattere di commento'''
        
        f=open(file,"r")
        lin=f.readlines()
        i=0
        area = [[]]
        flag = False
        remChar=""
        start=0
        if isinstance(nSkip, int):
            if nSkip>0:
                start=nSkip
            elif nSkip==-1:
                start=1
                self.label=lin[0].split()[n]
            elif nSkip==-2:
                def isNum(a):
                    res=True
                    for i in a:                        
                        try: float(i)
                        except: res=False
                    return res
                lin=[l for l in lin if isNum(l)]
        elif isinstance(nSkip, str):
            remChar=nSkip
        else:
            print "valore non riconosciuto per il parametro nskip"
            print "in load: ", nSkip

        #print area
        for l in lin[start:]:
            #print i,area[i]
            if l.split()==[]: #se c'e' riga vuota aggiunge una lista da riempire
                if flag:    #righe vuote consecutive contano come una sola
                    i=i+1
                    flag=False
                else:
                    flag=True
                    area.append([]) 
            else:
                if l.strip()[0] != remChar:area[i].append(float(l.split()[n-1]))
        f.close()
        if area[-1] ==[]:area.pop()
        t=area.pop()
        return t
    
class CplotInfo:
    '''classe contenente parametri per il plot. a ogni Ctrack da plottare si fa corrispondere
    un oggetto CplotInfo'''
    def __init__(self,label="Track"):
        self.type = "w l"
        self.using = [1,2]
        
class Ctrack:
    '''traccia contenente dati'''
    def __init__(self,data=[],label="Track"):
        self.data=data
        self.label=label
        
    def load(self,file,n,nSkip=0,index=-1):
        self.data=load(file,n,nSkip=0,index=-1)
    
    def __repr__(self):
        if len(self.data)!=0:s2="\n[%s,.. ..,%s]"%(self.data[0],self.data[-1])
        else: s2="[]"
        return "%s(%s el):%s"%(self.label,len(self.data),s2)        
    
    def __mul__(self,n):
        '''non e' tanto corretto che ritorni un array invece che un istanza di
        Ctrack, comunque e' sperimentale'''
        if isinstance(n,float):
            return [d * n for d in self.data]
        
class Cdata(list):
    '''gruppo di Ctrack contenenti dati, si puo' usare per generare file di dati e file di plot'''
    def __init__(self):
        self.xrange="[*:*]"
        self.yrange="[*:*]"
        self.xlabel="X"
        self.ylabel="Y"
        self.tracks=[]
        list.__init__(self)
        
    def __repr__(self):
        return "\n".join(["["]+[str(i).replace("\n"," ") for i in self]+["]"])


if  __name__=="__main__":
    test=Ctrack()
    test.load(
        "G:\\vince\\prog\\f_2005\\ufficiali\\prog_superplotGen1\\script\\HXMT_8_integrali_test\\area06_testNshell_plt.txt"
        ,2,-1)
    test2=Ctrack()
    test2.load(
        "G:\\vince\\prog\\f_2005\\ufficiali\\prog_superplotGen1\\script\\HXMT_8_integrali_test\\area06_testNshell_plt.txt"
        ,3,-1)
    t=Cdata()
    t.append(test)
    t.append(test2)
    print t
    d=Ctrack().load(
        "G:\\vince\\prog\\f_2005\\ufficiali\\prog_superplotGen1\\script\\HXMT_8_integrali_test\\area06_testNshell_plt.txt"
        ,4,-1)
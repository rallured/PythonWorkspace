from powLaw import pl4,multilayer
import random
from numpy import *

# schema proprieta' get e set ##############
class C(object):
    def __init__(self): self.__x = None
    def getx(self): return self.__x
    def setx(self, value): self.__x = value
    def delx(self): del self.__x
    x = property(getx, setx, delx, "I'm the 'x' property.")
# ###########################################    

class point:
    '''un punto nello spazio dei parametri. viene inizializzato con 
    una lista, la cui lunghezza determina ndim. se un elemento della
    lista e' valore singolo, lo assume, se e' lista genera casuale tra
    i primi due elementi. ha anche proprieta' flag a 0 di default'''
    def __init__(self,ran,flag=0):
        self.ndim=len(ran)
        self.flag=flag
        self.r=[]
        #print ran
        for i in ran:
            #print i,type(i),type(i)=="<type 'list'>"
            if isinstance (i,list):
                self.r.append(random.uniform(i[0],i[1]))
            elif not isinstance (i,list):
                self.r.append(i)
            else: errore
    def __getitem__(self,ind):
        return self.r[ind]
    def __setitem__(self,ind,val):
        self.r[ind]=val
    
    def writeOnFile(self,file):
        '''scrive il punto sul file dato'''
        for i in self.r:
            file.write(str(i)+"\t")
        file.write(str(self.flag))
        file.write("\n")      

    def __repr__(self):
        '''print per test'''
        s=""
        for i in self.r:
            s=s+ str(i)+"\t"
        s=s+ "flag: "+str(self.flag)
        return s

class pop:
    def __init__(self,npop):
        #crea un set di parametri
        self.npop=npop
        self.__pop=[]
        for i in range(npop):
            self.__pop.append(point(ran))
    def __getitem__(self,ind):
        return self.__pop[ind]
    def __setitem__(self,ind,p):
        if isinstance (p,list):
            p=point(p)
            print "*"
        if isinstance (p,point):            
            if p.ndim==self.__pop[-1].ndim:                
                self.__pop[ind]=copy.deepcopy(p)
            else:
                print "fornito un punto con %i dimensioni\
            invece che con %i dimensioni" %(p.ndim,self.__pop[-1].ndim)   
        else:
            print "il valore fornito per un elemento \
            della \npopolazione deve essere di tipo point o list!"
    def __repr__(self):
        s=""
        for i,p in enumerate(self.__pop):s=s+"%03i-\t" %(i)+str(p)
        return s
    def popWrite(self,file):
        for p in self.__pop:
            p.writeOnFile(popFile)
    def popRead(self,file):
        self.__pop=[]
        for p in file.readlines():
            pp=[float(i) for i in p.split("\t")[:-1]]
            f=p.split("\t")[-1]
            #print "flag=",f,pp
            self.__pop.append(point(pp,f))
            
if __name__=="__main__":
    ran=[[10,200],[-1,100],[0.01,1],[0.1,0.9]]
    npop=5
    p=pop(npop)

    #lo scrive su file
    popFile=open("population.txt","w")
    p.popWrite(popFile)
    popFile.close()

    #legge i valori da file, per ognuno:
    popFile=open("population.txt","r")
    #for p in popFile.readlines():
        #li trasforma in spessori
    p.popRead(popFile)

        #fitta i valori
        #se flaggato crea file e scrive spessori fittati e reali 
        #scrive i valori di partenza e fittato e fom
    popFile.close()
    ml2=multilayer(100,p[2])
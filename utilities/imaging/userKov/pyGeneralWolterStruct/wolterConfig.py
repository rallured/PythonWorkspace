import os
import string
'''manipolazione di configurazioni per geometria wolter'''

# accessori
class Callable:
    ''' da internet: http://aspn.activestate.com/ASPN/Cookbook/Python/Recipe/52304 (alex martelli)
    This is easy to solve with a simple tiny wrapper:
    serve come wrapper per creare metodo statico (si puo' chiamare senza istanziare classe)
    per compatibilita', permette chiamate come vecchie routine'''
    def __init__(self, anycallable):
        self.__call__ = anycallable
#fine accessori
        
class config:
    '''classe che gestisce le informazioni relative a una data configurazione geometrica
    wolter'''
    coating={"ir":"Ir monolayer","ml":"W/Si multilayer","m2":"Pt/C multilayer"}
    ff={"ff000":"0.0","ff00":"0.0","ff007":"0.07","ff015":"0.15"}
    diam={"max60":"60","max70":"70"}
    focale={"20m":"20","22m":"22.5","25m":"25","27m":"27.5","30m":"30"}
    def extractTitle(st):
        '''dalla stringa con il nome della directory estrae le informazioni
        sulla configurazione e restituisce una lista con coating, filling factor,
        focale, diametro massimo. (questa e' la vecchia routine)'''
        l=[string.lower(w) for w in st.split("_")]
        return [config.coating[l[0]],config.ff[l[1]],config.focale[l[2]],config.diam[l[3]]]
    extractTitle = Callable(extractTitle) #usa wrapper precedentemente definito

    def __init__(self,configString):
        self.config=configString

    def __str__(self):
        s=[]
        s.append("configurazione: %s" %self.config)
        s.append("coating: %s" %self.config)
        s.append("configurazione: %s" %self.config)
        s.append("configurazione: %s" %self.config)
        



        
def extractTitle(st):
    '''per compatibilita' con vecchie versioni: dalla stringa con il nome
    della directory estrae le informazioni sulla configurazione e restituisce
    una lista con coating, filling factor, focale, diametro massimo'''
    return config.extractTitle(st)

if __name__=="__main__":
    #esempio nome convenzionale: m2_ff015_20m_max70
    print "effettuo un test del modulo per le configurazioni:"
    print "prova sia la chiamata al metodo unbound della classe"
    print "che quella alla routine extractTitle per compatibilita'"
    print "output nel file testout.txt\n"
    print "Analizzo le configurazioni contenute nel file testConfList.txt:"
    
    testStrings=open("testConfList.txt","r").read().split()
    for s in testStrings: print s
    res=[]
    for rout in [extractTitle, config.extractTitle]:        
        res.append( "\nTEST con la routine %s\n" %rout)
        for st in testStrings:
            res.append( "----------------")
            res.append( "configurazione: %s" %s)
            res.append( "coating: %s\nfillingfactor: %s\ndiamMax: %s\nfocal length: %s" %tuple(rout(st)))
            res.append( "----------------\n")

    res="\n".join(res)
    print res
    open("testout.txt","w").write(res)
    
    
        


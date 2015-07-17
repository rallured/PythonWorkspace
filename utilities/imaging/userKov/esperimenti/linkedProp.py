'''test per i riferimenti e valori linkati:
vorrei costruire un oggetto 'telescopio', che contenga oggetti 'shell'
con varie proprieta' (es. diametro). vorrei poter accedere alle shell
come elementi numerati del telescopio, e avere una proprieta' diametri
di telescopio, contenente una lista dei valori delle proprieta' diametro
delle shell. in modo che se cambio un valore in telescopio.shell cambia
anche in telescopio.diametro. se assegno uno scalare a telescopio.diametro
lo riferisce a tutte le shell. si puo' fare? (con puntatori si potrebbe fare)'''

class shell(object): #questa e' una shell, che possiede la sola proprieta' d
    def __init__(self,diam):
        self.d=float(diam)
    
class telescope(list):
    '''questo e' un telescopio, cioe' una lista di oggetti shell. ha
    una proprieta' d: lista che contiene i diametri delle shell. voglio
    verificare la relazione tra telescopio[n].d e telescopio.d[n]'''

    '''def __setitem__(self,n,val):
        self[n]=val
    def __getitem__(self,n):
        return self[n]'''

    def __init__ (self,n):
        list.__init__(self)
        self.d=[]
        for i in range(n):
            self.append(shell(i*1.3)) 
            self.d.append(self[-1].d)
        '''questo da' gli stessi risultati
        for i in range(n):
            s=shell(i*1.3) #valori di fantasia ai diametri 
            self.append(s) 
            self.d.append(s.d)'''

if __name__=="__main__":
    def printCurrentProp(t):
        print "il secondo elemento t[1]=",t[1]
        print "la proprieta' t.d di t vale:",t.d
        print "la proprieta' d della seconda shell vale:", t[1].d
        print "--------------\n"
        
    t=telescope(5)
    print "\n\n--------"
    print "creo telescopio con 5 shell: t=telescope(5)"
    print "l'oggetto t e' :"
    print t
    printCurrentProp(t)
    print "ora cambio la proprieta' d di t con t.d[1]=7"
    t.d[1]=7
    print "il diametro della shell t[1].d vale",t[1].d
    print "NON SONO COLLEGATI!!"

    print "\n---- il contrario:"
    print "riinizializzo: t=telescope(5)"
    print "l'oggetto t e' :"
    print t
    print "NON RIINIZIALIZZA!!"
    printCurrentProp(t)
    print "ora cambio la proprieta' d della shell con t[1].d=5"
    t[1].d=5
    print "t.d vale",t.d
    print "ANCHE COSI' NON SONO COLLEGATI!!"
    
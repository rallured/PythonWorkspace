'''serie di test, composta da due test principali chiamati da main, su assegnazione
di valori a liste o loro elementi in cicli for e funzioni.'''

def addOne_RetList(list):
    '''aggiunge 1 a ogni elemento della lista, iterando sulla lista stessa
    list=[l+1 for l in list]. restituisce la lista stessa. se non gli si
    fa restituire valore non cambia niente. iterando su una copia list[:]
    non cambia niente.'''
    list=[l+1 for l in list]
    print "lista nella funzione: %s" %list
    return list

def replace2ndEl(list):
    '''rimpiazza il secondo elemento assegnando un nuovo valore a list[1].'''
    list[1]="elemento aggiunto"
    print "lista nella funzione: %s" %list
    return list

def replaceAllList(list):
    '''assegna un nuovo valore (di tipo stringa) alla lista.'''
    list=" Questo e' un nuovo valore stringa!!"
    print "lista nella funzione: %s" %list
    return list

def replaceConFor(list):
    '''aggiunge 1 a tutti gli elementi della lista usando un ciclo for con enumerate:
    for i,aa in enumerate (list):
     	list[i]=aa*2'''
    for i,aa in enumerate (list):
     	list[i]=aa*2
    print "lista nella funzione: %s" %list
    return list


    
def testPassListToFunk():
    print "\n\n----------------------------------"
    print "test su passaggio di liste ad una funzione. provo diverse",
    print "routines che agiscono in modo diverso su una lista.."
    print "-> vedi anche esempio su cicli for..\n"
    print "----------------------------------"
    for f in [addOne_RetList,replace2ndEl,replaceAllList,replaceConFor]:
        l=[1,2,3,4]
        print "provo la routine: %s" %f.__name__
        print f.__doc__
        print "lista input: %s"%l
        x=f(l)
        print "valore restituito: %s" %x
        print "lista dopo funzione: %s\n\n"%l
        raw_input()
        print "\n----------------------------------"
    print "\nCOMMENTO: le liste vengono modificate anche all'esterno se nella routine "
    print "si cambia uno dei valori. Se si modifica tutta la lista, le modifiche non si "
    print "rispecchiano alla fine."
    print "il comportamento e' quanto meno strano, cosa c'e' sotto?"
    print "(come conseguenza, se modifico la lista nella sub con un ciclo for"
    print "le modifiche rimangono!!"
#---------------------------------------------------------------------
def testListInCicloFor():
    '''testa il comportamento delle liste se si tenta di
    modificarne il valore nei cicli for '''

    def sep(a):
        '''stampa formattato e con separatore finale, il valore di a. mette pausa "premi un tasto"'''
        print "--> a=%s\n"%a 
        raw_input()
        print "----------------------------------"

    print "\n\n"+2*"----------------------------------\n"+testListInCicloFor.__doc__
    print "assegnazione:"
    a=range(5)    
    sep(a)
    
    print "la lista raddoppia se la modifico con:"
    print "a=[aa*2 for aa in a]"
    a=[aa*2 for aa in a]
    sep(a)
    print "invece non cambia se viene fatto dentro un ciclo for:"
    print "for aa in a:"
    print "\taa=aa*2"
    for aa in a:
        aa=aa*2
    sep(a)





if __name__=="__main__":
    testPassListToFunk()
    #testListInCicloFor()
def routine1():
    b='BBB'
    print "rout1 definisce b='BBB'"
    print 'print a, b:'
    print a,b

def routine2():
    print 'a e b non sono definiti in rout2'
    print 'print a, b:'
    print a,b
    print '(ha usato i valori del main)'
    
def routine3():
    print 'stampa una variabile non definita nel main'
    print '(anche se fosse definita in altra routine)'
    print 'print c'
    print "da' errore:"
    print c
    
a='A'
b='B'
print "main definisce a='A', b= 'B'"
print 'e chiama routine1()\n'
print '/-----------------------\\'
routine1()
print '\\----------------------/'
print 'tornato al main che chiama routine2()\n'
print '/-----------------------\\'
routine2()
print '\\----------------------/'
print 'tornato al main che chiama routine3()\n'
print '/-----------------------\\'
routine3()
print '\\----------------------/'


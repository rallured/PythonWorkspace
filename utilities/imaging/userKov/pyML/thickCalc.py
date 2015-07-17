def pl2th(a,b,c,N):
    try:
        iter(N)
        n1=N[0]
        n2=N[1]
    except:
        n1=1
        n2=N

    t=[a/((i+b)**c) for i in range (n1,n2+1)]
    return t

def th2pl(D,N,c):
    try:
        iter(N)
        n1=N[0]
        n2=N[1]
    except:
        n1=1
        n2=N
    d1,d2=D[0],D[1]

    x=(d1/d2)**(1/c)
    b=(n2-n1*x)/(x-1)
    return b

    

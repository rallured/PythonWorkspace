from pylab import *
import os

def readMis(filename):
    '''read a measure file in standard format, return
    x, y, error on y (from the first 3 columns) and time in seconds'''
    lines=open(filename,"r").read().split("\n")
    mf=0
    x=[]
    y=[]
    erry=[]
    time=[]
    for l in lines:
        if l:
            ll=l.split()
            x.append(float(ll[0]))
            y.append(float(ll[1]))
            erry.append(float(ll[2]))
            t=map(float,ll[3].split(":"))
            time.append(t[0]*3600+t[1]*60+t[2])
    return x,y,erry,time
    
misFile='14_ML338C_3600sec.dat'
normFile="_DB".join(os.path.splitext(misFile))

x1,y1,ye1,t1=map (array,readMis(misFile))
x2,y2,ye2,t2=map (array,readMis(normFile))
t1=t1-t1[0]
t2=t2-t2[0]
if x1!=x2:
    print "attenzione, le x sono diverse..."
    exit

b1start,b1stop,b2start,b2stop=map(float,raw_input("insert 4 beams").split())
#interpola il valore del fascio basandosi sul tempo
t1=t1-t1[0]
t2=t2-t1[0]
b1m=(b1stop-b1start)/t1[-1]
b1=b1m*t1+b1start
b2m=(b2stop-b2start)/(t2[-1]-t2[0])
b2=b2m*t2+b2start
reflex=(y1/y2)*(b2/b1) #normalizzazione con i valori corretti per il tempo

#output su file
outFile=open("_out".join(os.path.splitext(misFile)),"w")
for x,r in zip(x1,reflex):
    outFile.write("%s\t%s\n"%(x,r))
outFile.close()

plot (x1,reflex)
show()



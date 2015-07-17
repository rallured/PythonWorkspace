'''corregge i pesi calcolati (da densita' 8 a densita' 8.8)
per tutte le
directory elencate nel file tutti.txt
file modificati:
imp_OffAxis.txt  -->
pesi.txt  -->
'''

import os
import string
if __name__=="__main__":
    elenco=open("tutti.txt","r")
    dl=elenco.read().split()
    elenco.close()
    count=0
    count2=0
    for d in dl:
        flag=0
        d=d[:-1]
        c_in=open(os.path.join(d,"imp_OffAxis.txt"),"r")
        c=c_in.readlines()
        c_in.close()

        c_out=open(os.path.join(d,"imp_OffAxis.txt"),"w")
        for cc in c:
            if cc.split("=")[0]=="WallDensity":
                if float(cc.split("=")[1])==8:
                    print "sostituisco prontamente in ",d
                    count=count+1
                    cc=string.replace(cc,"8","8.8")
                    flag=1
                else:
                    print "cazzo! in %s WallDensity=%s" %(d,str(cc.split("=")[1]))
            c_out.write(cc)
        c_out.close()

        if flag==1:
            count2=count2+1
            c_in=open(os.path.join(d,"pesi.txt"),"r")
            c=c_in.readlines()
            c_in.close()

            c_out=open(os.path.join(d,"pesi.txt"),"w")
            
            for i in range(len(c)):
                #print i
                if i==0:
                    cc=" "+c[i].split()[0]+" "+c[i].split()[1]+"      "+str(float(c[i].split()[2])*1.1)+"\n"
                elif i>3:
                    cc=str(float(c[i])*1.1)+"\n"
                    #print c[i],cc
                else:
                    cc=c[i]
                #print i, cc
                c_out.write(cc)
            c_out.close()            
                

    print "corretti ", count, "file delle impostazioni e "
    print  count2, "file dei pesi"
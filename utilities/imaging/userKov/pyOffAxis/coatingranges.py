import os
import extractInfo
from generaPlot2 import extractTitle

fileLista="lista.txt"
os.chdir (r"C:\work\copiaOA\batch\batch4\simx2007")
tabFile=open("tab_coating.txt","w")
dicShCoat={0.0036:"Pt/C",0.00245:"Ni/C",0.1:"W/Si"}
aKeys=dicShCoat.keys()
aKeys.sort()
tabFile.write("ff\tfocal\tDmax\t"+"\t".join([dicShCoat[a] for a in aKeys])+"\n")

def shellCoatings():
    ranges={}
    ll=open(fileLista,"r").read().split("\n")
    for ff in ll:
        ang=map(float,extractInfo.shellAngles(ff))
        nshell=len(ang)
        print "contate %s shell" %nshell
        for k,v in dicShCoat.items():
            ranges[v]=[i<=k for i in ang].count(True)
        print ranges
        shellList=[ranges[dicShCoat[a]] for a in aKeys]
        shTot=shellList[0]
        for i in range(1,len(shellList)):
            shellList[i]=shellList[i]-shTot
            shTot=shTot+shellList[i]
        tabFile.write("\t".join(extractTitle(ff)[1:])+"\t"+"\t".join(map(str,shellList))+"\n")

shellCoatings()
tabFile.close()
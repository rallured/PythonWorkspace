import os


nshell=100
outDir="layer_err10perc"
parametri=[115.5,-0.9,0.27,0.35]
#parametri=[105,-0.9,0.27,0.35]
if not os.path.exists(outDir):
    os.mkdir(outDir)
    for i in range(1,nshell+1):
        f=open(os.path.join(outDir,"sceltigruppo%03i.dat"%i),"w")
        for j in parametri[:-1]: f.write("%23.16f\t"%j)
        f.write("%23.16f\t"%parametri[-1])
        f.close()
    
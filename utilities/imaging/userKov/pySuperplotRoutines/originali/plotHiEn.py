''' programma usa e getta, crea plot per le aree in asse
a gruppi '''

import os
import string

lof="scelti.txt"
mainDir="simX2006onAxis"

coat=["m2_","ml_"]
focal=["20m_","22m_","25m_","27m_","30m_"]
ff=["ff000_","ff007_","ff015_"]
dmax=["max60","max70"]

coatName={"m2_":"Pt/C multilayer","ml_":"W/Si multilayer"}
ffstring={"ff000_":"0.0","ff007_":"0.007","ff015_":"0.015"}
dmaxstring={"max60":"60 cm","max70":"70 cm"}

def fileToPlot(c,f,ff,d):
    file=os.path.join(mainDir,c+ff+f+d+"_sum","somma100.txt")
    return file

pf=open("plotHE.plt","w")
pf.write("set xlabel 'Energy (keV)'\n")
pf.write("set ylabel 'Effective area(cm^2)'\n")
for c in coat:
    for f in focal:
        outFile=c+f+"_HE.png"
        plotTitle=coatName[c]+", Focal Length= " + f[:-1]
        filesToPlot=[[fileToPlot(c,f,u,v),v,u] for v in dmax for u in ff ]
        pf.write("set title '"+plotTitle+"'\n")
        pf.write("plot '"+filesToPlot[0][0]+"' u 1:2 title 'Dmax="+dmaxstring[filesToPlot[0][1]]
              +", fil. fac.="+ffstring[filesToPlot[0][2]]+"' w l")
        for l in filesToPlot[1:]:
            v=filesToPlot.index(l)
            pf.write(",\\\n'"+filesToPlot[v][0]+"' u 1:2 title 'Dmax="+dmaxstring[filesToPlot[v][1]]
              +", fil. fac.="+ffstring[filesToPlot[v][2]]+"' w l")
        #print "*",filesToPlot,"*"
        pf.write("\nset terminal png\n")
        pf.write("set output '%s'\n" %outFile)
        pf.write("replot\n")
        pf.write("set term win\n")
        pf.write("set out\n")
        pf.write("\n\n")

pf.close()        
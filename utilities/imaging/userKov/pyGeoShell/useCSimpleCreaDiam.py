def plottaConfronto(nl,comblist):
    '''presa da versione di test, ma non testata in questo contesto. crea file di plot
    per confronto con vecchi risultati'''
    p=open("confrDiam.plt","w")
    for comb in comblist:
        ndf="F"+str(comb[0]).split(".")[0]+"_ff"+str(comb[1]).split(".")[1].zfill(3)+"_D"+str(comb[2])+".txt"
        nl2.FL,nl2.ff,nl2.Dmax=comb
        di,sp=nl2.diam2()
        #per plottare confronto con vecchi diametri
        p.write("plot '"+ndf+"' u 0:1 w lp")
        fcnf1="mlff"+str(comb[1]).split(".")[1].zfill(3)+str(comb[0]).split(".")[0]+"mmax"+str(comb[2])+".txt"
        fcnf2="m2ff"+str(comb[1]).split(".")[1].zfill(3)+str(comb[0]).split(".")[0]+"mmax"+str(comb[2])+".txt"
        if os.path.exists(os.path.join("af_files",fcnf1)):
           p.write(",\\\n'"+os.path.join("af_files",fcnf1)+"' u 0:1 w lp")
        if os.path.exists(os.path.join("af_files",fcnf2)):
           p.write(",\\\n'"+os.path.join("af_files",fcnf2)+"' u 0:1 w lp")
        '''if os.path.exists(os.path.join("frlike",ndf)):
           p.write(",\\\n'"+os.path.join("frlike",ndf)+"' u 0:1 w lp")
        fcnf3="ml_ff"+str(comb[1]).split(".")[1].zfill(3)+"_"+str(comb[0]).split(".")[0]+"m_max"+str(comb[2])+".txt"
        if os.path.exists(os.path.join("extractedDiam",fcnf3)):
           p.write(",\\\n'"+os.path.join("extractedDiam",fcnf3)+"' u 0:1 w lp")
        fcnf4="ir_ff"+str(comb[1]).split(".")[1].zfill(3)+"_"+str(comb[0]).split(".")[0]+"m_max"+str(comb[2])+".txt"
        if os.path.exists(os.path.join("extractedDiam",fcnf4)):
           p.write(",\\\n'"+os.path.join("extractedDiam",fcnf4)+"' u 0:1 w lp")'''
        p.write("\n")
        p.write("pause -1\n\n")

        p.write("set term png\n")
        p.write("set out '"+ndf+".png'\n")
        p.write("replot\n")
        p.write("set term win\n set out\n\n")
    p.close()                  

if __name__=="x__main__":
    nl2=CSimpleCreaDiam("imp_OffAxis.txt")

    focal=[i*2.5 for i in range(8,12)]
    ff=[0.00,0.07,0.15]
    D=[60,70]
    #comblist=[[f1,f2,f3] for f1 in focal for f2 in ff for f3 in D]

    comblist=[[20.0,0.0,60]]
    for comb in comblist:
        ndf="F"+str(comb[0]).split(".")[0]+"_ff"+str(comb[1]).split(".")[1].zfill(3)+"_D"+str(comb[2])+".txt"
        nl2.FL,nl2.ff,nl2.Dmax=comb
        nl2.calcGeo()
        di,sp=nl2.d,nl2.spes
        f= open(ndf,"w")
        for d in di[-1:0:-1]:
            f.write(str(d)+"\n")
        print ndf,nl2.FL,nl2.SH,nl2.Dmax,nl2.ff,nl2.nshell,nl2.wd,(nl2.thickVal[0],nl2.thickVal[1])
        print "\n"
        f.close()
    
    #plottaConfronto(nl2,comblist)

if __name__=="x__main__":

    fpFile="psf_data_02.txt"
    from userKov.pyGeneralRoutines.generalRoutines import loadCol

    angRT=loadCol(fpFile,9)
    angRT=[cc for cc in angRT if cc !=0]
    
    tel=CSimpleCreaDiam("imp_OffAxis.txt")
    print tel.d
    ad=tel.a1Distr(0,0.2)
    x=arange(100)*0.01*2*math.pi
    #plot(x,ad(x))
    #show()
    from userKov.pyOffAxis.focalPlane import hystogram, plotDist
    x2,y2=plotDist(fpFile)
    xx2,yy2=hystogram([ad(xx) for xx in x])
    plot(x2,y2,xx2,yy2)
    show()
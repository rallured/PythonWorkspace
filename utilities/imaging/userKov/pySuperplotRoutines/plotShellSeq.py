def plotShellSeq(folderlist,nshell,monolayerFlag,outfile='plotComp.plt'):
    '''
    given a list of folders containing the results from superplot,
    create a gnuplot script that
    plots the comparison betwee the effective area of each shell.
    '''
    l=['set datafile fortran']
    coloffset=[4 if (monoflag==1) else 2 for monoflag in monolayerFlag]
    for shn in range(1,nshell+1):
        l.append("plot\\")
        for fold,col in zip(folderlist[0:-1],coloffset[0:-1]):
            l.append("'"+fold+"\\%03i\\"%(shn)+"aeff%04i"%(shn)+".dat' u 1:%i w l,\\"%col)
        l.append("'"+folderlist[-1]+"\\%03i\\"%(shn)+"aeff%04i"%(shn)+".dat' u 1:%i w l"%coloffset[-1])
        l.append("pause -1")
        l.append("\n")
    #print "\n".join(l)

    fplt=open(outfile,'w')
    fplt.write("\n".join(l))
    fplt.close()

if __name__=="__main__":
##    plotShellSeq([r"D:\workOpt\aG2.1\hexitSatOpt1\run1_3plt\AeffxE2040",
##                                       r"D:\workOpt\aG2.1\hexitSatOpt1\run1_3plt\AeffxE4060",
##                                       r"D:\workOpt\aG2.1\hexitSatOpt1\run1_3plt\AeffxE6080",
##                                       r"D:\workOpt\aG2.1\hexitSatOpt1\run1_3plt\AeffxE1040"],3)
    plotShellSeq([r"E:\work\usb_F\Nuova\src_plt\NHXM_descoping\modulated_ptC50_1plt\upperlimit",
                  r"E:\work\usb_F\Nuova\src_plt\NHXM_descoping\modulated_mlphB_1plt\upperlimit",
                  r"E:\work\usb_F\Nuova\src_plt\NHXM_descoping\modulated_mlphB_1plt\upperlimit"],70,
                 [0,0,1])

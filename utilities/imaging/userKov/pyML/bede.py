#!/usr/bin/env python

from numpy import *
import os
import sys
import getopt

ver="### bede.py v1.1 5/5/07 Vincenzo Cotroneo OAB-INAF ###"
'''v1.1 aggiunto opzione cut'''

class bedeMeasure(object):
    '''misura al riflettometro, derivante da un file di misure'''
    def normalize(self):
        '''rinormalizza il fascio, usando il valore di self.beam che deve
        essere gia' stato impostato (vedi normalization)'''
        self.reflex=self.counts/self.beam
    def __init__(self,fileName,beam = -1):
        data=open(fileName,'r').read().split("\n") #legge tutte le righe
        self.info=data[0:31]    #le prime 30 sono commenti
        data=[map(float,d.split()) for d in data[31:] if d]  #le altre sono su due colonne
        data=array(data)
        self.angles= - data[:,0]
        self.counts= data[:,1]
        self.file=fileName
        if beam==-1:
            self.beam=2*self.counts[0]
            self.normalization="normalized to 2*counts(0)=%s"%(2*self.counts[0])            
        else:
            self.beam=float(beam)
            self.normalization="normalized to full beam=%s"%(self.beam)
        self.normalize()
        object.__init__(self)
    def write(self,outFileName=None,cutOff=None):
        if (not outFileName): outFileName=self.file+".dat"
        f=open(outFileName,"w")
        for a,r in zip(self.angles,self.reflex):
            if (not cutOff) or a > cutOff:            
                f.write("%s\t%s\n"%(a,r))
        f.close()
    def __repr__(self):
        l=[]
        sep="-------------------------"
        l.append("\n"+sep)
        l.append(object.__repr__(self))        
        l.append("%s values from file:\n'%s'"%(len(self.reflex),self.file))
        l.append(self.normalization)
        l.append(sep)
        l.append("%s\t%s\n%s\t%s\n...\n%s\t%s\n%s\t%s"%(self.angles[0],self.reflex[0],
                                             self.angles[1],self.reflex[1],
                                             self.angles[-2],self.reflex[-2],
                                             self.angles[-1],self.reflex[-1]))
        l.append(sep)
        return "\n".join(l)
        
            
if __name__=="__main__":
    '''Normalize measure data from BEDE X-ray reflectometer.
    \nusage:
    bede.py (-l listFile | dataFile, [options], [counts]).\n
    options:
    -l listFile:
        Use a file containing the list of datafile to process and (optionally)
        a 2nd and 3rd columns with full beam counts and name for output file.
        If fullBeamCounts is missing, prompt the user. If it is -1 then autocalculate.
    -a  Autocalculate the full beam using double of the first count value (same
        as giving a negative number of counts).
    -o outputFile:
        Name for the outputFile (default: add .dat to the original full name).
    -c angleMin:
        Cut the result file, giving only angles greater than angleMin.
    -t Test with a standard file (ignore all other options). Don't generate output.
        Show some results.
    '''

    def usage():
        s = ver+'''\nNormalize measure data from BEDE X-ray reflectometer.
    \nusage:
    bede.py (-l listFile | dataFile, [options], [counts]).\n
    options:
    -l listFile:
        Use a file containing the list of datafile to process and (optionally)
        a 2nd and 3rd columns with full beam counts and name for output file.
        If fullBeamCounts is missing, prompt the user. If it is -1 then autocalculate.
    -a  Autocalculate the full beam using double of the first count value (same
        as giving a negative number of counts).
    -o outputFile:
        Name for the outputFile (default: add .dat to the original full name).
    -c angleMin:
        Cut the result file, giving only angles greater than angleMin.
    -t Test with a standard file (ignore all other options). Don't generate output.
        Show some results.
    '''
        return s
        
    def norm(data,beam=None,outFile=None,cutOff=None):
        '''write a normalized file using the bedeMeasure class. beam is the
        number of counts of the full beam. if missing, ask to the user, if <0 autocalculate.
        the default output file has the same name of the datafile with extension .dat'''
        
        if not beam:
            print "full beam counts missing for datafile"
            print "'%s',"%data
            print "insert (or <=0 for autocalculation from 2*firstValue)"
            e=1
            while e:
                i=raw_input()
                try:
                    beam=float(i)
                    e=0
                except:
                    e=1
        elif beam <=0: beam=-1
        m=bedeMeasure(data,beam)
        m.write(outFile,cutOff)
        
    #----------- MAIN ---------------
    try:                                
        opts, args = getopt.gnu_getopt(sys.argv[1:], "o:l:c:ta", ["output=", "listFile=","test","cut=",
                                                                "autonorm"])
    except getopt.GetoptError:
        print usage()
        print
        print "given flags/values: ",opts
        print "given arguments: ", args
        sys.exit(-1)

    if not (opts or args):
        print usage()
        sys.exit()
    fileList=""
    auto=0
    beam=0
    cutOff=None
    outFile=None
    for o,a in opts:
        if o in ("-c","--cut"):
            cutOff=float(a)
        elif o in ("-l","--filelist"):
            if args:
                print "file name and filelist can't be set both."
                print "filename: %s\nfilelist: %s "%(args[0],a)
                sys.exit("filenameConflict")
            else:
                fileList=o
        elif o in ("-o","--output"):
            outFile=a
        elif o in ("-t","--test"):
            print "TEST VERSION!!"
            testFile=os.path.join(os.path.split(sys.argv[0])[0],"Test.x66")
            m=bedeMeasure(testFile)
            print m
            sys.exit()
        elif o in ("-a","--autonorm"):
            auto=1

    if fileList:
        lista=open(sys.argv[2],"r").read().split("\n")
        for s in lista:
            try: s= s.split()
            except:pass
            if auto: beam=-1
            elif len(s)>1: beam=float(s[1])
            if len(s)>2: outFile=s[2]
            norm(s[0],beam=beam,outFile=outFile,cutOff=cutOff)
    else:
        if auto:
            if len(args)>=2:
                print "beam counts can't be set at the same time to auto (-a flag) and to a value (%s)"%args[1]
                sys.exit("beamCountConflict")
            else: beam=-1
        else:
            if len(args)>=2 : beam=args[1]
        norm(args[0],beam=beam,outFile=outFile,cutOff=cutOff)

   
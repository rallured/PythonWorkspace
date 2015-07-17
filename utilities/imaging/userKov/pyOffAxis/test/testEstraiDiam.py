from pyOffAxis import estraiDiam 


def testParseTransRayTrace():
    print "\n------------------------------------"
    print "testo la routine parseTransRayTrace:"
    print estraiDiam.parseTransRayTrace.__doc__
    v,l=estraiDiam.parseTransRayTrace(os.path.curdir,"testTRT.txt")
    print "lunghezza lista labels: %s" %len(l)
    print "lunghezza lista valori: %s" %len(v)
    print "labels %s" %l
    
def testEstraiDiamDic():
    print "\n------------------------------------"
    print "testo la routine estraiDiamDic:"
    print estraiDiam.estraiDiamDic.__doc__
    listaTestDir=["testdir1","testdir2"]
    print estraiDiam.estraiDiamDic (listaTestDir)

def testEstraiDiam():
    print "\n------------------------------------"
    print "testo la routine estraiDiam:"
    print estraiDiam.estraiDiam.__doc__
    print "\n"
    print estraiDiam.estraiDiam ("testDirList.txt")    
    

    
if __name__=="__main__":
    testParseTransRayTrace()
    testEstraiDiamDic()
    testEstraiDiam()



import os

if  __name__=="__main__":
    bat=open("ps2bmp.bat","w")
    lista=[os.path.splitext(f)[0] for f in os.listdir("integrali_1") if os.path.splitext(f)[1]==".ps"]
    for f in lista:
        bat.write("gs -sDEVICE=bmp256 -sOutputFile=%s.bmp %s.ps\n"%(f,f))
    bat.close()
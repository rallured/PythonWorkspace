def isNum(a):
    res=True
    for i in a:                        
        try: float(i)
        except: res=False
    return res

print isNum([1,2,3])
print isNum([1,2,"dd"])
print isNum([2,"2.2","3d5"])  
print isNum([])  
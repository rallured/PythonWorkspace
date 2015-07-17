a=[["la Giorgia","Rodolfo"],[" arriva"," caga"],[" sul cesso"," in ufficio"],[" e scorreggia","."]]
b=(0,1)
c=(1)

print "itero su b"
noncisono=[i for i in range(4) if i not in b ]
stato=[1,2,3,4]
for i in enumerate(b):
    stato[i]=a[i]
    
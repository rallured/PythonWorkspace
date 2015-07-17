pro testIDL



ndim=5
file='risultati00001.dat'

;RESTORE,"risultati.tpl"
dati=READ_ASCII(FILEPATH(file,$
    SUBDIR =['..\..\..\work\pythonAnalyzer']))
;PRINT, dati
;PLOT, dati.FIELD1,dati.FIELD6,PSYM=1



selData=dati.field01[0:ndim+1,ndim+1:*:2*(ndim+1)]

;plot, selData[0,*],seldata[1,*],PSYM=1

weights = CLUST_WTS(seldata[*,*],n_iterations=100,/double)
Result = CLUSTER( seldata[*,*], Weights,/double)


;result[1:10]=0
;result[11:25]=1
;result [26:30]=2
;result [31:*]=3

print, result[*]

plot, selData[0,*],seldata[5,*],/nodata
plots, selData[0,*],seldata[5,*],PSYM=1,color=(result+1)*5000

;plot, selData[0,*],seldata[5,*],/nodata
;for i = 0,3 do begin
;	subarr = where(result eq i)
;	oplot, selData[0,subarr],seldata[5,subarr],PSYM=1,color=i*3000
;endfor

END
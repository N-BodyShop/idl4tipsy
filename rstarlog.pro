function rstarlog,file,VERBOSE = verbose
;;; RSTARLOG:  Tipsy starlog reader for IDL
;;; Author:  GS modified from original rtipsy by James Wadsley
;;; 
if (N_PARAMS() eq 0) then begin
  print, "rstarlog.pro  Reads tipsy starlog files: "
  print
  print, "Usage: "
  print, "        sl = rstarlog(filename [,/BIG] [,/VERBOSE]"
  print
  print, "Input parameters: "
  print, "  filename  filename string"
  print, "  /VERBOSE  print messages (optional)"
  print, "Return values:"
  print, "  sl    tipsy starlog struct array"
  print, "Please read rstarlog.pro for the structure definitions"
  print
  print, "Example: "
  print, "  sl = rstarlog('mysimulation.starlog')"
  print, "  plot, sl.rform[0], sl.rform[1], psym=3"
  return,-1
endif

if (file_test(file) eq 0) then file = file+".starlog"

; Open XDR file for reading
openr,lun,file,/xdr,/get_lun

; starlog record structure
;record = {iOrderStar:0L, iOrderGas:0L, timeform:0.d, rform:dblarr(3), vform:dblarr(3), massForm:0.0d, rhoform:0.0d, Tempform:0.0d}
record = {iOrderStar:0L, iOrderGas:0L, timeform:0.d, x:0.d,y:0.d,z:0.d,vx:0.d, $
	vy:0.d, vz:0.d, massForm:0.0d, rhoform:0.0d, Tempform:0.0d}

iSize = 0L
readu,lun,iSize

; Error check the integer we just read from the starlog file to see if it 
; matches our record.
if (n_bytes(record) ne iSize) then begin
  print, "Record size in IDL ("+string(n_bytes(record))+") does not " $
	  + "match size read from file ("+string(iSize)+")"
  return,-1
endif

; Find out how many records there are
fs = fstat(lun)
n=(fs.size-4)/iSize

if (keyword_set(verbose)) then print,"Reading in ",n," stars."

cats = replicate(record,n)
readu,lun,cats
  
close,lun
return,cats[uniq(cats.iorderstar,sort(cats.iorderstar))]

end

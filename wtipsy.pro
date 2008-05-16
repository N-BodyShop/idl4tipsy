PRO wtipsy,outfile,header,catg,catd,cats,NATIVE=native
;This program takes structures like those read in via rtipsy, and
;writes them to native binary output or standard binary output if
;/STANDARD is set.  This assumes you're working on a linux box

if (N_PARAMS() eq 0) then begin
  print, "wtipsy.pro -- Writes tipsy files with structures as read in"
  print, "using rtipsy.pro."
  print, ""
  print, "Usage: "
  print, "        wtipsy, outfilename ,header ,g ,d ,s, [/STANDARD]"
  print, ""
  print, "Input parameters: "
  print, "  outfilename  string containing name of output file"
  print, "  g,d,s     gas, dark and star structures"  
  print, "            if no gas or stars just include a dummy variable"    
  print, "Please read rtipsy.pro for the structure definitions"
  print, "  /NATIVE will write the output in native binary format"
  return
endif

; Standard output needs an 8 byte pad (dummy) in the header
;header = { time:double(0.0), n:ngas+ndark, ndim:3L, ngas:ngas, $
;	           ndark:ndark, nstar:nstar}
; For reference
;catg = replicate({mass: 1.,x: 1.,y : 1., z:1.,vx:1.,vy:1.,vz:1.,dens:1.,tempg:1.,h : 1. , zmetal : 1., phi : 1.},header.ngas)
;catd = replicate({mass: 1.,x: 1.,y : 1., z:1.,vx:1.,vy:1.,vz:1.,eps: 1.,phi: 1.},header.ndark)
;cats = replicate({mass: 1.,x: 1.,y : 1., z:1.,vx:1.,vy:1.,vz:1.,metals:1.,tform:1.,eps: 1.,phi: 1.},header.nstar)

IF (keyword_set(native)) THEN BEGIN
;NATIVE
    OPENW,lun,outfile,/get_lun
    WRITEU,lun,header
ENDIF ELSE BEGIN
;STANDARD
    OPENW,lun,outfile,/xdr,/get_lun
    dummy=1L
    WRITEU,lun,header,dummy
ENDELSE

if(keyword_set(catg)) then writeu,lun,catg
if(keyword_set(catd)) then writeu,lun,catd
if(keyword_set(cats)) then writeu,lun,cats

CLOSE,1
FREE_LUN, 1

END

PRO wtipsy,outfile,header,catg,catd,cats,STANDARD=standard
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
  print, "  /STANDARD will write the output in standard format if you're using a i386 box"
  return
endif

; For reference
;catg = replicate({mass: 1.,x: 1.,y : 1., z:1.,vx:1.,vy:1.,vz:1.,dens:1.,tempg:1.,h : 1. , zmetal : 1., phi : 1.},header.ngas)
;catd = replicate({mass: 1.,x: 1.,y : 1., z:1.,vx:1.,vy:1.,vz:1.,eps: 1.,phi: 1.},header.ndark)
;cats = replicate({mass: 1.,x: 1.,y : 1., z:1.,vx:1.,vy:1.,vz:1.,metals:1.,tform:1.,eps: 1.,phi: 1.},header.nstar)

IF (keyword_set(standard) EQ 0) THEN BEGIN
;NATIVE
    OPENW,1,outfile
    WRITEU,1,header,catg,catd,cats
ENDIF ELSE BEGIN
;STANDARD
    OPENW,1,outfile,/xdr
    dummy=1L
    WRITEU,1,header,dummy,catg,catd,cats
ENDELSE

CLOSE,1
FREE_LUN, 1

END

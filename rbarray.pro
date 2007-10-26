function rbarray,file,LONG=long,DOUBLE=double, VERBOSE = verbose,
		 SWAP=swap,MAXSIZE=maxsize
;;; RTIPSY:  Tipsy reader for IDL
;;; Author:  James Wadsley
;;; 
if (N_PARAMS() eq 0) then begin
  print, "rbarray.pro  Reads single element binary array  files detecting the format: "
  print, "big endian, little endian, padded (standard) or non-padded header "
  print
  print, "Usage: "
  print, "        rbarray( filename [,/LONG][,/DOUBLE] [,/VERBOSE]"
  print, "                 [,/SWAP][,MAXSIZE=maxsize])"
  print
  print, "Input parameters: "
  print, "  filename  filename string"
  print, "  /VERBOSE  print messages (optional)"
  print, "  /LONG     sets format to integer"
  print, "  /DOUBLE   sets format to double precision"
  print, "  /SWAP     forces endian swapping"
  print, "  maxsize   sets test condition for endian swapping"
  print, "Return values:"
  print, "  The array"
  print, "Please read rbarray.pro for the structure definitions"
  print
  print, "Example: "
  print, "  array = rbarray( 'mysimulation.iord',/LONG)"
  print, "  print, n_elements(array)"
  return,-1
endif

if ( keyword_set(maxsize) eq 0 ) then maxsize = 100000000L
if ( keyword_set(long) ) then type='long' $
else if (keyword_set(double)) then type='double' $
else type='float'
;;; Note: IDL structures are never paddded 
arrsize= 0
dummy = arrsize

close,1
openr,1,file

;readu,1,arrsize
; Check how long long integers are (32 or 64 bit)
if ( arrsize eq 0 ) then begin
  point_lun,1,0
  ;arrsize=LONG64(0)
  arrsize=0L
  readu,1,arrsize
  verylong=1
endif else readu,1,dummy

; swap_endian if /swap is set.  Kind of clever trick to recognize when you
; want to do this 
endianswap = 0
if ( keyword_set(swap) OR ( arrsize lt 0 OR arrsize gt MAXSIZE ) ) then begin
  endianswap = 1
  arrsize=swap_endian(arrsize)
  if (keyword_set(verbose)) then print,"SWAP_ENDIAN"
endif

; create dummy array that will hold the contents of the file
CASE type OF
'float': blah=fltarr(arrsize)
'double': blah=dblarr(arrsize)
ELSE: blah=LONARR(arrsize)
ENDCASE

; read the file
readu,1,blah

if (endianswap eq 1) then blah=swap_endian(blah)

close,1
return,blah

end

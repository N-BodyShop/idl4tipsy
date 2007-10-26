function rbvector,file,TYPE=type, TIME = time,VERBOSE = verbose,SWAP=swap,MAXSIZE=maxsize
;;; RTIPSY:  Tipsy reader for IDL
;;; Author:  James Wadsley
;;; 
if (N_PARAMS() eq 0) then begin
  print, "rbvector.pro  Reads 3 element vector binary array tipsy files detecting the format: "
  print, "big endian, little endian, padded (standard) or non-padded header "
  print
  print, "Usage: "
  print, "        rbvector( filename [,TYPE=type] [,TIME=time] [,MAXSIZE=maxsize] [,/VERBOSE] [,/SWAP])"
  print
  print, "Input parameters: "
  print, "  filename  filename string"
  print, "  time      desired output time (optional)"
  print, "  type      file format: double, float[default], long (optional)"
  print, "  maxsize   number used in checking endianness, def=1e7 (optional)"
  print, "  /VERBOSE  print messages (optional)"
  print, "  /SWAP     swap endianness (optional)"
  print, "Return values:"
  print, "  The array"
  print, "Please read rbvector.pro for the structure definitions"
  print
  print, "Example: "
  print, "  array = rbvector( '/net/mega-2/stinson/volumes/22/22.00256.rform')"
  print, "  print, n_elements(array)"
  return,-1
endif

if ( keyword_set(maxsize) eq 0 ) then maxsize = 10000000L
if ( keyword_set(type) eq 0 ) then type='float'
;;; Note: IDL structures are never paddded 
arrsize= 0L

close,1
openr,1,file

Loop:  

readu,1,arrsize

; Check how long long integers are (32 or 64 bit)
if ( arrsize eq 0 ) then begin
  point_lun,1,0
  arrsize=LONG64(0)
  readu,1,arrsize
  verylong=1
endif

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
'long': BEGIN
x=LONARR(arrsize)
y=LONARR(arrsize)
z=LONARR(arrsize)
  END
'double': BEGIN
x=DBLARR(arrsize)
y=DBLARR(arrsize)
z=DBLARR(arrsize)
  END
ELSE: BEGIN
x=FLTARR(arrsize)
y=FLTARR(arrsize)
z=FLTARR(arrsize)
  END
ENDCASE

; read the file
readu,1,x
readu,1,y
readu,1,z
blah={x:x,y:y,z:z}

if (endianswap eq 1) then blah=swap_endian(blah)

;;; Loop over output times if requested
if (keyword_set(time)) then begin
  if (abs(time-header.time) gt 1e-3) then begin
    on_ioerror, ReadError
    goto, Loop
  endif
endif

close,1
return,blah

ReadError:
print,"RBVECTOR ERROR: Output time not found ",time
on_ioerror,NULL

close,1
return, -1

end

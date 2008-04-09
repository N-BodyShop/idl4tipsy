pro galactics2tipsy,output=output,softening=softening,BINARY=binary, $
	NODARK=nodark,VERBOSE = verbose
;;; Galactics2tipsy:  convert Dubinski ICs to tipsy
;;; Author:  GS modified from original rtipsy by James Wadsley
;;; 
if (N_PARAMS() LT 0) then begin
  print, "galactics2tipsy.pro  Converts Dubinski ICs to tipsy: "
  print
  print, "Usage: "
  print, "        galactics2tipsy,filename[,tipsyfile] [,/VERBOSE]"
  print
  print, "Input parameters: "
  print, "  output     output tipsy filename string"
  print, "  softening  what to set the gravitational softening"
  print, "  /BINARY    read in binary formatted files"
  print, "  /NODARK    don't put dark particles into tipsy file"
  print, "  /VERBOSE   print messages (optional)"
  print
  print, "Example: "
  print, "  galactics2tipsy,'galaxy','galaxy.std'"
  return
endif

tempdir = '/home/stinson/idl'
if(keyword_set(binary) EQ 0) then restore,tempdir+'/galacticstemp.sav'
if(keyword_set(nodark) EQ 0) then begin
print,'reading halo'
if(keyword_set(binary) EQ 0) then halo = read_ascii('halo',template=temp) $
else begin
  ; Open galactics file for reading
  openr,lun,'halo',/get_lun
  fs = fstat(lun)
  ; Each record consists of 7 4-byte items (mass + postion(3) + velocity(3))
  nDark=(fs.size)/28
  halo = replicate({mass: 1.,x: 1.,y : 1., z:1.,vx:1.,vy:1.,vz:1.}, nDark)
  readu,lun,halo
  close,lun
endelse
endif

print,'reading disk'
if(keyword_set(binary) EQ 0) then disk = read_ascii('disk',template=temp) $
else begin
  openr,lun,'disk',/get_lun
  fs = fstat(lun)
  ; There are 2 header numbers, nDisk and time (see above)
  ; Then each record consists of 7 4-byte items (mass + postion + velocity)
  nDisk=(fs.size-8)/28
  nd = 0L
  time = 0.0
  readu,lun,nd,time
  if (nd ne nDisk) then printf,"Galactics header doesn't match file size."
  disk = replicate({mass: 1.,x: 1.,y : 1., z:1.,vx:1.,vy:1.,vz:1.}, nDisk)
  readu,lun,disk
  close,lun
endelse

print,'reading bulge'
if(keyword_set(binary) EQ 0) then bulge = read_ascii('bulge',template=temp) $
else begin
  openr,lun,'bulge',/get_lun
  fs = fstat(lun)
  ; Each record consists of 7 4-byte items (mass + postion(3) + velocity(3))
  nBulge=(fs.size)/28
  bulge = replicate({mass: 1.,x: 1.,y : 1., z:1.,vx:1.,vy:1.,vz:1.}, nBulge)
  readu,lun,bulge
  close,lun
endelse

; Find out how many records there are
nStar = n_elements(disk.mass) + n_elements(bulge.mass)
if (keyword_set(nodark) EQ 0) then begin
  nDark = n_elements(halo.mass)
  n=nDark + nStar
endif else begin
  n = nStar
  nDark = 0L
endelse

if(keyword_set(time) EQ 0) then time = 0.0

header = { time:double(time), n:n, ndim:3L, ngas:0L, ndark:nDark, nstar:nStar, dummy:1L }
if(keyword_set(nodark) EQ 0) then begin
  catd = replicate({mass: 0.,x: 0.,y : 0., z:0.,vx:0.,vy:0.,vz:0., $
                  eps: softening,phi: 0.},header.ndark)
endif
cats = replicate({mass: 0.,x: 0.,y : 0., z:0.,vx:0.,vy:0.,vz:0., $
                  metals:0.02,tform:0.,eps: softening,phi: 1.},header.nstar)
;catg = replicate({mass: 1.,x: 1.,y : 1., z:1.,vx:1.,vy:1.,vz:1., $
;                  dens:1.,tempg:1.,h : 1. , zmetal : 1., phi : 1.}, $
;	          header.ngas)

if(keyword_set(nodark) EQ 0) then begin
catd.mass = halo.mass;*1e4
catd.x = halo.x
catd.y = halo.y
catd.z = halo.z
catd.vx = halo.vx;*1e2
catd.vy = halo.vy;*1e2
catd.vz = halo.vz;*1e2
endif

cats.mass = [disk.mass,bulge.mass];*1e4
cats.x = [disk.x,bulge.x]
cats.y = [disk.y,bulge.y]
cats.z = [disk.z,bulge.z]
cats.vx = [disk.vx,bulge.vx];*1e2
cats.vy = [disk.vy,bulge.vy];*1e2
cats.vz = [disk.vz,bulge.vz];*1e2

if (keyword_set(output) eq 0) then output = 'galaxy.std'
print,'writing '+output
openw,lun,output,/get_lun,/xdr
if(keyword_set(nodark)) then writeu,lun,header,cats $
else writeu,lun,header,catd,cats
close,lun

end

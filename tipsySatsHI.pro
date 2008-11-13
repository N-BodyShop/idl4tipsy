pro tipsysatshi, infile, groups, outfile=outfile, h=h,g=g,d=d,s=s, gethi=gethi, $
               ascii=ascii

if not keyword_set(outfile) then outfile = infile

grpfile = infile+'.amiga.grp'
gtpfile = infile+'.amiga.gtp'
grp = read_ascii_array(grpfile)

;cmp_file = filename+'.cmp'

rtipsy, infile, h,g,d,s
if keyword_set(gethi) then hi = (read_ascii_array(infile+'.HI'))

;stop
;print,"Read tipsy file."


;*******************************************
;INDEX GAS,DM & STARS FROM PARTICULAR HALO
;*******************************************
;for i = 1,ngroups do begin

for m = 0, n_elements(groups)-1 do begin

i = groups[m]

ind = where(grp eq i,comp=indcomp)
inds = ind(where(ind ge h.ngas+h.ndark))
indg = ind(where(ind lt h.ngas))
indd = ind(where(ind ge h.ngas and ind lt h.ngas+h.ndark))
;tempdata = read_ascii_array(cmp_file)
;components = tempdata[h.ngas+h.ndark:h.ndark+h.nstar+h.ngas-1]
;stop

stars = s[inds-h.ngas-h.ndark]
gas = g[indg]
dark = d[indd-h.ngas] 
;***************************************************************************
; get CENTRE OF MASS from amiga (first entry in gtp file): REPOSITION
;**************************************************************************
;rtipsy, gtpfile, h1,g1,d1,s1

;cx= s1[i-1].x
;cy= s1[i-1].y
;cz= s1[i-1].z

foo = min(s.phi, cmindex)
cx = s[cmindex].x
cy = s[cmindex].y
cz = s[cmindex].z

gas.x = gas.x - cx
gas.y = gas.y - cy
gas.z = gas.z - cz
dark.x = dark.x - cx
dark.y = dark.y - cy
dark.z = dark.z - cz

stars.x =  stars.x - cx
stars.y =  stars.y - cy
stars.z =  stars.z - cz


;***********************************************
;eliminate proper motion of the halo/galaxy
;***********************************************
propx=mean(stars.vx)
propy=mean(stars.vy)
propz=mean(stars.vz)

stars.vx=stars.vx-propx
stars.vy=stars.vy-propy
stars.vz=stars.vz-propz
gas.vx=gas.vx-propx
gas.vy=gas.vy-propy
gas.vz=gas.vz-propz
dark.vx=dark.vx-propx
dark.vy=dark.vy-propy
dark.vz=dark.vz-propz

;**************************************************************
; find rotation curve and transtarslate in XY plane using align.pro
;*************************************************************
dist_units=1e5*h.time
dist_units=5.0e5*h.time
limit=5./dist_units            ;use stars within this limit (kpc/dist_units)
align,stars,dark,gas,limit


;sout = where(components ge 1 and components le 4)

;radius = sqrt(stars.x^2+stars.y^2+stars.z^2)
;limit = max(radius)
;gout = where (sqrt(gas.x*gas.x + gas.y*gas.y +gas.z*gas.z) lt limit)
;dout = where (sqrt(dark.x*dark.x + dark.y*dark.y +dark.z*dark.z) lt limit)


num = h 
num.nstar = n_elements(stars)
num.ngas = n_elements(gas)
num.ndark = n_elements(dark)
num.n = num.nstar + num.ngas + num.ndark
;stop
fileout = outfile + '.halo.' + strtrim(string(i),2) + '.std'
;fileout='MW'+strmid(filename,5,6,/reverse_offset)+'_'+strtrim(string(i),2)+'.asc'
if keyword_set(ascii) then $
  wtipsy, fileout, num, gas, dark, stars $
  else wtipsy,fileout,num,gas,dark,stars,/standard

if keyword_set(gethi) then begin
    fileout = outfile + '.halo.' + strtrim(string(i),2) + '.HI'
    openw,lunhi,fileout,/get_lun
    printf,lunhi,num.n
    for j=0L,n_elements(ind)-1 do begin
        printf,lunhi,hi[ind[j]]
    endfor
    close,lunhi
endif


endfor
;stop
end


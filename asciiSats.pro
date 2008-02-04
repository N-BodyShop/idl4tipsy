; cuts out AMIGA halos and writes them to tipsy ascii files 
; just input filename and the number of halo/galaxies you want to output... 
; galaxies are sequentially ordered from the most massive
; you will need to modify the name of the outfile you want to write to
pro asciisats, filename, ngroups, h=h,g=g,d=d,s=s

grpfile = filename+'.amiga.grp'
gtpfile = filename+'.amiga.gtp'
grp = read_lon_array(grpfile)
rtipsy, filename, h,g,d,s

;*******************************************
;INDEX GAS,DM & STARS FROM PARTICULAR HALO
;*******************************************
for i = 1,ngroups do begin

ind = where(grp eq i,comp=indcomp)
inds = ind[where(ind ge h.ngas+h.ndark)]-h.ngas-h.ndark
indg = ind[where(ind lt h.ngas)]
indd = ind[where(ind ge h.ngas and ind lt h.ngas+h.ndark)]-h.ngas

stars = s[inds]
gas = g[indg]
dark = d[indd] 
;***************************************************************************
; get CENTRE OF MASS from amiga (first entry in gtp file): REPOSITION
;**************************************************************************
rtipsy, gtpfile, h1,g1,d1,s1

cx= s1[i-1].x
cy= s1[i-1].y
cz= s1[i-1].z

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
; do this if you want to align your galaxies
;*************************************************************
;dist_units=2.85714e4*h.time ; this is simulation dependant!
;limit=5./dist_units            ;use stars within this limit (kpc/dist_units)
;align,stars,dark,gas,limit


num = h 
num.nstar = n_elements(stars)
num.ngas = n_elements(gas)
num.ndark = n_elements(dark)
num.n = num.nstar + num.ngas + num.ndark


fileout='036'+strmid(filename,5,6,/reverse_offset)+'_'+strtrim(string(i),2)+'.asc'  ; change the out file naming convention to something relevant
wtipsy_asc,fileout,num,gas,dark,stars
endfor
stop
end


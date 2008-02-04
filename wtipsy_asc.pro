;write a tipsy format ascii file 
pro wtipsy_asc,outfile,h,g,d,s

ndark=h.n-h.nstar-h.ngas
;***********************************************************************
;write the asci file in tipsy format that can be read by galaxy_scope
;***********************************************************************
;outfile='grp'+ntostr(group)+'.asc'
openw,1,outfile
;openw,2,markfile

printf,1,h.n,h.ngas,h.nstar
printf,1,h.ndim
printf,1,h.time
;*****mass
j=0
for i=0L, h.ngas-1 do begin
printf,1,g[i].mass ;[,format='(F10.3)'] or whatever the format is
j=j+1
;printf,2,j
endfor
for i=0L,ndark-1 do begin
printf,1,d[i].mass 
j=j+1
;printf,2,j
endfor
for i=0L,h.nstar-1 do begin
printf,1,s[i].mass 
j=j+1
;printf,2,j
endfor
;*****position
for i=0L,h.ngas-1 do begin
printf,1,g[i].x 
endfor
for i=0L,ndark-1 do begin
printf,1,d[i].x 
endfor
for i=0L, h.nstar-1 do begin
printf,1,s[i].x 
endfor
for i=0L, h.ngas-1 do begin
printf,1,g[i].y 
endfor
for i=0L, ndark-1 do begin
printf,1,d[i].y 
endfor
for i=0L, h.nstar-1 do begin
printf,1,s[i].y 
endfor
for i=0L, h.ngas-1 do begin
printf,1,g[i].z 
endfor
for i=0L, ndark-1 do begin
printf,1,d[i].z 
endfor
for i=0L, h.nstar-1 do begin
printf,1,s[i].z 
endfor
;****** velocity
for i=0L, h.ngas-1 do begin
printf,1,g[i].vx 
endfor
for i=0L, ndark-1 do begin
printf,1,d[i].vx 
endfor
for i=0L, h.nstar-1 do begin
printf,1,s[i].vx 
endfor
for i=0L, h.ngas-1 do begin
printf,1,g[i].vy 
endfor
for i=0L, ndark-1 do begin
printf,1,d[i].vy 
endfor
for i=0L, h.nstar-1 do begin
printf,1,s[i].vy 
endfor
for i=0L, h.ngas-1 do begin
printf,1,g[i].vz 
endfor
for i=0L, ndark-1 do begin
printf,1,d[i].vz 
endfor
for i=0L, h.nstar-1 do begin
printf,1,s[i].vz 
endfor
;;******** softening length
for i=0L, ndark-1 do begin
printf,1,d[i].eps
endfor
for i=0L, h.nstar-1 do begin
printf,1,s[i].eps 
endfor
;************density, temp, smoothing
for i=0L, h.ngas-1 do begin
printf,1,g[i].dens 
endfor
for i=0L, h.ngas-1 do begin
printf,1,g[i].tempg 
endfor
for i=0L, h.ngas-1 do begin
printf,1,g[i].h 
endfor
;**** metals
for i=0L, h.ngas-1 do begin
printf,1,g[i].zmetal 
endfor
for i=0L, h.nstar-1 do begin
printf,1,s[i].metals 
endfor
;**** formation time
for i=0L, h.nstar-1 do begin
printf,1,s[i].tform 
endfor
;******** potential
for i=0L, h.ngas-1 do begin
printf,1,g[i].phi
endfor
for i=0L, ndark-1 do begin
printf,1,d[i].phi 
endfor
for i=0L, h.nstar-1 do begin
printf,1,s[i].phi
endfor
close,1


return
end


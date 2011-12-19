; This program takes the output from AMIGA and creates .grp, .gtp and 
; .stat files similar to that produced by SKID. 
;
; Author: Alyson Brooks
; Last updated: July 26, 2009

pro v277_grp_stat, tipsyfile, boxsize=boxsize, munit = munit, vunit = vunit, h0 = h0, multiplefile = multiplefile

; INPUTS
;
; tipsyfile is the name of the tipsy file for the run.  If you are not running 
; 	this code in the simulation directory where tipsyfile is located, then 
; 	tipsyfile must include the full directory path
;
; IMPORTANT: this code assumes that the output file name from AMIGA is simply the 
;	tipsy file name + addenda from AMIGA.  If this is not the case, this code 
;	will not find your AMIGA results.  I searches for the .AHF_halos and 
;	.AHF_particles files. 
;
; boxsize: This should be the comoving length unit (boxsize) in kpc/h.   
;
; munit: the mass unit of the simulation.  If no input is given here, munit will 
; 	default to 3.171d15 
;
; vunit: the velocity unit of the simulation.  If no input is given here, munit will 
; 	default to 691.0 km/s
;
; h0: the value of h used in the simulation.  If no input is given, h0 will 
;	default to 0.73
;
; multiple: if you wish to run this code on more than one input, multiplefile is a 
;	file that contains a line for each input.  In this case, all inputs must 
;	be given (i.e., munit, vunit, and h0 are not optional inputs) for each 
;	input file

if (n_params() eq 0) then begin
  print, "mpi_grp_stat.pro -- Creates a .grp, .gtp, and .stat file from AMIGA outputs"
  print, " "
  print, "Usage: "
  print, "mpi_grp_stat, tipsyfile, boxsize, munit = munit, vunit = vunit, h0 = h0, multiplefile = multiplefile"
  print, "Boxsize in kpc/h"
  print, "Munit in M_sun, v_unit in km/s, h0 = H0/100.  "
  print, "Please read mpi_grp_stat.pro for more info on the inputs"
  print, " "
  print, " "
endif

IF keyword_set(multiplefile) THEN readcol, multiplefile, tipsyfile, format='a'
IF keyword_set(h0) eq 0 then h0 = 0.73 

; start a big loop over the number of input simulations
FOR j=0,n_elements(tipsyfile)-1 DO BEGIN

IF keyword_set(munit) eq 0 then begin
name = strmid(tipsyfile[j],5,7)
if name eq 'cosmo25' or name eq '.cosmo2' then begin
  boxsize = 18250.
  munit=2.310e15
  vunit = 630.28
endif
if name eq 'cosmo50' or name eq 'osmo50c' or name eq 'smo50cm' or name eq '50cmb.2' then begin
  boxsize = 36500. 
  munit=1.84793e16
  vunit = 1261.05
endif
;if name ne 'cosmo25' and name ne 'cosmo50' and name ne 'osmo50c' and name ne 'smo50cm' and name ne '.cosmo2' then begin
;  munit = 3.171d15
;  vunit = 691.0
;  h0 = 0.7
;  boxsize = 20000.
;endif  
ENDIF
lunit = boxsize/h0
print, 'Mass Unit = ', munit

rtipsy, tipsyfile[j], h,g,d,s
command = "ls "+tipsyfile+"*halos"
spawn,command,halosfile  ; return the files without the .amiga.grp
command = "ls "+tipsyfile+"*particles"
spawn,command,particlefile  ; return the files without the .amiga.grp


; The format of the particlefile is a long array.  The very first value is the 
; total number of halos in the file.  For each halo, the number 
; of particles in the halo is listed, and that many next lines contain the 
; index of the particles that belong to that halo.  When the first halo is 
; done, then the number of particles in the next halo is given and that many 
; next lines contain the indices to the particles in that halo, and so on.  So, 
; we must grab all the indices without treating npart of each halo as an index. 
readcol, halosfile, npart, nvpart, xc, yc, zc, vxc, vyc, vzc, mvir, rvir, vmax, rmax, sigv, spin, format='l,f,f,f,f,f,f,f,f,f,f,f,f,f'

; **************************************************************************
; First, create the .grp file
; **************************************************************************
; RDFLOAT is not working for the new particles file that has two columns (a 
; second column on particle type).  So, below I've resorted to using readcol
; because it actually works, but it is SLOooooooW.  There must be a better way!

;if h.n lt 4.2e6 then begin 
; rdfloat, particlefile, index, skipline=1
; index = long(index)
;endif else 
readcol, particlefile, index, format='l', skipline=1
init=0
startind = lonarr(n_elements(npart))
endind = lonarr(n_elements(npart))
for i=0L, n_elements(npart)-1 do begin
  addind = i+1 
  startind[i] = addind+init
  endind[i] = startind[i]+npart[i]-1 
  init = init+npart[i]
endfor

; Now we have the starting and ending point of each halo's particles in the 
; particlefile.  Let's put them into an array consecutively rather than 
; separated by npart.
; We want to remove halos that have been contaminated by low res dm particles
; here.  NOTE that there can be mild contamination in halos with a significant 
; number of baryonic particles (which have low particle masses).  We sometimes 
; want to know about these, as they are not exclusively in the res region.  
; Below we pick a certain fraction of the mass that is allowed to be in this 
; low res form (10%).  
; The .stat file contains a column to identify these contaminated halos.

particles = lonarr(n_elements(index))
grp = lonarr(n_elements(index))
mindm = min(d.mass)
for i=0L, n_elements(startind)-1 do begin
  halo = index[startind(i):endind(i)]
  dark_ind = where(halo LT h.ndark, ndark)
  if ndark NE 0 then begin
    dark_mass = total(d[halo(dark_ind)].mass)
    test2 = where(d[halo(dark_ind)].mass GT 1.5*mindm, ntest2)
    if ntest2 NE 0 then begin
      heavy_dark_mass = total(d[halo(dark_ind(test2))].mass)
      if (heavy_dark_mass/dark_mass GE 0.10) then continue else particles[startind(i):endind(i)] = index[startind(i):endind(i)]
    endif else particles[startind(i):endind(i)] = index[startind(i):endind(i)]
  endif else continue ;NOTE this means we DO NOT allow output of halos that contain only baryons
  grp[startind(i):endind(i)] = i+1
endfor
nonzero = where(particles ne 0)
particles = particles(nonzero)
grp = grp(nonzero) 
;AMIGA reorders the particles to dark, star, gas.  Put back in gas, dark, star order:
if h.ngas ne 0 then begin
gind = where(particles ge h.ndark+h.nstar+1, ng)
dind = where(particles lt h.ndark, nd)
sind = where(particles ge h.ndark and particles lt h.ndark+h.nstar, ns)
sparticles = lonarr(nd+ns+ng)
sgrp = lonarr(nd+ns+ng)
sparticles[0:ng-1] = particles[gind]-h.ndark-h.nstar+1
sgrp[0:ng-1] = grp[gind]
sparticles[ng:ng+nd-1] = particles[dind]+h.ngas
sgrp[ng:ng+nd-1] = grp[dind]
if h.nstar ne 0 then begin
 sparticles[nd+ng:nd+ng+ns-1] = particles[sind]+h.ngas 
 sgrp[nd+ng:ng+nd+ns-1] = grp[sind]
endif
endif else begin
sparticles = particles
sgrp = grp
endelse

; Now put the particle indices in ascending order, and get the grp number 
; that each particle belongs to.
ntot = long(h.n)
sortind = sort(sparticles)
sortedgrp = sgrp(sortind)
totarray = lonarr(ntot) 
ind = sparticles(sortind)
for i=0L,n_elements(sortind)-1 do totarray[ind(i)] = sortedgrp(i) 


; Totarray now contains the grp that each particle belongs to, and we 
; can write it to a .grp file.

; NOTE: some particles are assigned by AMIGA to more than one halo!  This
; means that some of the ind values that go into totarray get written over
; after an intial assignment of grp number.  This is ok if the halos are
; read in by order of mass, as that guarantees that the reassigned value
; is a substructure grp number.  However, AMIGA currently orders the halos
; by npart rather than my mass.  However, it is likely that the oddball
; halos that are out of mass order are lowres halos that accidentally made
; it into the catalog.  So, for now I will not rearrange the input by mass.
; This may become necessary at a later time, though.

; Write the .grp file
get_lun, lun
filename = strcompress(tipsyfile[j]+'.amiga.grp')
openw, lun, filename
printf, lun, ntot
for i=0L,n_elements(totarray)-1 do printf, lun, long(totarray[i])
close, lun
free_lun, lun

; **************************************************************************
; Now, create the .stat file
; **************************************************************************
sorted = sortedgrp(sort(sortedgrp))
unique = sorted(uniq(sorted))
ngrps = n_elements(unique)
amigaind = unique-1

stat = replicate({grp:0L, sat:strarr(1), contam:strarr(1), false:strarr(1), n_part:0L, n_gas:0L, n_star:0L, n_dark:0L, m_vir:dblarr(1), r_vir:dblarr(1), gas_mass:dblarr(1), star_mass:dblarr(1), dark_mass:dblarr(1), v_max:dblarr(1), r_max:dblarr(1), sig_v:dblarr(1), x_c:dblarr(1), y_c:dblarr(1), z_c:dblarr(1), vx_c:dblarr(1), vy_c:dblarr(1), vz_c:dblarr(1), amiga_orig_id:0L}, ngrps)

; Take care of the easy stuff first. 
; By definition, there should not be grp 0 
;for i=0L,ngrps-1 do stat[i].grp = i 
;for i=0L,ngrps-1 do stat[i].amiga_orig_id = amigaind[i]
for i=0L,ngrps-1 do stat[i].grp = unique[i] 
stat.n_part = npart(amigaind)
for i=0L,ngrps-1 do stat[i].m_vir = double(mvir(amigaind[i])/h0) 
for i=0L,ngrps-1 do stat[i].r_vir = double(rvir(amigaind[i])/h0)
for i=0L,ngrps-1 do stat[i].v_max = double(vmax(amigaind[i]))
for i=0L,ngrps-1 do stat[i].r_max = double(rmax(amigaind[i])/h0)
for i=0L,ngrps-1 do stat[i].sig_v = double(sigv(amigaind[i]))
for i=0L,ngrps-1 do stat[i].x_c = double(xc(amigaind[i])/h0)
for i=0L,ngrps-1 do stat[i].y_c = double(yc(amigaind[i])/h0)
for i=0L,ngrps-1 do stat[i].z_c = double(zc(amigaind[i])/h0)
for i=0L,ngrps-1 do stat[i].vx_c = double(vxc(amigaind[i]))
for i=0L,ngrps-1 do stat[i].vy_c = double(vyc(amigaind[i]))
for i=0L,ngrps-1 do stat[i].vz_c = double(vzc(amigaind[i]))
; Now figure out info on gas and stars, if they exist!
IF h.ngas+h.nstar GT 0 THEN BEGIN 
  gmass = double(g.mass[0])
  if h.nstar GT 0 then smass = double(s.mass[0])
  dmass = double(d.mass[0])
  FOR i=0L, ngrps-1 do begin
    particles_considered = index[startind(amigaind[i]):endind(amigaind[i])]
    gas_ind = where(particles_considered GT h.ndark+h.nstar+1)
    if gas_ind[0] NE -1 then stat[i].n_gas = n_elements(gas_ind) else stat[i].n_gas = 0
    if gas_ind[0] NE -1 then stat[i].gas_mass = total(gmass[particles_considered(gas_ind)-h.ndark-h.nstar+1]) else stat[i].gas_mass = 0
  IF h.nstar GT 0 then begin
    star_ind = where(particles_considered GE h.ndark and particles_considered LT h.ndark+h.nstar)
    if star_ind[0] NE -1 then stat[i].n_star = n_elements(star_ind) else stat[i].n_star = 0 
    if star_ind[0] NE -1 then stat[i].star_mass = total(smass[particles_considered(star_ind)-h.ndark]) else stat[i].star_mass = 0 
  ENDIF else stat[i].star_mass = 0
    dark_ind = where(particles_considered LT h.ndark)
    if dark_ind[0] NE -1 then begin
	stat[i].n_dark = n_elements(dark_ind)
	stat[i].dark_mass = total(dmass[particles_considered(dark_ind)])
	test2 = where(dmass[particles_considered(dark_ind)] GT mindm)
	if test2[0] NE -1 then stat[i].contam = 'contam' else stat[i].contam = 'clean'
    endif else begin
	stat[i].n_dark = 0
	stat[i].dark_mass = 0
        stat[i].contam = 'clean'
        stat[i].false = 'false'
    endelse	
  ENDFOR
ENDIF ELSE BEGIN
  dmass = double(d.mass[0])
  FOR i=0L, ngrps-1 do begin
    particles_considered = index[startind(amigaind[i]):endind(amigaind[i])]
    stat[i].gas_mass = 0
    stat[i].n_gas = 0
    stat[i].star_mass = 0
    stat[i].n_star = 0
    dark_ind = where(particles_considered LT h.ndark)
    stat[i].n_dark = n_elements(dark_ind) 
    stat[i].dark_mass = total(dmass[particles_considered(dark_ind)]) 
    test2 = where(dmass[particles_considered(dark_ind)] GT mindm)
    if test2[0] NE -1 then stat[i].contam = 'contam' else stat[i].contam = 'clean'
  ENDFOR
ENDELSE

; This next is designed to detect if the main halo has been broken down into smaller 
; false halos at the center, a common problem.
pot_false = where(stat.r_vir lt 5.0 and stat.n_gas ne 0 and stat.n_star ne 0, nfalse)
;pot_false = where(stat.r_vir lt 1.0 and stat.n_dark lt 100, nfalse)
if nfalse ne 0 then begin 
    stat[pot_false].false = 'false' 
    stat[0].false = 'false?'
endif

;Find satellites
test_xc = stat.x_c-stat[0].x_c[0]
test_yc = stat.y_c-stat[0].y_c[0]
test_zc = stat.z_c-stat[0].z_c[0]
test_xc = test_xc[1:n_elements(test_xc)-1]
test_yc = test_yc[1:n_elements(test_xc)-1]
test_zc = test_zc[1:n_elements(test_xc)-1]
distance = ((test_xc^2.+test_yc^2.+test_zc^2.)^0.5)*1000. ;change from Mpc to kpc
sat = where(distance le stat[0].r_vir[0],comp=non)
stat[sat+1].sat = 'yes'
stat[non+1].sat = 'no'
stat[0].sat = 'no'


; Write the .stat file
get_lun, lun
filename = strcompress(tipsyfile[j]+'.amiga.stat')
openw, lun, filename
 printf, lun, format='(A5,2x,A9,2x,A8,2x,A8,2x,A8,2x,A16,2x,A10,2x,A16,2x,A16,2x,A16,2x,A10,2x,A9,2x,A10,2x,A9,2x,A9,2x,A9,2x,A10,2x,A10,2x,A10,2x,A7,2x,A10,2x,A7)', 'Grp', 'N_tot', 'N_gas', 'N_star', 'N_dark', 'Mvir(M_sol)', 'Rvir(kpc)', 'GasMass(M_sol)', 'StarMass(M_sol)', 'DarkMass(M_sol)', 'V_max', 'R@V_max', 'VelDisp', 'Xc', 'Yc', 'Zc', 'VXc', 'VYc', 'VZc', 'Contam', 'Satellite?', 'False?' ;'ID_A'

for i=0L,ngrps-1 do printf, lun, format='(I5,2x,I9,2x,I8,2x,I8,2x,I8,2x,A16,2x,D10.4,2x,A16,2x,A16,2x,A16,2x,D10.5,2x,D9.5,2x,D10.5,2x,D9.5,2x,D9.5,2x,D9.5,2x,D10.5,2x,D10.5,2x,D10.5,2x,A7,2x,A7,2x,A7)', stat[i].grp, stat[i].n_part, stat[i].n_gas, stat[i].n_star, stat[i].n_dark, stat[i].m_vir, stat[i].r_vir , stat[i].gas_mass*munit, stat[i].star_mass*munit, stat[i].dark_mass*munit, stat[i].v_max, stat[i].r_max, stat[i].sig_v, stat[i].x_c, stat[i].y_c, stat[i].z_c, stat[i].vx_c, stat[i].vy_c, stat[i].vz_c, stat[i].contam, stat[i].sat, stat[i].false ;stat[i].amiga_orig_id

close, lun
free_lun, lun

; **************************************************************************
; Now, create the .gtp file
; **************************************************************************

statfile = strcompress(tipsyfile[j]+'.amiga.stat')
readcol, statfile, grp, npart, ngas, nstar, ndark, mvir, rvir, gmass, smass, dmass, vmax, rmax, sigv, xc, yc, zc, vxc, vyc, vzc, format='l,l,l,l,l,d,d,d,d,d,d,d,d,d,d,d,d,d,d'
  
star = replicate({mass:0.0, x:0.0, y:0.0, z:0.0, vx:0.0, vy:0.0, vz:0.0, metals:0.0, tform:0.0, eps:0.0, phi:0.0},n_elements(grp))
  for i=0L,n_elements(grp)-1 do begin
	star[i].mass = (gmass[i]+smass[i]+dmass[i])/munit
	star[i].x = (xc[i]*1000./(lunit))-0.5
	star[i].y = (yc[i]*1000./(lunit))-0.5
	star[i].z = (zc[i]*1000./(lunit))-0.5
	star[i].vx = vxc[i]/vunit
	star[i].vy = vyc[i]/vunit
	star[i].vz = vzc[i]/vunit
	star[i].eps = rvir[i]/lunit
  endfor
  star[*].metals = 0.0  
  star[*].phi = 0.0  
  star[*].tform = h.time
  header = h
  header.n = n_elements(grp)
  header.nstar = n_elements(grp)
  header.ngas = 0
  header.ndark = 0
  outfile = strcompress(tipsyfile[j]+'.amiga.gtp')
  wtipsy, outfile, header, gas, dark, star, /standard

ENDFOR

end


pro shrinkspherecenter,minparticles=minparticles,shrinkfactor=shrinkfactor,doplot=doplot,verbose=verbose,dostop=dostop

; for all halos in each file, find the center using shrinking spheres
; from power et al 2003?
; based on Andrew Pontzen's python code, sort of

; minparticles - stop iterating the shrinking sphere when there are
;                less than minparticles left.  default is 30.

; shrinkfactor - how much to shrink the sphere by in each step.
;                default is 90%.

; doplot keyword will plot particles in each concentric sphere in a
; different color so you can check if it's working right.

; dostop keyword will stop after each halo is plotted so you can
; actually look at it if you want

; verbose keyword prints COM and particle info for each iteration


readcol,'files.list',file,format='a',/silent
nfiles = n_elements(file)

if not keyword_set(minparticles) then minparticles = 30
if not keyword_set(shrinkfactor) then shrinkfactor = 0.9d
if keyword_set(doplot) then colors = 20

for ii=0,nfiles-1 do begin
    print,'*** '+file[ii]+' ****'
    rtipsy,file[ii],h,g,d,s
    grp = read_lon_array(file[ii]+'.amiga.grp')
    halocenterfilename = file[ii]+'.shrinkcenters'
    ; ID the BHs, I don't really care about the other halo centers
;    bh = where(s.tform lt 0.0,nbh)
;    bhgrp = grp[h.ngas+h.ndark+bh]
;    uniquebhgrp = bhgrp[uniq(bhgrp,sort(bhgrp))]
    ; changed my mind, let's do it for all halos.
    uniquegrp = grp[uniq(grp,sort(grp))]
    ; ignore halo 0
    notzero = where(uniquegrp ne 0,nnotzero)
    if nnotzero eq 0 then continue
    uniquegrp = uniquegrp[notzero]
    nhalos = n_elements(uniquegrp)
    centerx = fltarr(nhalos)
    centery = fltarr(nhalos)
    centerz = fltarr(nhalos)

    for i=0,nhalos -1 do begin
        if keyword_set(verbose) then print,'*** HALO',uniquegrp[i],' *****'
        ; find particles in this halo
        gasind = where(grp[0:h.ngas-1] eq uniquegrp[i],ngas)
        if ngas ne 0 then gas = g[gasind]
        darkind = where(grp[h.ngas:h.ngas+h.ndark-1] eq uniquegrp[i],ndark)
        if ndark ne 0 then dark = d[darkind]
        starind = where(grp[h.ngas+h.ndark:h.n-1] eq uniquegrp[i],nstar)
        if nstar ne 0 then star = s[starind]
        if ngas ne 0 AND ndark ne 0 AND nstar ne 0 then begin
            masses = double([gas.mass,dark.mass,star.mass])
            x = double([gas.x,dark.x,star.x])
            y = double([gas.y,dark.y,star.y])
            z = double([gas.z,dark.z,star.z])
        endif ; some tiny halos may not have any gas.  but they may have DM and BH
        if ngas eq 0 AND ndark ne 0 AND nstar ne 0 then begin
            masses = [dark.mass,star.mass]
            x = [dark.x,star.x]
            y = [dark.y,star.y]
            z = [dark.z,star.z]
        endif 
        if ngas ne 0 AND ndark ne 0 AND nstar eq 0 then begin
            masses = [gas.mass,dark.mass]
            x = [gas.x,dark.x]
            y = [gas.y,dark.y]
            z = [gas.z,dark.z]
        endif 

        cx = total(masses*x)/total(masses)
        cy = total(masses*y)/total(masses)
        cz = total(masses*z)/total(masses)
        radius = sqrt((x-cx)^2. + (y-cy)^2. + (z-cz)^2.)
        maxrad = max(radius)
        npart = n_elements(x)
        if keyword_set(verbose) then print,'COM:',cx,cy,cz
        if keyword_set(verbose) then print,'radius: ',maxrad,' nparticles: ',npart
        if keyword_set(verbose) then print,'ngas:',ngas,' ndark:',ndark,' nstar:' ,nstar
        if keyword_set(verbose) then print,'**********************'
        if keyword_set(doplot) then plot,x,y,psym=3,/ynozero
        ; now start shrinking the sphere
        while npart gt minparticles do begin
            maxrad = maxrad*shrinkfactor
            if ngas ne 0 then gasradius = radius[0:ngas-1] $
            else gasradius = maxrad*2. ; if no gas then ignore gas.
            if nstar ne 0 then starradius = radius[ngas+ndark:ngas+nstar+ndark-1] $
            else starradius = maxrad*2. ; if no star then ignore star.
            if ndark ne 0 then darkradius = radius[ngas:ngas+ndark-1] $
            else darkradius = maxrad*2.; if no dark then ignore dark
            ; reset radii, particle identities
            oldradius = radius
            newind = where(radius lt maxrad,npart)
	    if npart eq 0 then continue ; exit loop and record info
            masses = masses[newind]
            radius = radius[newind]
            newgasind = where(gasradius lt maxrad,ngas)
            newdarkind = where(darkradius lt maxrad,ndark)
            newstarind = where(starradius lt maxrad,nstar)
            if ngas ne 0 AND nstar ne 0 AND ndark ne 0 then begin
                gas = gas[newgasind]
                dark = dark[newdarkind]
                star = star[newstarind]
                x = [gas.x,dark.x,star.x]
                y = [gas.y,dark.y,star.y]
                z = [gas.z,dark.z,star.z]
            endif
            ; if center has no gas particles
            if ngas eq 0 AND nstar ne 0 AND ndark ne 0 then begin
                dark = dark[newdarkind]
                star = star[newstarind]
                x = [dark.x,star.x]
                y = [dark.y,star.y]
                z = [dark.z,star.z]
            endif
            ; if center has no star particles
            if nstar eq 0 AND ngas ne 0 AND ndark ne 0 then begin
                gas = gas[newgasind]
                dark = dark[newdarkind]
                x = [gas.x,dark.x]
                y = [gas.y,dark.y]
                z = [gas.z,dark.z]
            endif
            ; if center has no dark particles
            if nstar ne 0 AND ngas ne 0 AND ndark eq 0 then begin
                gas = gas[newgasind]
                star = star[newstarind]
                x = [gas.x,star.x]
                y = [gas.y,star.y]
                z = [gas.z,star.z]
            endif
            ; if center has only dark particles
            if nstar eq 0 AND ngas eq 0 then begin
                dark = dark[newdarkind]
                x = dark.x
                y = dark.y
                z = dark.z           
            endif
            ; if center has only star particles
            if ngas eq 0 AND ndark eq 0 then begin
                star = star[newstarind]
                x = star.x
                y = star.y
                z = star.z      
            endif
            ; if center has only gas particles
            if nstar eq 0 AND ndark eq 0 then begin
                gas = gas[newgasind]
                x = gas.x
                y = gas.y
                z = gas.z      
            endif

            if n_elements(x) ne npart then print,'Problem!'
            if n_elements(x) ne npart then stop
            if keyword_set(doplot) then oplot,x,y,psym=3,color=colors
            ; new COM
            cx = total(masses*x)/total(masses)
            cy = total(masses*y)/total(masses)
            cz = total(masses*z)/total(masses)
            ; new radius
            radius = sqrt((x-cx)^2. + (y-cy)^2. + (z-cz)^2.)
            maxrad = max(radius)
            if keyword_set(verbose) then print,'COM:',cx,cy,cz
            if keyword_set(verbose) then print,'radius: ',maxrad,' nparticles: ',npart
            if keyword_set(verbose) then print,'ngas:',ngas,' ndark:',ndark,' nstar:' ,nstar
            if keyword_set(verbose) then print,'**********************'
            if keyword_set(doplot) then colors = colors+10
        endwhile

        ; record final data for each halo
        centerx[i] = cx
        centery[i] = cy
        centerz[i] = cz
        if keyword_set(dostop) then  stop
    endfor

    openw,lun,halocenterfilename,/get_lun
    printf,lun,"; haloid  x  y  z"
    for j=0,nhalos-1 do printf,lun,uniquegrp[j],centerx[j],centery[j],centerz[j]
    close,lun
    free_lun,lun

endfor

end

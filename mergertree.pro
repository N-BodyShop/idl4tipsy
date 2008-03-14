pro mergertree, simname, gid, zs=zs, grps=grps, tipsys=tipsys, maxz=maxz, $
	dKpcUnit=dKpcUnit, halfmassz=halfmassz

  omegam = 0.24
  circlesize = 4.5
  if(keyword_set(dKpcUnit) eq 0) then dKpcUnit = 50e3/0.73
  ;You have: 3 (73 km / s / Mpc)^2 / (8 pi G)
  ;You want: sunmass / (kpc)^3
  ;        * 147.84795  (critical density in above units)
  dMsolUnit = 147.84795*(dKpcUnit)^3.
  standardmass = 1e12/dMsolUnit
  if(keyword_set(maxz) eq 0) then maxz = 4

  if(n_elements(grps) LE 0) then begin
    files = file_search(simname+".0????.data.grp")
    files = files[sort(files)]
    fileind = where((indgen(n_elements(files)) MOD 2) EQ 0)
    files = files[fileind]
    nfs = n_elements(files)
    print,"Reading ",files[nfs-1]
    newgroups = rbarray(files[nfs-1],type='long')
    grps = lonarr(n_elements(newgroups),nfs)
    grps[*,nfs-1] = newgroups
    readgrps = 1
  endif else begin
    nfs = n_elements(grps[0,*])
    newgroups = grps[*,nfs-1]
    readgrps = 0
  endelse

  if(n_elements(zs) LE 0) then begin
    tipsys = file_search(simname+".0????")
    tipsys = tipsys[sort(tipsys)]
    tipsys = tipsys[fileind]
    rtipsy,tipsys[nfs-1],h,/justhead
    z = (1.-h.time)/h.time
    zs = fltarr(nfs)
    zs[0] = z
    pmass = omegam / h.n
    readzs = 1
  endif else begin
    z = zs[nfs-1]
    pmass = omegam / n_elements(grps[*,0])
    readzs = 0
  endelse
  
  ; particles that are member of the group we're concerned with
  ind = where(newgroups EQ gid, nind)
  newgroupmass = nind*pmass

  openw,lun,'g'+strtrim(gid,2)+tipsys[nfs-1]+'.txt',/get_lun
  paperps,file='g'+strtrim(gid,2)+tipsys[nfs-1]+'.eps',charsize=2
  plotsym,0,circlesize * newgroupmass / standardmass
  plot,[z+7e-2],[50],psym=8,xtit='z',yrange=[-20,200], $
    xrange=[1.5e-2*maxz,1.1*maxz], $
    xstyle=1,xmargin=[4.2,0.8],ymargin=[3,0.2],/xlog
  halfmassyet = 0
  for i=nfs-2,0,-1 do begin
    if (readgrps GT 0) then begin
      print,"Reading ",files[i]
      oldgroups = rbarray(files[i],type='long')
      grps[*,i] = oldgroups
    endif else oldgroups = grps[*,i]
    if (readzs GT 0) then begin
      rtipsy,tipsys[i],h,/justhead
      z = (1.-h.time)/h.time
      zs[i] = z
    endif else z = zs[i]
    ; the group ids at this time that contain particles
    ; that are part of gid
    oldnewgroups = oldgroups[ind]
    oldgids = oldnewgroups[uniq(oldnewgroups, sort(oldnewgroups))]
    printf,lun,"z:",z,"   file:  ",tipsys[i]
    printf,lun,"oldgids:", oldgids
    neogids = n_elements(oldgids)
    gchildmass = fltarr(neogids)
    gtotmass = fltarr(neogids)
    height = 0.
    for j=0,neogids-1 do begin
      ; all the particles that are part of oldgids[j]
      oldind = where(oldgroups EQ oldgids[j],noind)
      ; the subset oldgids[j] that will become part of gid
      match, ind, oldind, imatch, oimatch
      childgrps = oldind[oimatch]
      ; particles that aren't part of gid, but are part of the oldgid
      ;strippedois = 
      gchildmass[j] = n_elements(childgrps)*pmass
      gtotmass[j] = noind*pmass
    endfor
    ; Plot in mass order.  Group 0 has to be at bottom, then from
    ; most to least massive going up.
    ;0 group plotting
    zeromass = gchildmass[where(oldgids EQ 0,complement=nozero)]
    plotsym,0,circlesize*zeromass[0]/standardmass
    oplot,[z],[0],psym=8
    height = height+20

    if(nozero[0] NE -1) then begin
      ; sort the rest by mass
      spind = nozero[reverse(sort(gchildmass[nozero]))]
      if(halfmassyet eq 0 AND $
         gchildmass[spind[0]]/newgroupmass LT 0.5) then begin
	    halfmassz = z
	    halfmassyet = 1
      endif
      for j=0,neogids-2 do begin
        plotsym,0,circlesize *gchildmass[spind[j]]/ standardmass
        oplot,[z],[height],psym=8
        xyouts,[0.9*z],[height-10],strtrim(oldgids[spind[j]],2),charsize=1
        height = height + 20
      endfor
    endif
    ; print out arrays hopefully in the same order as oldgids are printed
    printf,lun,"group tot mass: ",gtotmass
    printf,lun,"mass to be in gid: ",gchildmass
  endfor
  close,lun
  free_lun,lun

end

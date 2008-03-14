function calcidlstats,file

idlstatfile = file+'.idl.stat'
period = 1.
rmaxneighbor = 0.04 ; ~2-3 Mpc in 50 h^-1 Mpc box

print,'reading ',file
rtipsy,file,h,g,d,s
print,'reading ',file+'.data.grp'
grp = rbarray(file+'.data.grp',type='long')
ngrps = long(max(grp))
gmass = fltarr(ngrps)
glambda = fltarr(ngrps)
gid = fltarr(ngrps)
record={gid:0L,mass:0.0,lambda:0.0,maxr:0.0,maxneighbormassr:0.0,$
        comx:0.0,comy:0.0,comz:0.0,comvx:0.0,comvy:0.0,comvz:0.0, $
        jx:0.0,jy:0.0,jz:0.0,jtot:0.0,nneighbors:0,maxneighbormass:0.0}
s = replicate(record, ngrps)
print,'calculating stats'
for i=1L,ngrps do begin
  s[i-1].gid = i
  if (i MOD 100 EQ 0) then print,i
  ind = where(grp eq i)
  ps = d[ind]
  mass = total(ps.mass)
  s[i-1].mass = mass*4.7526e16
  ; group reference point
  relx = ps[0].x
  rely = ps[0].y
  relz = ps[0].z
  ; center to reference point frame
  posind = where(ps.x - relx GT 0.5*period, nposind)
  if(nposind GT 0) then ps[posind].x = ps[posind].x - period
  negind = where(ps.x - relx LE 0.5*period, nnegind)
  if(nnegind GT 0) then ps[negind].x = ps[negind].x - period
  s[i-1].comx = total(ps.x*ps.mass)/mass
  ;y's
  posind = where(ps.y - rely GT 0.5*period, nposind)
  if(nposind GT 0) then ps[posind].y = ps[posind].y - period
  posind = where(ps.y - rely LE 0.5*period, nnegind)
  if(nnegind GT 0) then ps[negind].y = ps[negind].y - period
  s[i-1].comy = total(ps.y*ps.mass)/mass
  ;z's
  posind = where(ps.z - relz GT 0.5*period, nposind)
  if(nposind GT 0) then ps[posind].z = ps[posind].z - period
  posind = where(ps.z - relz LE -0.5*period, nnegind)
  if(nnegind GT 0) then ps[negind].z = ps[negind].z + period
  s[i-1].comz = total(ps.z*ps.mass)/mass

  s[i-1].comvx = total(ps.vx*ps.mass)/mass
  s[i-1].comvy = total(ps.vy*ps.mass)/mass
  s[i-1].comvz = total(ps.vz*ps.mass)/mass

  x = ps.x -s[i-1].comx
  y = ps.y -s[i-1].comy
  z = ps.z -s[i-1].comz
  r = sqrt(x*x + y*y + z*z)
  s[i-1].maxr = max(r)
  vx = ps.vx -s[i-1].comvx
  vy = ps.vy -s[i-1].comvy
  vz = ps.vz -s[i-1].comvz

  jx = total(ps.mass * (y*vz - z*vy))
  jy = total(ps.mass * (z*vx - x*vz))
  jz = total(ps.mass * (x*vy - y*vx))

  jtot = sqrt(jx*jx + jy*jy + jz*jz)
  s[i-1].jx = jx
  s[i-1].jy = jy
  s[i-1].jz = jz
  s[i-1].lambda = jtot / sqrt(5./3. * mass*mass*mass * max(r))

  if(s[i-1].comx GT 0.5*period) then s[i-1].comx = s[i-1].comx - period
  if(s[i-1].comx LE -0.5*period) then s[i-1].comx = s[i-1].comx + period
  if(s[i-1].comy GT 0.5*period) then s[i-1].comy = s[i-1].comy - period
  if(s[i-1].comy LE -0.5*period) then s[i-1].comy = s[i-1].comy + period
  if(s[i-1].comz GT 0.5*period) then s[i-1].comz = s[i-1].comz - period
  if(s[i-1].comz LE -0.5*period) then s[i-1].comz = s[i-1].comz + period
endfor

print,'Finding neighbors'
for i=1L,ngrps do begin
  if (i MOD 1000 EQ 0) then print,i
  dx = s.comx - s[i-1].comx
  dy = s.comy - s[i-1].comy
  dz = s.comz - s[i-1].comz
  posind = where(dx GT 0.5*period,nposind)
  if(nposind GT 0) then dx[posind] = dx[posind] - period
  negind = where(dx GT 0.5*period,nnegind)
  if(nnegind GT 0) then dx[negind] = dx[negind] + period
  posind = where(dy GT 0.5*period,nposind)
  if(nposind GT 0) then dy[posind] = dy[posind] - period
  negind = where(dy GT 0.5*period,nnegind)
  if(nnegind GT 0) then dy[negind] = dy[negind] + period
  posind = where(dz GT 0.5*period,nposind)
  if(nposind GT 0) then dz[posind] = dz[posind] - period
  negind = where(dz GT 0.5*period,nnegind)
  if(nnegind GT 0) then dz[negind] = dz[negind] + period
  r = sqrt(dx*dx + dy*dy + dz*dz)
  ind = where(r LT rmaxneighbor AND r GT 0,nind)
  s[i-1].nneighbors = nind
  if(s[i-1].nneighbors GT 0) then begin
    s[i-1].maxneighbormass = max(s[ind].mass)
    mind = min(ind[where(s[ind].mass eq max(s[ind].mass))])
    s[i-1].maxneighbormassr = r[mind]
  endif
endfor
print,'writing .idl.stat file'
openw,lun,file+'.idl.stat',/get_lun,/xdr
writeu,lun,s;i,mass*4.7526e16,lambda,max(r),comx,comy,comz,comvx,comvy,comvz,jx,jy,jz,jtot
close,lun
free_lun,lun

print,'plotting mass, lambda'
paperps,file=file+'.mass.eps',charsize=2
  hm = histogram(s.mass,bin=min(s.mass),locations=massbins)
  plot,massbins,hm,/xlog,xrange=[min(s.mass),max(s.mass)], $
	  /ylog,yrange=[0.9,max(hm)],psym=10, $
	  xmargin=[6.0,0.8],ymargin=[3.0,0.5], $
	  xtit='Mass [M'+sunsymbol()+']',ytit = 'N / dM'
ppsclose

paperps,file=file+'.lambda.eps',charsize=2
  hl = histogram(s.lambda,bin=min(s.lambda),locations=lambdabins)
  plot,lambdabins,hl,/xlog,xrange=[min(s.lambda),max(s.lambda)], $
	  /ylog,yrange=[0.9,max(hl)],psym=10, $
	  xmargin=[6.0,0.8],ymargin=[3.0,0.5], $
	  xtit=textoidl(' \lambda'),ytit = 'N / d'+textoidl(' \lambda')
ppsclose

return,s;{gid:gid,mass:gmass,lambda:glambda}

END

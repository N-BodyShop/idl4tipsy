function apocenters,g,d,s

; Finds apocenters of just STARS! from what's passed in.
;
; g is the gas particles
; d is the dark particles
; s is the star particles

print,"Finding stellar apocenters"
   phis = [g.phi,d.phi,s.phi]
   allrs= [sqrt(g.x^2. + g.y^2. +g.z^2.), $
            sqrt(d.x^2. + d.y^2. +d.z^2.), $
            sqrt(s.x^2. + s.y^2. +s.z^2.)]
   rsqr = s.x^2. + s.y^2. +s.z^2.
   rs = sqrt(rsqr)
   vsqr = (s.vx*s.x + s.vy*s.y +s.vz*s.z)/rs
   vs = sqrt(vsqr)

   nea = n_elements(allrs)-1
   isr = sort(allrs)
   apo = fltarr(n_elements(s))
   for i =0L,n_elements(s)-1 do begin
    Larr =crossp([s[i].x,s[i].y,s[i].z],[s[i].vx,s[i].vy,s[i].vz])
    Lsqr=transpose(Larr)#Larr
    phi_effi = s[i].phi + Lsqr[0]/(2.0*rsqr[i])
    goal = 0.5*vs[i]^2. + phi_effi
    temp = where(allrs[isr] EQ rs[i])
    lowi = temp[0]
    highi = nea
    if(2*lowi LT highi) then highi = 2*lowi
    loval = phis[isr[lowi]] + Lsqr[0]/(2.*allrs[isr[lowi]]^2.)
    hival = phis[isr[highi]] + Lsqr[0]/(2.*allrs[isr[highi]]^2.)
   addstep = highi - lowi
    while(hival LT goal) do begin
      if ( highi + addstep GE nea ) then begin
        highi = nea - 1
        break
      endif else highi = highi + addstep ;$
      hival = phis[isr[highi]] + Lsqr[0]/(2.*allrs[isr[highi]]^2.)
    endwhile
      ;if(2*highi LT nea) then highi = 2*highi $
      ;else highi =nea
    halfi = floor((highi-lowi)/2.)
    newval = phis[isr[halfi]] + Lsqr[0]/(2.*allrs[isr[halfi]]^2.)
    frdiff =abs((newval - goal)/goal)
    its = 0
    while (frdiff GT 1e-3) do begin
     if ( hival LT loval) then begin ; should never do this part
      if ( newval GT goal ) then lowi = halfi $
      else highi = halfi
     endif else begin
      if ( newval GT goal ) then highi = halfi $
      else lowi = halfi
     endelse
      halfi = floor((highi-lowi)/2)
      halfi = lowi + halfi
      newval = phis[isr[halfi]] + Lsqr[0]/(2.*allrs[isr[halfi]]^2.)
      frdiff =abs((allrs[isr[halfi]] - allrs[isr[lowi]])/allrs[isr[halfi]])
     if(its GT 17) then begin
      if (highi - halfi LT 3) then continue
      print,'Not converging'
      STOP
     endif
     its = its+1
    endwhile
    apo[i] = allrs[isr[halfi]]
; Keplerian system only
    ;GM = -s[i].phi * rs[i]
    ;e = Lsqr[0]/GM
    ;a = Lsqr[0]/GM/(1.-e^2.)
    ;apo[i] = a*(1+e)
   endfor

return,apo

end

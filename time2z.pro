function z2Hub, z
  a=1.0/(1.0+z)
  H = 70.0 * 1.02528 * 10.0^(-12.0) ;in years^-1
  omega_m=0.3
  omega_v=0.7
  curvature=1-omega_m-omega_v
  return,H*sqrt(omega_m*a+curvature*a*a+omega_v*a*a*a*a)/(a*a)
end

function time2z,t
  ageofuniverse=1.346773e10
  lbt=ageofuniverse-t
  znew = (lbt/ageofuniverse)^(2.0/3.0)
  zold = fltarr(n_elements(lbt))
  for i=0,n_elements(lbt)-1 do begin
    it=0
    while (abs(znew[i]-zold[i])/znew[i] GT 1e-7) do begin
      f = lbt[i] - lookback(znew[i])
      fprime = 1.0/(znew[i]*z2Hub(znew[i]))
      zold[i]=znew[i]
      znew[i] = znew[i] + f/fprime
      it=it+1
      if ((it GT 20) OR (znew[i] GT 99)) then break
    endwhile
  endfor

  return,znew
end

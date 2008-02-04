pro plots_xvy, x, y, z, $
               _extra=_extra, zrange=zrange, aplotkey=aplotkey

  if n_params() eq 0 then begin
    print,"PLOTS_XVY, x, y, z, zrange=zrange, _extra=_extra"
    print,"   Makes an x versus y plot with symbols color-coded"
    print,"   according to value of z (via range 'range')."
    return
  endif

  nel = n_elements(redshift)
  if n_elements(magnitude) ne nel or n_elements(param3) ne nel then begin
    print,"PLOTS_XVY: size(x) = size(y) = size(z) please!"
    return
  endif

  if not keyword_set(zrange) then begin
    zrange=[min(z),max(z)]
  endif

  if n_elements(aplotkey) eq 0 then aplotkey=1
  
  ncol = 255
  if aplotkey then aplot, 1, x, y, /nodata, _extra=_extra $
  else plot,x,y,/nodata,_extra=_extra
  plots, x, y, color=(round(interpol([0,ncol-1],zrange,z))<ncol-1>0),noclip=0,_extra=_extra

end

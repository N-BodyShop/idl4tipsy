; PRO MOMENTS
;
; Creates moment maps from a data cube. Use read_cube_fits to read in
; a fits file and corresponding header structure.
;
; INPUTS:
;
;    CUBE: An array of the 3d HI data cube (read in by read_cube_fits)
;
;    HEADER: The header structure (as read in by read_cube_fits)
;
;
; OUTPUTS:
;
;    MOM0: Output 2d array of the mom0 map, in atoms/cm^2
;
;    MOM1: Output 2d array of the moment 1 map (km/s)
;
;    MOM2: Output 2d array of the moment 2 map (km/s)
;
;    XAXIS: The scaled xaxis
;
;    YAXIS: The scaled yaxis
;
;
; OPTIONAL INPUTS:
;
;    PLOT: Set this in order to plot to the screen. If both PLOT and
;          OUTFILE are set, it defaults to PLOT (plots to screen)
;
;    THRESHOLD: A sensitivity threshold for column density in
;               atoms/cm^2. If set, all moment maps are evaluated
;               only above the cutoff threshold given.
; 
;    OUTFILE: Set to plot to a file. If both PLOT and OUTFILE are set,
;             it defaults to PLOT (plots to screen)
;
;    SELECT: Doesn't do anything yet. I'm thinking about implementing
;            a clickable interface in order to select region of
;            interest.
;
;    NL0: Number of contour levels for the column density (mom0)
;         map. Defaults to 10.
;
;    NL1: Number of contour levels for the velocity field (mom1)
;         map. Defaults to 10.
;
;    NL2: Number of contour levels for the velocity dispersion (mom2)
;         map. Defaults to 10.
;
;    NAME: A name to append on the plot title.
;
;
;
; To go from atoms/cm^2 use:
; c = constants()
; foo = mom0 * c.mh / c.msun * (c.kpc / 1000.)^2
;
pro moments, cube, header, mom0, mom1, mom2, xaxis=xaxis, $
             yaxis=yaxis, $
             doplot=doplot, threshold=threshold, outfile=outfile,$
             select=select, plotthresh=plotthresh,nl0=nl0, nl1=nl1,$
             nl2=nl2, name=name, maxthreshold=maxthreshold, $
             nc0=nc0, nc1=nc1, nc2=nc2, docontour=docontour, fits=fits,$
             xrange=xrange, yrange=yrange, bigiel=bigiel

if not (keyword_set(cube) and keyword_set(header)) then begin
    print,"moments, cube, header, [mom0, mom1, mom2, xaxis, yaxis, $"
    print,"         plot=plot, threshold=threshold, outfile=outfile, $"
    print,"         nl0=nl0, nl1=nl1, nl2=nl2, name=name]"
    return
endif
if not keyword_set(cube) then begin
    print,"CUBE must be set to a data cube array."
    return
endif else if n_elements(size(cube,/dimensions)) ne 3 then begin
    print,"CUBE must be a three-dimensional array."
endif
if not keyword_set(header) then begin
    print,"HEADER keyword must be set."
    return
endif
if not keyword_set(nl0) then nl0 = 10.
if not keyword_set(nl1) then nl1 = 10.
if not keyword_set(nl2) then nl2 = 10.
if not keyword_set(threshold) then threshold = 0.
if (nl0 mod 1) ne 0 then begin
    print,"nl0 must be an integer."
    return
endif
if (nl1 mod 1) ne 0 then begin
    print,"nl1 must be an integer."
    return
endif
if (nl2 mod 1) ne 0 then begin
    print,"nl2 must be an integer."
endif
if not keyword_set(threshold) then threshold=1d18

; set up axes
;xaxis = dindgen(header.naxis1) * header.cdelt1 + header.crval1
;yaxis = dindgen(header.naxis2) * header.cdelt2 + header.crval2
;vaxis = dindgen(header.naxis3) * header.cdelt3 + header.crval3

getkpc, header, xaxis=xaxis, yaxis=yaxis, vaxis=vaxis

; sum along velocity axis for first moment.
mom0 = total(cube, 3)

; figure out where mom0 equals 0
zero = where(mom0 le 0, complement=nonzero)



; set up other moments
vel = double(cube)
vel2 = double(cube)
for i=0,n_elements(vaxis)-1 do begin
    vel[*,*,i] = cube[*,*,i] * vaxis[i]
    vel2[*,*,i] = cube[*,*,i] * (vaxis[i])^2
endfor

; sum along velocity axis for moments 1 and 2
mom1 = total(vel, 3)
mom2 = total(vel2, 3)
; normalize mom1
if (zero[0] ne -1) then begin
    mom1[nonzero] = mom1[nonzero] / mom0[nonzero]
    mom2[nonzero] = mom2[nonzero] / mom0[nonzero]
    mom1[zero] = 0.0
    mom2[zero] = 0.0
endif else begin 
    mom1 = mom1 / mom0
    mom2 = mom2 / mom0
endelse

mom2 = sqrt(mom2 - mom1^2)

; scale the moment 0 map to proper units:
; mom0 (atoms/cm^2) = x Msun / boxarea * (1 kpc in cm)^2 * x grams/Msun *
;    hydrogen atom in g
mom0 = double(mom0) / (header.cdelt1 * header.cdelt2) * $
  (1.0d / double(3.08e21))^2 * double(1.989e33) / double(1.67e-24)

; artificially chop off center of moment 0 maps
if keyword_set(maxthreshold) then begin
    high = where(mom0 gt maxthreshold)
    if high[0] ne -1 then mom0[high] = maxthreshold
endif

vals = where(mom0 lt threshold)
if (vals[0] ne -1) then begin
    mom0[vals] = !values.d_nan
    mom1[vals] = !values.d_nan
    mom2[vals] = !values.d_nan
endif


; set up for ranges
if keyword_set(xrange) then begin
    xind = where((xaxis ge xrange[0]) and (xaxis le xrange[1]))
    mom0 = mom0[xind,*]
    mom1 = mom1[xind,*]
    mom2 = mom2[xind,*]
    xaxis = xaxis[xind]
endif

if keyword_set(yrange) then begin
    yind = where((yaxis ge yrange[0]) and (yaxis le yrange[1]))
    mom0 = mom0[*,yind]
    mom1 = mom1[*,yind]
    mom2 = mom2[*,yind]
    yaxis = yaxis[yind]
endif

if keyword_set(fits) then begin
;    hmom1
;    hmom2
;    hmom3
    mwrheader = makeheader(header)
    foo = size(fits, /structure)
    if (foo.type_name eq 'int') then begin
        if (fits eq 1) then begin
            mwrfits, mom0, 'mom0.fits', mwrheader, /create
            mwrfits, mom1, 'mom1.fits', mwrheader, /create
            mwrfits, mom2, 'mom2.fits', mwrheader, /create
        endif
    endif else begin
        mwrfits, mom0, fits + '.mom0.fits', mwrheader, /create
        mwrfits, mom1, fits + '.mom1.fits', mwrheader, /create
        mwrfits, mom2, fits + '.mom2.fits', mwrheader, /create
    endelse
endif

if keyword_set(outfile) or keyword_set(doplot) then begin

    if not keyword_set(doplot) then begin
        set_plot,'ps'
        device,filename=outfile + '.mom0.ps',xsize=7,ysize=5.7,/inches
    endif else begin
        window, 4, xsize=700,ysize=570
    endelse
    loadct,0
    levs0 = dblarr(nl0)
    levs0[nl0-1] = max(mom0[where(finite(mom0))]) * .95
    for i=0, nl0-2 do begin
        levs0[nl0-2-i] = levs0[nl0-1-i] / sqrt(2.)
    endfor
    levs0 = levs0[where(levs0 gt threshold)]
    nl0 = n_elements(levs0)
    cvec = floor((indgen(nl0) / (-1. * nl0) + 1.)*250)
    titl = 'HI Column Density (atoms / cm' + textoidl('^2') + ')'
    if keyword_set(name) then titl = titl + ', ' + name

    contour, mom0, xaxis, yaxis, levels=levs0,xtitle='kpc',ytitle='kpc',$
      title=titl,position=[0.1,0.1,0.75,0.90],/fill, xstyle=1,ystyle=1,$
      c_colors=cvec, xrange=[-10,10],yrange=[-10,10]
    if keyword_set(docontour) then $
      contour, mom0, xaxis, yaxis, levels=levs0, /overplot
    if keyword_set(bigiel) then $
      loadct, 39
      contour, mom0, xaxis, yaxis, levels=1.13216e21, c_colors=250, /overplot
    loadct, 0
    colorbar, divisions=6, ncolors=cvec[0]-cvec[n_elements(cvec)-1], $
      range=[levs0[0],levs0[n_elements(levs0)-1]],$
      /vertical, position=[0.1,0.80,0.90,0.85],format='(%"%7.1e")',$
      /invert,/right;, xrange=[-10,10], yrange=[-10,10]

    if not keyword_set(doplot) then begin 
        device,/close
        device,filename=outfile + '.mom1.ps',/color,bits=8,$
          xsize=7, ysize=5.7, /inches, /encapsulated
    endif else stop

    loadct,39
    cvec = indgen(nl1) / nl1 * (250) + 20
    vmax = max(mom1[where(finite(mom1))])
    vmin = min(mom1[where(finite(mom1))])
    levs1 = dindgen(nl1) / double(nl1) * (vmax - vmin) + vmin
    titl = 'HI Velocity Field (km/s)'
    if keyword_set(name) then titl = titl + ', ' + name

    contour, mom1, xaxis, yaxis, levels=levs1, xtitle='kpc', ytitle='kpc', $
      title=titl, position=[0.1,0.1,0.75,0.90],/cell_fill, xstyle=1,ystyle=1,$
      c_colors=cvec
;    contour, mom1, xaxis, yaxis, levels=levs1, /overplot
    colorbar, /right, /vertical, position=[0.1,0.80,0.90,0.85], $
      format='(%"%7.3g")', range=[levs1[0], levs1[n_elements(levs1)-1]], $
      divisions=6, bottom=cvec[0],ncolors=cvec[n_elements(cvec)-1] - cvec[0]

    if not keyword_set(doplot) then begin
        device,/close
        device,filename=outfile + '.mom2.ps', /color,bits=8,$
          xsize=7, ysize=5.7, /inches, /encapsulated
    endif else stop

    loadct, 39
    cvec = dindgen(nl2) / nl2 * 250 + 20
    vsqmax = max(mom2[where(finite(mom2))])
    vsqmin = min(mom2[where(finite(mom2))])
    levs2 = dindgen(nl2) / double(nl2) * (vsqmax - vsqmin) + vsqmin
    titl = 'HI Velocity Dispersion (km/s)'
    if keyword_set(name) then titl = titl + ', ' + name

    contour, mom2, xaxis, yaxis, levels=levs2, xtitle='kpc', ytitle='kpc', $
      title=titl, position=[0.1,0.1,0.75,0.90],/cell_fill, xstyle=1,ystyle=1,$
      c_colors=cvec
    colorbar, /right, /vertical, position=[0.1,0.80,0.90,0.85], $
      format='(%"%7.3g")', range=[levs2[0], levs2[n_elements(levs1)-1]],$
      divisions=6, bottom=cvec[0],ncolors=cvec[n_elements(cvec)-1] - cvec[0]

    if not keyword_set(doplot) then begin
        device,/close
        set_plot,'x'
    endif else stop


endif

;return,mom0

end


function makeheader, cubehead



sxaddpar, header, 'CTYPE1', 'KPC-LIN'
sxaddpar, header, 'CUNIT1', 'KPC'
sxaddpar, header, 'CRVAL1', cubehead.crval1
sxaddpar, header, 'CDELT1', cubehead.cdelt1
sxaddpar, header, 'CRPIX1', cubehead.crpix1

sxaddpar, header, 'CTYPE2', 'KPC-LIN'
sxaddpar, header, 'CUNIT2', 'KPC'
sxaddpar, header, 'CRVAL2', cubehead.crval2
sxaddpar, header, 'CDELT2', cubehead.cdelt2
sxaddpar, header, 'CRPIX2', cubehead.crpix2

return,header

end

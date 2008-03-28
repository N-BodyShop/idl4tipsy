;fit an ellipse at some user defined surface brightness contour, and output the
;ratio of the semi-minor to semi-major axis
; you need to input filename,the desired surface brightness contour,the band #
; from filters list, the relevant camera angle #,fov and # of pixels
pro ellipse,file,SBlimit,band,cam,fov,npix

image =mrdfits(file,cam); 
gal = image(*,*, band[0])
makex, gal,x,y  ;set up coordinate grid
pixKpc = fov/npix

; the next bit changes to units of M/arcsec^2
filters= mrdfits(file[0],12) ; the 12 needs adjusting depending on no.cameras
Lunit = filters[band].L_lambda_to_L_n ; internal units to W/m^2
Fo=3.631e-23
units = 4.35e10 ; sr>>> arcsec^2
sbfactor = Lunit/units/Fo
gal=-2.5*alog10(sbfactor*abs(gal))

; now define the region of interest and use fit_ellipse
ind=where(gal LT SBlimit)
ell=gal*0
ell[ind]=1  
data=fit_ellipse(ind,xsize=npix,ysize=npix,SemiAxes=semiaxes)
;show data & best-fit ellipse
tvscl,ell
plots,data[0,*],data[1,*],color=128,/device,thick=2
print,'ellipticity =',10*(semiaxes[0]-semiaxes[1])/semiaxes[0]
stop
end


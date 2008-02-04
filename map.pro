;make a surface brightness map with color bar in your favourite band
pro map
file = 'directory/broadband.fits' ; tell it where to look!

cam=16 ; face on with dust needs adjusting 
band = 3 ; chose your band (look in the  filters file)
SBmin = 31 ; max and min surface brightness levels
SBmax = 23

begplot,name="SBmaps.ps",/color,xsize=5.,ysize=5 ; size is inches

image =mrdfits(file,cam); 
gal = image(*,*, band)
gal-abs(gal)
makex, gal,x,y  ;set up coordinate grid; this makex was written by Peter Yoachimm and is also "checked in"
pixKpc = 25./600.  ; fov/npixels needs to be adjusted acordingly
filters= mrdfits(file,12) ; the 12 needs adjusting depending on no.cameras
Lunit = filters[band].L_lambda_to_L_n ; internal units to W/m^2
Fo=3.631e-23
units = 4.35e10 ; sr>>> arcsec^2
sbfactor = Lunit/units/Fo
gal=sbfactor*gal
x=pixKpc*x
y=pixKpc*y

;; define a square as your plotting symbol
usersym, [-.5,-.5,.5,.5], [-.5,.5,.5,-.5], /fill
symsize = symsize

;; you may need to get plotting code: plots_xvy.pro
loadct,3
plots_xvy, x, y, -2.5*alog10(gal), zrange=[SBmin,SBmax],psym=8,symsize=symsize,xrange=[-15,15],yrange=[-15,15],/xstyle,/ystyle,position=[0.15,0.2,0.9,0.95],aplot=0,charsize=0.7,xtit="kpc",ytit="kpc"

loadct,3
colorbar,range=[SBmin,SBmax],color=0,charsize=0.75,position=[0.2,0.05,0.5,0.15],format='(I2.0)',divisions=(SBmin-SBmax),xtitle=textoidl("M arcsec^{-2}"),/noerase
endplot
stop
end


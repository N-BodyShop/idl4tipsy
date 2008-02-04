pro makex, im, x,y ;take an image and make axis for models
xn=n_elements(im(*,0))
yn=n_elements(im(0,*))
xr=dindgen(xn)-xn/2.+.5
yc=dindgen(yn)-yn/2.+.5
x=xr#(yc*0+1)
y=(xr*0+1)#yc

return
end
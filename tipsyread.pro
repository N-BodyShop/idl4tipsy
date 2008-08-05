pro tipsyread

;read the output from tipsy's "writebox" into a data structure
; and saves it as a fits file.

;if (N_params() lt 1) then begin
;    message, 'syntax: type tipsyread,file=filename ', /INFO
;    return
;endif

;file='168.box'
readcol,'boxes.list',boxes,format='a',/silent

for i=0,n_elements(boxes)-1 do begin

file=boxes[i]

a=read_ascii(file,count=count,record_start=3,delimiter=',') ; everything else
readcol,file,numline=3,format='l,l,l',one,two,three,/silent ; #s of particles
readcol,file,numline=2,skipline=1,ttime,/silent ;ndim and time

nparticles=one[0]
ngas=two[0]
nstar=three[0]
ndim=ttime[0]
time=ttime[1]
ndark=nparticles-nstar-ngas


;stop

marker=7*nparticles+ndark
marker2=marker+nstar
marker3=marker2+4*ngas

;darksoft is taken out since there is no dark in the tube!!!

mass=a.field1[0:nparticles-1]
x=a.field1[nparticles:(2*nparticles)-1]
y=a.field1[(2*nparticles):(3*nparticles-1)]
z=a.field1[(3*nparticles):(4*nparticles-1)]
vx=a.field1[(4*nparticles):(5*nparticles-1)]
vy=a.field1[(5*nparticles):(6*nparticles-1)]
vz=a.field1[(6*nparticles):(7*nparticles-1)]
;darksoft=a.field1[(7*nparticles):(marker-1)]
starsoft=a.field1[(marker):(marker2-1)]
den=a.field1[marker2:(marker2+ngas)-1]
temp=a.field1[(marker2+ngas):(marker2+(2*ngas))-1]
smooth=a.field1[(marker2+(2*ngas)):(marker2+(3*ngas))-1]
gasmetal=a.field1[(marker2+(3*ngas)):(marker2+(4*ngas))-1]
starmetal=a.field1[marker3:marker3+nstar-1]
tform=a.field1[marker3+nstar:(marker3+(2*nstar))-1]
pe=a.field1[marker3+(2*nstar):(marker3+(2*nstar)+nparticles)-1]

;stop

data={mass:mass,x:x,y:y,z:z,vx:vx,vy:vy,vz:vz,starsoft:starsoft,den:den,temp:temp,smooth:smooth,gasmetal:gasmetal,starmetal:starmetal,tform:tform,pe:pe}

;stop
mwrfits,data,file+strtrim('.fits'),/create
;stop
endfor


return
end


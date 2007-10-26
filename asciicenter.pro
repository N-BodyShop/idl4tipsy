pro asciicenter, filename, h=h,g=g,d=d,s=s,dist_units=dist_units

if(keyword_set(dist_units) EQ 0) then dist_units = 1.
if(keyword_set(h) EQ 0 AND keyword_set(s) EQ 0) then rtipsy, filename, h,g,d,s


cx = total(s.x*s.mass)
cy = total(s.y*s.mass)
cz = total(s.z*s.mass)
mt = total(s.mass)

cx= cx/mt
cy= cy/mt
cz= cz/mt

g.x = g.x - cx
g.y = g.y - cy
g.z = g.z - cz
d.x = d.x - cx
d.y = d.y - cy
d.z = d.z - cz

s.x =  s.x - cx
s.y =  s.y - cy
s.z =  s.z - cz


;***********************************************
;eliminate proper motion of the halo/galaxy
;***********************************************
propx=mean(s.vx)
propy=mean(s.vy)
propz=mean(s.vz)

s.vx=s.vx-propx
s.vy=s.vy-propy
s.vz=s.vz-propz
g.vx=g.vx-propx
g.vy=g.vy-propy
g.vz=g.vz-propz
d.vx=d.vx-propx
d.vy=d.vy-propy
d.vz=d.vz-propz

radius = dist_units*sqrt(s.x^2+s.y^2+s.z^2)
limit = max(radius)/dist_units ; kpc
gout = where (sqrt(g.x*g.x + g.y*g.y +g.z*g.z) lt limit)
dout = where (sqrt(d.x*d.x + d.y*d.y +d.z*d.z) lt limit)


num = h 
num.nstar = n_elements(s)
num.ngas = n_elements(gout)
num.ndark = n_elements(dout)
num.n = num.nstar + num.ngas + num.ndark

print,"Creating "+filename+".asc"
wtipsy_asc,filename+'.asc',num,g[gout],d[dout],s

end


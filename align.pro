pro align,gas,dark,stars,limit

;***************************************************************
; align the particles so angular momentun vector parallel to z axis
; limit sets the radius of within which you use gas to align 
;***************************************************************

J=0.
Jxs=0.
Jys=0.
Jzs=0.
Js=0.
mt=0.

;*****************************
;if you want to use young stars
;******************************
;r = where(sqrt(stars.x*stars.x+stars.y*stars.y+stars.z*stars.z) lt limit and stars.tform gt 0.25)
;mt = total(stars[r].mass)
;Jxs = total(stars[r].mass*(stars[r].y*stars[r].vz-stars[r].z*stars[r].vy))
;Jys = total(stars[r].mass*(stars[r].z*stars[r].vx-stars[r].x*stars[r].vz))
;Jzs = total(stars[r].mass*(stars[r].x*stars[r].vy-stars[r].y*stars[r].vx))

;****; if you want to use gas instead of stars
;r = where(sqrt(gas.x*gas.x+gas.y*gas.y+gas.z*gas.z) gt limit/4. and
;sqrt(gas.x*gas.x+gas.y*gas.y+gas.z*gas.z)  lt limit)
r = where(sqrt(gas.x*gas.x+gas.y*gas.y+gas.z*gas.z)  lt limit)
mt = total(gas[r].mass)
Jxs = total(gas[r].mass*(gas[r].y*gas[r].vz-gas[r].z*gas[r].vy))
Jys = total(gas[r].mass*(gas[r].z*gas[r].vx-gas[r].x*gas[r].vz))
Jzs = total(gas[r].mass*(gas[r].x*gas[r].vy-gas[r].y*gas[r].vx))
;*****************************

Js=sqrt(Jxs*Jxs+Jys*Jys+Jzs*Jzs)

if Js gt 0. then begin
           jjx=Jxs/Js
           jjy=Jys/Js
           jjz=Jzs/Js
           costh=jjz
           sinth=sqrt(1.0-jjz*jjz)
if  sinth gt 0.0 then begin
              sinph=jjy/sinth
              cosph=jjx/sinth
endif
endif 
if Js le 0.  then begin
           cosph = 1.0
           sinph = 0.0
endif

        ax=costh*cosph
        bx=costh*sinph
        cx=-sinth
        ay=-sinph
        by=cosph
        cy=0.0
        az=sinth*cosph
        bz=sinth*sinph
        cz=costh


        print, ax,bx,cx
        print, ay,by,cy
        print, az,bz,cz
; /**** translate star particles (and change units) ****/

           txs=stars.x
           tys=stars.y
           tzs=stars.z
           stars.x=(ax*txs+bx*tys+cx*tzs)
           stars.y=(ay*txs+by*tys+cy*tzs)
           stars.z=(az*txs+bz*tys+cz*tzs)

           txs=stars.vx
           tys=stars.vy
           tzs=stars.vz
           stars.vx=(ax*txs+bx*tys+cx*tzs)
           stars.vy=(ay*txs+by*tys+cy*tzs)
           stars.vz=(az*txs+bz*tys+cz*tzs)
; /**** translate gas particles  (and change units) ****/

           txs=gas.x
           tys=gas.y
           tzs=gas.z
           gas.x=(ax*txs+bx*tys+cx*tzs)
           gas.y=(ay*txs+by*tys+cy*tzs)
           gas.z=(az*txs+bz*tys+cz*tzs)

           txs=gas.vx
           tys=gas.vy
           tzs=gas.vz
           gas.vx=(ax*txs+bx*tys+cx*tzs)
           gas.vy=(ay*txs+by*tys+cy*tzs)
           gas.vz=(az*txs+bz*tys+cz*tzs)
; translate dark matter particles  (and change units)
           txs=dark.x
           tys=dark.y
           tzs=dark.z
           dark.x=(ax*txs+bx*tys+cx*tzs)
           dark.y=(ay*txs+by*tys+cy*tzs)
           dark.z=(az*txs+bz*tys+cz*tzs)

           txs=dark.vx
           tys=dark.vy
           tzs=dark.vz
           dark.vx=(ax*txs+bx*tys+cx*tzs)
           dark.vy=(ay*txs+by*tys+cy*tzs)
           dark.vz=(az*txs+bz*tys+cz*tzs)

return
end

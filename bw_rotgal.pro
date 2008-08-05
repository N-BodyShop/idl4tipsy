pro bw_rotgal, INFILE = infile, RBOX = rbox, OUTBOX = outbox, $
               ANGLE = angle, LUNIT = lunit, TYPE = type, HALOID = haloid, $
               MULTIPLE = multiplefile

; NAME:
;   bw_rotgal
;
; PURPOSE:
;     To create a TIPSY file of some or all of a simulation snapshot,
;     rotated such that the angular momentum vector of the primary
;     galaxy is oriented along the z-direction, with an optional
;     angular offset.  
;
; CALLING SEQUENCE:
;  IDL> bw_rotgal, infile = infile
;
; REQUIRED INPUTS:
;     infile = Name of TIPSY snapshot file
;      e.g. 'MW1.1024g1bwK.00512'
;
; OPTIONAL INPUTS:
;
;     rbox = In comoving kpc, the halfwidth of a cube within which particles are
;     included in the angular momentum calculation that will orient
;     the disk.  Generally want
;     this to be small.  Bars can screw you up, so always visually
;     check output with tipsy. Default = 1 comoving kpc.
;
;     outbox = In physical units, the DIAMETER of the box you would
;     like to output in the rotated coordinate system. Default =
;     entire box. If you want the whole box, do not set outbox.
;                                                                
;
;     angle = Sets an angular offset in the rotated coordinate
;     system. Angle of 0 would yield an output box with the galaxy
;     angular momentum vector in the z-direction (edge-on), and an
;     angle of 90 would yield a face on system.  Default - 0.  Only
;                                                          accepts
;                                                          angle > 0.
;
;     lunit = Comoving length unit of the simulation.  Used to convert the rbox
;     and inbox to simulation coordinates. Default = 2.85714e4 kpc 
;
;     type = type of particle used to calculate angular momentum of galaxy.
;     Acceptable inputs are: gas, star, or baryon. Default = gas.
;
;     haloid = The AMIGA grp number assigned to the halo that is being 
;     rotated.  Default = 1 (typically the most massive galaxy in the 
;     simulation. 
;
;     multiplefile = if you wish to orient and cut multiple halos from 
;     one simulation, this file contains the optional inputs for each 
;     halo.  This allows rtipsy to be called just once for the simulation, 
;     drastically reducing runtime.  In this case, the input file MUST 
;     contain a column for EACH of the optional inputs.
;
; OUTPUTS:
;     A rotated (sub)set of the input snapshot in TIPSY format in
;     simulation units.  
;       infile+".size.angle.haloid.rot", where "size" is outbox and "angle"
;       is the number of degrees between disk angular momentum vector
;       and z vector, and haloid is the grp number of the halo.  
;     e.g. MW1.1024g1bwK.00512.025.00.01.rot is a 25 kpc diamter box 
;     around the MW main halo with the angular momentum vector 
;     aligned with the z-direction.
;
;     The size 999 corresponds to a rotation of the whole box.
;     e.g. MW1.1024g1bwK.00512.999.00.rot
;
;     A markfile that contains the indices of all particles in the
;     output box.  The markfile has the format:
;       infile+".size.angle.rot.mark"
;     A markfile is not output for a whole box rotation
;
; EXTERNAL ROUTINES:
;     rtipsy
;     wtipsy
;     az_norm
;     az_dotp
;     readcol
;     
; EXTERNAL DATA:

; COMMENTS:
;   LIMITATIONS = Only accepts angle > 0.  CURRENTLY DOESN"T ACCEPT
;   NON-ZERO ANGLE.  DOES NOT YET OUTPUT A MARKFILE.  The "angle" 
;   input is thus currently ignored. 
;   
;   Currently the code relies on the halo finder being amiga, the halo
;   finder naming convention being infile+.amiga.gtp, and the first
;   element of the .gtp file being the primary halo of the sim box
;
;   It can be helpful to output only a subset of the snapshot in
;   rotation form for purposes of image creation and running the
;   galaxy scope.
;
; REVISION HISTORY:
;   2007-Fall       Charlotte Christensen (optimization to allow it to
;   deal with larger halos)
;   2006-October-28 Alyson Brooks (added ability to choose the haloid 
; 		number, and made "outbox" work) 
;   2006-October-13 Beth Willman, version 2.0
;   2006-August     Adi Zolotov, version 1.0 writes rotation matric

;.read tipsy files
print,infile

rtipsy,STRCOMPRESS(infile+".amiga.gtp"),x,x,x,halo
rtipsy,infile,h,g,d,s

IF n_elements(multiplefile) ne 0 then readcol, multiplefile, rbox, outbox, angle, lunit, type, haloid, format='d,d,d,d,a,l'

IF n_elements(rbox) eq 0 then rbox = 1.0 ;kpc
IF n_elements(outbox) eq 0 then outbox = 999 ;kpc
IF n_elements(type) eq 0 then type = ['gas']
IF n_elements(haloid) eq 0 then haloid = [1]
IF n_elements(angle) eq 0 then angle = [0]

FOR j=0,n_elements(haloid)-1 do begin
    haloind = haloid[j]-1   ; the grp number -1 to index the gtp halos

    IF type[j] eq 'gas' then type[j] = 1
    IF type[j] eq 'star' then type[j] = 2
    IF type[j] eq 'baryon' then type[j] = 3

    IF (outbox[j] LT 10)   THEN outbox_str = '000'+STRN(FIX(outbox[j]))
    IF (outbox[j] LT 100   AND outbox[j] GE 10) THEN outbox_str = '00'+STRN(FIX(outbox[j]))
    IF (outbox[j] LT 1000  AND outbox[j] GE 100) THEN outbox_str = '0'+STRN(FIX(outbox[j]))
    IF (outbox[j] GE 1000) THEN outbox_str = STRN(FIX(outbox[j]))

    IF (angle[j] LT 10) THEN angle_str = '0'+STRN(FIX(angle[j])) ELSE angle_str = STRN(FIX(angle[j]))

    IF (haloid[j] LT 10)   THEN haloid_str = '00'+STRN(FIX(haloid[j]))
    IF (haloid[j] LT 100   AND haloid[j] GE 10) THEN haloid_str = '0'+STRN(FIX(haloid[j]))
    IF (haloid[j] GE 100) THEN haloid_str = STRN(FIX(haloid[j]))

    outfile = STRING(infile + '.' + outbox_str + '.' + angle_str + '.' + haloid_str + '.rot')
    print,j,haloid[j],outfile

;if haloid[j] lt 10.0 and haloid[j] ge 1.0 then begin
;    if angle lt 10 and angle ge 0 then $
;      outfile = STRING(infile+'.000'+STRN(FIX(outbox[j]))+'.0'+STRN(angle)+'.000'+STRN(haloid[j])+'.rot')
;    if angle ge 10 and angle lt 100 then $
;      outfile = STRING(infile+'.000'+STRN(FIX(outbox[j]))+'.'+STRN(angle)+'.000'+STRN(haloid[j])+'.rot')
;endif

;if haloid[j] lt 100.0 and haloid[j] ge 10.0 then begin
;    if angle lt 10 and angle ge 0 then $
;      outfile = STRING(infile+'.00'+STRN(FIX(outbox[j]))+'.0'+STRN(angle)+'.00'+STRN(haloid[j])+'.rot')
;    if angle ge 10 and angle lt 100 then $
;      outfile = STRING(infile+'.00'+STRN(FIX(outbox[j]))+'.'+STRN(angle)+'.00'+STRN(haloid[j])+'.rot')
;endif
;if haloid[j] lt 1000.0 and outbox[j] ge 100.0 then begin
;    if outbox[j] ge 100.0 then begin
;        if angle lt 10 and angle ge 0 then $
;          outfile = STRING(infile+'.0'+STRN(FIX(outbox[j]))+'.0'+STRN(angle)+'.0'+STRN(haloid[j])+'.rot')
;        if angle ge 10 and angle lt 100 then $
;          outfile = STRING(infile+'.0'+STRN(FIX(outbox[j]))+'.'+STRN(angle)+'.0'+STRN(haloid[j])+'.rot')
;    endif
;ENDIF
;if outbox[j] lt 10000.0 and outbox[j] ge 1000.0 then begin
;    if outbox[j] ge 1000.0 then begin
;        if angle lt 10 and angle ge 0 then $
;          outfile = STRING(infile+'.'+STRN(FIX(outbox[j]))+'.0'+STRN(angle)+'.'+STRN(haloid[j])+'.rot')
;        if angle ge 10 and angle lt 100 then $
;          outfile = STRING(infile+'.'+STRN(FIX(outbox[j]))+'.'+STRN(angle)+'.'+STRN(haloid[j])+'.rot')
;    endif
;ENDIF

OPENR, 3, outfile, ERROR = err
CLOSE, 3

;This checks to see if the files has already been generated.  If it
;hasn't continue
IF (err ne 0) then begin

;if wholebox eq 1 then begin
;    if angle lt 10 and angle ge 0 then $
;      outfile = STRING(infile+'.999.0'+STRN(angle,F=I1)+'.rot')
;    if angle ge 10 and angle lt 100 then $
;      outfile = STRING(infile+'.999.'+STRN(angle,F=I2)+'.rot')
;endif
markfile = STRCOMPRESS(outfile+'.mark')

print,'Outbox Number: ',outbox[j]
print, 'Your rotated galaxy will be: ', outfile
print, 'The markfile containing the indices of particles in the output box: ', markfile

; simulation units
IF n_elements(lunit) eq 0 then lunit=2.857d4 else lunit=lunit[0]; kpc per sim unit
; Put rbox into simulation units
rbox1 = rbox[j]/lunit

;Select box of particles to use
;IF type[j] eq 'gas' then begin
IF type[j] eq 1 then begin
    dist_x=g.X-halo[haloind].x
    dist_y=g.Y-halo[haloind].y
    dist_z=g.Z-halo[haloind].z
    pcuts=where((abs(dist_x) lt rbox1) and $
              (abs(dist_y) lt rbox1) and $
                (abs(dist_z) lt rbox1),ncuts)
    IF (N_ELEMENTS(pcuts) LT 5) THEN stop = 1 ELSE stop = 0
    IF (pcuts[0] NE -1) THEN box=g[pcuts]
    print, "Number of gas particles in L box:", ncuts
ENDIF
;IF type[j] eq 'star' then begin
IF type[j] eq 2 then begin
    dist_x=s.X-halo[haloind].x
    dist_y=s.Y-halo[haloind].y
    dist_z=s.Z-halo[haloind].z
    pcuts=where((abs(dist_x) lt rbox1) and $
              (abs(dist_y) lt rbox1) and $
              (abs(dist_z) lt rbox1),ncuts)
    IF pcuts[0] ne -1 then box=s[pcuts]
    print, "Number of star particles in L box:", ncuts
ENDIF
;IF type[j] eq 'baryon' then begin
IF type[j] eq 3 then begin
    dist_x=g.X-halo[haloind].x
    dist_y=g.Y-halo[haloind].y
    dist_z=g.Z-halo[haloind].z
    pcuts1=where((abs(dist_x) lt rbox1) and $
              (abs(dist_y) lt rbox1) and $
              (abs(dist_z) lt rbox1),ncuts)
    print, "Number of gas particles in L box:", ncuts
    dist_x=s.X-halo[haloind].x
    dist_y=s.Y-halo[haloind].y
    dist_z=s.Z-halo[haloind].z
    pcuts2=where((abs(dist_x) lt rbox1) and $
              (abs(dist_y) lt rbox1) and $
              (abs(dist_z) lt rbox1),ncuts)
    print, "Number of star particles in L box:", ncuts
    box=[g(pcuts1),s(pcuts2)]
ENDIF

IF (stop ne 1) then begin 
    print,'Enough gas particles to continue'
    x=box.x-halo[haloind].x
    y=box.y-halo[haloind].y
    z=box.z-halo[haloind].z
    px=(box.vx-halo[haloind].vx)*box.mass
    py=(box.vy-halo[haloind].vy)*box.mass
    pz=(box.vz-halo[haloind].vz)*box.mass

    L_box= dblarr(3)
    for ii=0L,n_elements(box)-1L do L_box= L_box $
      +crossp([x[ii],y[ii],z[ii]],[px[ii],py[ii],pz[ii]])
    zhat=L_box/az_norm(L_box)

;begplot, 'rotation.ps',xsize=7.0,ysize=7.0
;xrange=[-1,1]*rbox
;plot,x,y,psym=3,title='before',xrange=xrange,yrange=xrange
;oplot, [0.0,zhat[0]],[0.0,zhat[1]]
;plot,x,z,psym=3,title='before',xrange=xrange,yrange=xrange
;oplot, [0.0,zhat[0]],[0.0,zhat[2]]
;plot,y,z,psym=3,title='before',xrange=xrange,yrange=xrange
;oplot, [0.0,zhat[1]],[0.0,zhat[2]]

    r_nth=[x[1],y[1],z[1]]
    xtilda=r_nth-zhat*az_dotp(zhat,r_nth)
    xhat=xtilda/az_norm(xtilda)
    yhat=crossp(zhat,xhat)

    IF KEYWORD_SET(verbose) then begin
;print out to make sure everything is ok
        splog,xhat
        splog,yhat
        splog,zhat

        splog, az_dotp(xhat,xhat)
        splog, az_dotp(xhat,yhat)
        splog, az_dotp(xhat,zhat)
        splog, crossp(xhat,yhat)-zhat
        splog, crossp(yhat,zhat)-xhat
        splog, crossp(zhat,xhat)-yhat
   ENDIF

        axis=dblarr(3,3)
        axis[*,0]=xhat
        axis[*,1]=yhat
        axis[*,2]=zhat

;rotate positions of gas

;If cutting a smaller box
        h1 = h
        g1 = g
        s1 = s
        d1 = d
        IF outbox[j] ne 999 then begin
            g1.x = g.x-halo[haloind].x
            g1.y = g.y-halo[haloind].y
            g1.z = g.z-halo[haloind].z
            s1.x = s.x-halo[haloind].x
            s1.y = s.y-halo[haloind].y
            s1.z = s.z-halo[haloind].z
            d1.x = d.x-halo[haloind].x
            d1.y = d.y-halo[haloind].y
            d1.z = d.z-halo[haloind].z
        ENDIF

;rr_g=dblarr(h.ngas,3)
;rr_g[*,0]=g1.x
;rr_g[*,1]=g1.y
;rr_g[*,2]=g1.z
;rotated_g=dblarr(h.ngas,3)
;rotated_g=axis##rr_g
;g1x=rotated_g[*,0]
;g1y=rotated_g[*,1]
;g1z=rotated_g[*,2]

;print,'Gas Arrays Equal? ',ARRAY_EQUAL(g1x,g1.x),' ',ARRAY_EQUAL(g1y,g1.y),' ',ARRAY_EQUAL(g1z,g1.z)

;rotate velocity of gas
;rr_g_v=dblarr(h.ngas,3)
;rr_g_v[*,0]=g1.vx
;rr_g_v[*,1]=g1.vy
;rr_g_v[*,2]=g1.vz
;rotated_g_v=dblarr(h.ngas,3)
;rotated_g_v=axis##rr_g_v
;g1.vx=rotated_g_v[*,0]
;g1.vy=rotated_g_v[*,1]
;g1.vz=rotated_g_v[*,2]

;rotate dark positions
;rr_d=dblarr(h.ndark,3)
;rr_d[*,0]=d1.x
;rr_d[*,1]=d1.y
;rr_d[*,2]=d1.z
;rotated_d=dblarr(h.ndark,3)
;rotated_d=axis##rr_d
;d1x=rotated_d[*,0]
;d1y=rotated_d[*,1]
;d1z=rotated_d[*,2]

;rotate velocity of dark
;rr_d_v=dblarr(h.ndark,3)
;rr_d_v[*,0]=d1.vx
;rr_d_v[*,1]=d1.vy
;rr_d_v[*,2]=d1.vz
;rotated_d_v=dblarr(h.ndark,3)
;rotated_d_v=axis##rr_d_v
;d1.vx=rotated_d_v[*,0]
;d1.vy=rotated_d_v[*,1]
;d1.vz=rotated_d_v[*,2]

; rotate star positions
;rr_s=dblarr(h.nstar,3)
;rr_s[*,0]=s1.x
;rr_s[*,1]=s1.y
;rr_s[*,2]=s1.z
;rotated_s=dblarr(h.nstar,3)
;rotated_s=axis##rr_s
;s1x=rotated_s[*,0]
;s1y=rotated_s[*,1]
;s1z=rotated_s[*,2]

;rotate velocity of stars
;rr_s_v=dblarr(h.nstar,3)
;rr_s_v[*,0]=s1.vx
;rr_s_v[*,1]=s1.vy
;rr_s_v[*,2]=s1.vz
;rotated_s_v=dblarr(h.nstar,3)
;rotated_s_v=axis##rr_s_v
;s1.vx=rotated_s_v[*,0]
;s1.vy=rotated_s_v[*,1]
;s1.vz=rotated_s_v[*,2]

;write tipsy file:
        print, 'Writing output'
        IF outbox[j] eq 999 THEN BEGIN

            g1x = (axis##TEMPORARY([[g1.x],[g1.y],[g1.z]]))[*,0]
            g1y = (axis##TEMPORARY([[g1.x],[g1.y],[g1.z]]))[*,1]
            g1.z = (axis##TEMPORARY([[g1.x],[g1.y],[g1.z]]))[*,2]
            g1.x = g1x
            g1.y = g1y

            g1vx = (axis##TEMPORARY([[g1.vx],[g1.vy],[g1.vz]]))[*,0]
            g1vy = (axis##TEMPORARY([[g1.vx],[g1.vy],[g1.vz]]))[*,1]
            g1.vz = (axis##TEMPORARY([[g1.vx],[g1.vy],[g1.vz]]))[*,2]
            g1.vx = g1vx
            g1.vy = g1vy

            d1x = (axis##TEMPORARY([[d1.x],[d1.y],[d1.z]]))[*,0]
            d1y = (axis##TEMPORARY([[d1.x],[d1.y],[d1.z]]))[*,1]
            d1.z = (axis##TEMPORARY([[d1.x],[d1.y],[d1.z]]))[*,2]
            d1.x = d1x
            d1.y = d1y

            d1vx = (axis##TEMPORARY([[d1.vx],[d1.vy],[d1.vz]]))[*,0]
            d1vy = (axis##TEMPORARY([[d1.vx],[d1.vy],[d1.vz]]))[*,1]
            d1.vz = (axis##TEMPORARY([[d1.vx],[d1.vy],[d1.vz]]))[*,2]
            d1.vx = d1vx
            d1.vy = d1vy

            s1x = (axis##TEMPORARY([[s1.x],[s1.y],[s1.z]]))[*,0]
            s1y = (axis##TEMPORARY([[s1.x],[s1.y],[s1.z]]))[*,1]
            s1.z = (axis##TEMPORARY([[s1.x],[s1.y],[s1.z]]))[*,2]
            s1.x = s1x
            s1.y = s1y

            s1vx = (axis##TEMPORARY([[s1.vx],[s1.vy],[s1.vz]]))[*,0]
            s1vy = (axis##TEMPORARY([[s1.vx],[s1.vy],[s1.vz]]))[*,1]
            s1.vz = (axis##TEMPORARY([[s1.vx],[s1.vy],[s1.vz]]))[*,2]
            s1.vx = s1vx
            s1.vy = s1vy

            wtipsy,outfile,h,g,d,s,/STANDARD 
        
        ENDIF ELSE BEGIN
            indicies = where(((g1.x)^2.+(g1.y)^2.+(g1.z)^2.)^0.5*lunit LE outbox[j])
            if indicies[0] ne -1 then begin 
                g1 = g1(indicies)
                g1x = (axis##TEMPORARY([[g1.x],[g1.y],[g1.z]]))[*,0]
                g1y = (axis##TEMPORARY([[g1.x],[g1.y],[g1.z]]))[*,1]
                g1.z = (axis##TEMPORARY([[g1.x],[g1.y],[g1.z]]))[*,2]
                g1.x = g1x
                g1.y = g1y

                g1vx = (axis##TEMPORARY([[g1.vx],[g1.vy],[g1.vz]]))[*,0]
                g1vy = (axis##TEMPORARY([[g1.vx],[g1.vy],[g1.vz]]))[*,1]
                g1.vz = (axis##TEMPORARY([[g1.vx],[g1.vy],[g1.vz]]))[*,2]
                g1.vx = g1vx
                g1.vy = g1vy
                h1.ngas = n_elements(g1.x)
            endif else begin 
                h1.ngas = 0
                g1 = g1[{mass:0, x:0, y:0, z:0, vx:0, vy:0, vz:0, dens:0, tempg:0, h:0, zmetal:0, phi:0}]
            endelse
            indicies = where(((s1.x)^2.+(s1.y)^2.+(s1.z)^2.)^0.5*lunit LE outbox[j]) 
            if indicies[0] ne -1 then begin             
                s1 = s1(indicies)
                s1x = (axis##TEMPORARY([[s1.x],[s1.y],[s1.z]]))[*,0]
                s1y = (axis##TEMPORARY([[s1.x],[s1.y],[s1.z]]))[*,1]
                s1.z = (axis##TEMPORARY([[s1.x],[s1.y],[s1.z]]))[*,2]
                s1.x = s1x
                s1.y = s1y

                s1vx = (axis##TEMPORARY([[s1.vx],[s1.vy],[s1.vz]]))[*,0]
                s1vy = (axis##TEMPORARY([[s1.vx],[s1.vy],[s1.vz]]))[*,1]
                s1.vz = (axis##TEMPORARY([[s1.vx],[s1.vy],[s1.vz]]))[*,2]
                s1.vx = s1vx
                s1.vy = s1vy
                h1.nstar = n_elements(s1.x)
            endif else begin 
                h1.nstar = 0
                s1 = [{mass:0, x:0, y:0, z:0, vx:0, vy:0, vz:0, metals:0, tform:0, eps:0, phi:0}]
            endelse
            indicies = where(((d1.x)^2.+(d1.y)^2.+(d1.z)^2.)^0.5*lunit LE outbox[j])
            if indicies[0] ne -1 then begin 
                d1 = d1(indicies)
                d1x = (axis##TEMPORARY([[d1.x],[d1.y],[d1.z]]))[*,0]
                d1y = (axis##TEMPORARY([[d1.x],[d1.y],[d1.z]]))[*,1]
                d1.z = (axis##TEMPORARY([[d1.x],[d1.y],[d1.z]]))[*,2]
                d1.x = d1x
                d1.y = d1y

                d1vx = (axis##TEMPORARY([[d1.vx],[d1.vy],[d1.vz]]))[*,0]
                d1vy = (axis##TEMPORARY([[d1.vx],[d1.vy],[d1.vz]]))[*,1]
                d1.vz = (axis##TEMPORARY([[d1.vx],[d1.vy],[d1.vz]]))[*,2]
                d1.vx = d1vx
                d1.vy = d1vy
                h1.ndark = n_elements(d1.x)
            endif else begin
                h1.ndark = 0
                d1 = [{mass:0, x:0, y:0, z:0, vx:0, vy:0, vz:0, eps:0, phi:0}]
            endelse
            h1.n = h1.ngas+h1.ndark+h1.nstar

            wtipsy, outfile, h1,g1,d1,s1, /standard
        ENDELSE
    endif
ENDIF ELSE print,'File ',outfile,' already exists'
ENDFOR

;oldrange= xrange
;xrange= (oldrange+xcenter)/MW1_length
;yrange= (oldrange+ycenter)/MW1_length
;zrange= (oldrange+zcenter)/MW1_length
;;plot,g[gascuts].X,g[gascuts].Y,psym=3,xrange=xrange,yrange=yrange
;;plot,g[gascuts].X,g[gascuts].Z,psym=3,xrange=xrange,yrange=zrange
;;plot,g[gascuts].Y,g[gascuts].Z,psym=3,xrange=yrange,yrange=zrange

;rotcentvel=axis##[vxcenter,vycenter,vzcenter]
;center= axis##[xcenter,ycenter,zcenter]

;str_cen=create_struct('centerx',0.D,'centery',0.D,$
;                      'centerz',0.D,'cenvelx',0.D,$
;                      'cenvely',0.D,'cenvelz',0.D)
;censtr=replicate(str_cen,6)
;censtr.centerx=center[0]
;censtr.centery=center[1]
;censtr.centerz=center[2]
;censtr.cenvelx=rotcentvel[0]
;censtr.cenvely=rotcentvel[1]
;censtr.cenvelz=rotcentvel[2]


;dx= (g[gascuts].X*MW1_length-center[0]) ; kpc
;dy= (g[gascuts].Y*MW1_length-center[1]) ; kpc
;dz= (g[gascuts].Z*MW1_length-center[2]) ; kpc
;;plot,dx,dy,psym=3,xrange=oldrange,yrange=oldrange
;;plot,dx,dz,psym=3,xrange=oldrange,yrange=oldrange
;;plot,dy,dz,psym=3,xrange=oldrange,yrange=oldrange
;;endplot

RETURN
end

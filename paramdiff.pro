;Takes two gasoline parameter files (pfile1 and pfile2) and compares their contents.  
;Use the v(erbose) option to show all the parameters, not just the different ones.
;Charlotte Christensen
;5/22/12
;
;Note that this program assumes that the '=' file in the parameter
;file is seporated by spaces on both sides
;
PRO PARAMDIFF,pfile1,pfile2,v = v

cm_per_kpc = 3.08568021d21
gm_per_msol = 1.98892d33
amu_per_gm = 6.022142d23
gm_per_H = 1.673534d-24
H_per_gm = 5.9753790e+23

READCOL, pfile1,param1,equals1,value1,FORMAT='A,A,A',/SILENT
READCOL, pfile2,param2,equals2,value2,FORMAT='A,A,A',/SILENT

;only look through items that are not commented and are properly formatted
a=STRMATCH(equals1,'=')
b=STRMATCH(param1,'*#*')
ind=WHERE(a EQ 1 AND b EQ 0)
param1=param1[ind]
value1=value1[ind]

a=STRMATCH(equals2,'=')
b=STRMATCH(param2,'*#*')
ind=WHERE(a EQ 1 AND b EQ 0)
param2=param2[ind]
value2=value2[ind]

ind1 = fltarr(N_ELEMENTS(param1)) - 1
ind2 = fltarr(N_ELEMENTS(param2)) - 1
print,''
print,' -------------------------- ',pfile1,pfile2,FORMAT='(T1,A,T30,A,T58,A)'
FOR i1 = 0, N_ELEMENTS(param1) - 1 DO BEGIN

    i2 = where(param2 eq param1[i1])
    ind1[i1] = i2
    IF (i2[0] NE -1) THEN BEGIN ;If the parameter appears in the second file
        IF (N_ELEMENTS(i2) gt 1) THEN BEGIN
            print,'Duplicate: ',param1[i1]
            i2 = i2[0]
        ENDIF
        ind2[i2] = i1
        IF STRMATCH(value1[i1],value2[i2]) THEN BEGIN 
            IF KEYWORD_SET(v) THEN print,' + ',param1[i1],FORMAT= '(T1,A,T4,A,T30,A)',value1[i1]
        ENDIF ELSE BEGIN
            print,' - ',param1[i1],FORMAT= '(T1,A,T4,A,T30,A, T58,A)',value1[i1],value2[i2]
        ENDELSE
    ENDIF
ENDFOR
print,' --------------------------- ',pfile1,'---------',FORMAT='(T1,A,T30,A,T58,A)'
FOR i1 = 0, N_ELEMENTS(param1) - 1 DO BEGIN
    IF (ind1[i1] EQ -1) THEN BEGIN
        print,' - ',param1[i1],FORMAT= '(T1,A,T4,A,T30,A)',value1[i1]
    ENDIF
ENDFOR

print,' ------------------------------------------------------- ',pfile2
FOR i2 = 0, N_ELEMENTS(param2) - 1 DO BEGIN
    IF (ind2[i2] EQ -1) THEN BEGIN
        print,' - ',param2[i2],FORMAT= '(T1,A,T4,A,T30,A)',value2[i2]
    ENDIF
ENDFOR

END

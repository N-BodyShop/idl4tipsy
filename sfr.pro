PRO SFR, s, meansfr,binsize=binsize, OVERPLOT=overplot,move=move, $
       NOPLOT=noplot,YMAX=ymax,LINESTYLE=linestyle,MINT=mint,maxt=maxtime, $
       MASSUNIT=massunit,TIMEUNIT=timeunit,MASSFORM=massform,sfh=sfr, $
       tarray=timearray,dob=dob,AGE=age,XLOG=xlog,ylog=ylog,nbins=nbins, $
       _EXTRA=_extra
;plots the star formation rate for a simulation

IF (keyword_set (timeunit) EQ 0) THEN timeunit=4.0435875e10
IF (keyword_set(massunit) EQ 0) THEN massunit=13.6e16

maxti=MAX(s.tform)
ind=WHERE(s.tform GT 0.0)
if ( keyword_set(move) EQ 1 ) then  tform=s[ind].tform*timeunit + move $
else if (keyword_set(age) EQ 1) then tform=(maxti-s[ind].tform)*timeunit $
else tform=s[ind].tform*timeunit

IF (keyword_set(massform) EQ 0) THEN mass=s[ind].mass*massunit $
ELSE mass=massform[ind]*massunit

IF (keyword_set(mint) EQ 0) THEN mintime=MIN(tform) ELSE mintime=mint
IF (keyword_set(maxt) EQ 0) THEN maxtime=MAX(tform)

IF (keyword_set(xlog) EQ 1) then begin
  IF (keyword_set(nbins) EQ 0) then nbins=50.
  if(mintime LT maxtime/1e2) then mintime=maxtime/1e2
  logmintime=alog10(mintime)
  logbinsize=(alog10(maxtime)-logmintime)/nbins
  binsize = 10.^(findgen(nbins)*logbinsize +logmintime)
  timearray = mintime+binsize
  midtimearray = mintime+10.^(findgen(nbins)*logbinsize+0.5*logbinsize +logmintime)
endif else begin
IF (keyword_set(binsize) EQ 0) THEN binsize=5.e7
nbins=FIX((maxtime-mintime)/binsize)+1
timearray=FINDGEN(nbins)*binsize+mintime
midtimearray=FINDGEN(nbins)*binsize+mintime+binsize*0.5
binsize=fltarr(nbins)+binsize
endelse

sfr=FLTARR(nbins)
b=sfr

FOR i=0,nbins-1 DO BEGIN
    inbin=WHERE(tform GE timearray[i] AND tform LT (timearray[i]+binsize[i]),ninbin)
   IF (ninbin EQ 0) THEN sfr[i]=0 $
   ELSE sfr[i]=TOTAL(mass[inbin])/binsize[i]

   IF (KEYWORD_SET(dob) AND i GT 1) $
     THEN b[i]=sfr[i]/mean(sfr[0:i-1])$
   ELSE b[i]=0
ENDFOR

;STOP

IF (KEYWORD_SET(dob)) THEN BEGIN
  ytitle="b" 
  sfr = b
ENDIF else ytitle='SFR [M'+sunsymbol()+' yr!u-1!n]'

IF (KEYWORD_SET(ymax) EQ 0) THEN ymax=MAX(sfr)
 
IF (KEYWORD_SET(noplot) EQ 0) THEN BEGIN
;IF (keyword_set(xlog) eq 0) then xlog=0    
IF (keyword_set(age) eq 1) then xtit="Stellar Age [Gyr]" ELSE xtit="Time [Gyr]"
ymin = 0
if (keyword_set(ylog)) then begin
  ymax = 2*MAX(sfr) 
  ymin=min(sfr[where(sfr GT 0)])
endif else ylog=0

; ---------PLOTTING--------------------
IF KEYWORD_SET(overplot) THEN $
oplot, midtimearray/1e9,sfr,psym=10,_EXTRA=_extra $
ELSE plot,midtimearray/1e9,sfr,psym=10, yrange=[ymin,ymax], $
       xtitle=xtit,ytitle=ytitle,ylog=ylog,_EXTRA=_extra 
ENDIF ; noplot

IF N_PARAMS() EQ 2 THEN BEGIN
    mint=mintime+5.e8
    maxt=mintime+15.e8
    ind=WHERE(midtimearray GE mint AND midtimearray LE maxt)
    meansfr=MEAN(sfr[ind])
    IF (KEYWORD_SET(noplot) EQ 0) THEN BEGIN
       ;if ( ymax GT max(sfr)) then dashmax=ymax
       plots,[mint/1e9,mint/1e9],[0,ymax],linestyle=2
       plots,[maxt/1e9,maxt/1e9],[0,ymax],linestyle=2
    ENDIF
ENDIF

END

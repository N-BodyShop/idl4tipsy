;restore,'skidstattemp.sav'
;s = read_ascii('fmt.01024.stat',template=temp)
prefix = 'fmtss'
step = '00512'
file = prefix+'.'+step
idlstatfile = file+'.idl.stat'
minlambda = 0.03
maxlambda = 0.05
minmass = 5e11
maxmass = 1e12

if (file_test(idlstatfile) eq 0) then begin
  s = calcidlstats(file)
endif else begin
  s = ridlstat(idlstatfile)
endelse

print,'picking halos'
gids = s[where(s.lambda GT minlambda AND s.lambda LT maxlambda AND  $
	       s.mass GT minmass AND s.mass LT maxmass AND $
	       s.maxneighbormass LT minmass)].gid

neg = n_elements(gids)

openw,lun,'halfmasszs.dat',/get_lun
for i =0,neg-1 do begin
  if (gids[i] eq 0) then continue
  print,'Making merger tree for gid ',gids[i]
  txtfile = 'g'+strtrim(gids[i],2)+file+'.txt'
 ; if(file_test(txtfile) LT 1) then 
  mergertree,prefix,gids[i],zs=zs,grp=grps,tipsys=tipsys,halfmassz=halfmassz
  printf,lun,gids[i],halfmassz
endfor
close,lun
free_lun,lun

END

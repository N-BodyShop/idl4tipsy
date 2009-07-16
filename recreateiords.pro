
timeu=38.7676659d9  ; Double precision required
;timeu = dKpcUnit / sqrt(G dMsolUnit / dKpcUnit)
ngasorig=961397
; for f in `ls sim.0????`; do ${TIPSY_TOOLS}/printstdhead < $f; done>heads
readcol,'heads',n,ngas,nstar,ndim,a,format='l,l,l,l,f'
ndark = n-ngas-nstar
; ls sim.0???? > out
readcol,'outs',outfile,format='a'
nouts = n_elements(n)
print,'reading starlog'
sl=rstarlog('g1536')
lastiord = rbarray('g1536.01024.iord',type='long')
lastgasiord = lastiord[0:ngas[nouts-1]-1]

gasstart = lonarr(ngasorig)+1
nprev=ngasorig
for i=0L,nouts-1 do begin
    print,outfile[i]
    ; come up with all possible surviving gas particles
    ; duplicates sorted out via uniq later

    t = (lookbacka(1e-7)-lookbacka(a[i]))/timeu
    latergas = sl[where(sl.timeform gt t+3e-7)].iordergas
    latergas = latergas[uniq(latergas,sort(latergas))]

    giords = [lastgasiord, latergas];, ngs]

    iords = lonarr(n[i])
    giords = giords[uniq(giords,sort(giords))]
    iords[0:ngas[i]-1]=giords
    iords[ngas[i]:ngas[i]+ndark[i]-1] = lindgen(ndark[i])+ngasorig
    iords[ngas[i] + ndark[i]:n[i]-1]=lindgen(nstar[i])+min(sl.iorderstar)
    openw,lun,outfile[i]+'.iord',/get_lun,/xdr
    writeu,lun,n[i],iords
    free_lun,lun
    nprev=ngas[i]
endfor

end

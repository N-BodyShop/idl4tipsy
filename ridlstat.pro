function ridlstat,idlstatfile
  print,'reading idlstatfile'
  record={gid:0L,mass:0.0,lambda:0.0,maxr:0.0,maxneighbormassr:0.0,$
	  comx:0.0,comy:0.0,comz:0.0,comvx:0.0,comvy:0.0,comvz:0.0, $
	  jx:0.0,jy:0.0,jz:0.0,jtot:0.0,nneighbors:0,maxneighbormass:0.0}
  openr,lun,idlstatfile,/get_lun,/xdr
  fs = fstat(lun)
  s = replicate(record, fs.size/n_bytes(record))
  readu,lun,s
  close,lun
  free_lun,lun
  return,s
END

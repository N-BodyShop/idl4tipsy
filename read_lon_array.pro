; like read_ascii_array only reads a long integer rather than float
FUNCTION read_lon_array,filename

  openr,lun,filename,/get_lun
  readf,lun,arrsize
  readcol,filename,arrsize,format='L',numline=1
  array = intarr(arrsize)
  readf,lun,array
  close,lun
  free_lun,lun
  return, array

END

; bin2ascii
; Charlotte Christensen, 1/22/10

;This program reads a tipsy array in binary format and outputs the
;data in ascii format

;If your filename is XXXXXXX.HI, your ascii file will be XXXXXXX.HI.ascii

;This code relies on the binary template, binArray2Ascii_template.sav

;Example:
;  >bin2ascii,'gas_merger0.1_single.00060.H2'

pro bin2ascii,filename
RESTORE, 'binArray2Ascii_template.sav'
array = read_binary(filename,template = temp)
openw,1,filename+'.ascii'
printf,1,array.num
printf,1,TRANSPOSE(array.data)
close,1
end

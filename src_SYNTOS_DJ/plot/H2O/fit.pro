PRO fit
set_plot, 'ps'
device,/encapsulated,/color,xsize=30,ysize=25,filename='fit.eps'
loadct,39

N1 = 553
N2 = 3532


data = fltarr(3,N1)
openr,lun,'2MASSJ00583814-1747311.txt',/get_lun
skip_lun,lun,15,/lines
data = fltarr(3,N1)
readf,lun,data
close,/all
data = data[*,where(data[0,*] GE 0.7,/NULL)]
data = data[*,where(data[0,*] LE 2.5,/NULL)]
a = size(data)
N1 = a[2]


model = fltarr(5,N2)
openr,lun,'t2900g45_Schwenke_spec.dat',/get_lun
skip_lun,lun,1,/lines
readf,lun,model
close,/all
model2 = fltarr(5,N2)
for i=0,N2-1 do begin
	model2[*,i] = model[*,N2-1-i]
endfor
model = model2
result = interpol(model[2,*],model[0,*],data[0,*])
model = result


x = reform(data[0,*],N1)
y = reform(data[1,*],N1)
err = reform(data[2,*],N1)
z = reform(model,N1)

functargs = {x:x, y:y, err:err, z:z}
start_param = [3e8]

k_fit = mpfit('fitfunc',start_param,functargs=functargs,perror=perror)

c1 = k_fit[0]
print, c1, ' +/- ', perror[0]
plot, x, z/c1, xrange=[0.7,2.5]
oplot, x, y,color=240
device,/close
end

function fitfunc, param, x=x, y=y, err=err, z=z
	model = z/param[0]
	return, (y-model)
	
end
  
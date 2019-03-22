set_plot, 'ps'
device,/encapsulated,/color,xsize=25,ysize=35,filename='h2o.eps'
loadct,39

N1 = 3532
N2 = 553
SCAN = fltarr(5,N1)
HITEMP = fltarr(5,N1)
Schwenke = fltarr(5,N1)
data = fltarr(3,N2)


openr,lun,'t2900g45_SCAN_spec.dat',/get_lun
skip_lun,lun,1,/lines
readf,lun,SCAN
close,/all
scale1 = 3.16086e8 

openr,lun,'t2900g45_HITEMP_spec.dat',/get_lun
skip_lun,lun,1,/lines
readf,lun,HITEMP
close,/all
scale2 = 3.45e8

openr,lun,'t2900g45_Schwenke_spec.dat',/get_lun
skip_lun,lun,1,/lines
readf,lun,Schwenke
close,/all
scale3 = 3.05230e8

openr,lun,'2MASSJ00583814-1747311.txt',/get_lun
skip_lun,lun,15,/lines
readf,lun,data
close,/all


!p.multi=[0, 1, 4]
!x.thick = 5
!y.thick = 5

plot, SCAN[0,*], SCAN[1,*]/SCAN[2,*],/nodata,color=0,xrange=[0.7,2.5],xstyle=1,yrange=[0,1],ystyle=1,charsize=2.5,XTitle='Wavelength [!4' + STRING("154B) + '!Xm]',YTitle='Normalised flux'
oplot, SCAN[0,*], SCAN[2,*]/SCAN[1,*],color=60,thick=3
oplot, Schwenke[0,*], Schwenke[2,*]/Schwenke[1,*],color=160,thick=3
oplot,HITEMP[0,*],HITEMP[2,*]/HITEMP[1,*],color=240,thick=3
xyouts,2.4,0.40,'SCAN',color=60,charsize=1.5,alignment=1
xyouts,2.4,0.25,'HITEMP',color=240,charsize=1.5,alignment=1
xyouts,2.4,0.1,'Partridge & Schwenke',color=160,charsize=1.5,alignment=1


plot,  SCAN[0,*], SCAN[2,*]/scale1, color=0, /nodata,xrange=[0.7,2.5],xstyle=1,ystyle=1,yrange=[0,1.1],charsize=2.5,XTitle='Wavelength [!4' + STRING("154B) + '!Xm]',YTitle='Flux'
oplot,  SCAN[0,*], SCAN[2,*]/scale1, color=0,thick=3
oplot, data[0,*], data(1,*),color=240,thick=3
xyouts,2.4,0.85,'SCAN',charsize=1.5,alignment=1
xyouts,0.9,0.1,'2MASSJ00583814-1747311',color=240,charsize=1.0

plot,  HITEMP[0,*], HITEMP[2,*]/scale2, color=0, /nodata,xrange=[0.7,2.5],xstyle=1,ystyle=1,yrange=[0,1.1],charsize=2.5,XTitle='Wavelength [!4' + STRING("154B) + '!Xm]',YTitle='Flux'
oplot,  HITEMP[0,*], HITEMP[2,*]/scale2, color=0,thick=3
oplot,  HITEMP[0,*], HITEMP[2,*]/3e8, color=0,thick=3,linestyle=1
oplot, data[0,*], data(1,*),color=240,thick=3
xyouts,2.4,0.85,'HITEMP',charsize=1.5,alignment=1
xyouts,0.9,0.1,'2MASSJ00583814-1747311',color=240,charsize=1.0

plot,  Schwenke[0,*], Schwenke[2,*]/scale3, color=0, /nodata,xrange=[0.7,2.5],xstyle=1,ystyle=1,yrange=[0,1.1],charsize=2.5,XTitle='Wavelength [!4' + STRING("154B) + '!Xm]',YTitle='Flux'
oplot,  Schwenke[0,*], Schwenke[2,*]/scale3, color=0,thick=3
oplot, data[0,*], data(1,*),color=240,thick=3
xyouts,2.4,0.85,'Partridge & Schwenke',charsize=1.5,alignment=1
xyouts,0.9,0.1,'2MASSJ00583814-1747311',color=240,charsize=1.0



device,/close

end
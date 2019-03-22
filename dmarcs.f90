program dmarcs

integer :: idriftok, io
character :: line*27
real*8  :: Tcormx

do iii=1,2
  print *, 'Iteration ', iii
  print *
  print *, 'Calling MARCS.'
  call system('./marcs_dj_2')
  call system('rm fort*')

!  print *, 'Checking convergence.'
!  open(unit=1,file='Tcorr.txt',status='old')
!  read(1,'(f10.1)') Tcormx
!  if(Tcormx .le. 10.) exit 
!  close(1)

  print *, 'Calling DRIFT.'
  call system('cp drift2marcs.dat d2m.save')
  call system('./static_weather8 > drift.out')
  
  print *, 'Checking output.'
  open(unit=2,file='drift.out',readonly)
  idriftok = 0
  do
    read(2,'(a27)',iostat=io) line(1:27)
    if(io .lt. 0) exit
    if(index(line,' regular end of integration')) then
      idriftok = 1
      exit
    end if
  end do
  close(2)
  
  if(idriftok .ne. 1) then
    print *, 'DRIFT did not converge, using old DRIFT file'
    call system('cp d2m.save drift2marcs.dat')
  end if
  
  call system('rm restart.dat')
  call system('cp arcivaab.dat arcivaaa.dat')
end do

print *, 'DRIFT-MARCS converved successfully! :)'



end

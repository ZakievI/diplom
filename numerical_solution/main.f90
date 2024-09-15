program main
  use mod
  call testing()
  h=5.0d0
  R_1=d1
  R_2=(h+R_1)/2
  nj=100
  ds=pi/nj
  ds2=ds*h/R_1
  F_m=1
  call main_1(2) !1-постоянен, 2-1/r^^3, 3-Ar+b1
  !call main_2
  !call main_3(2) !1-квадртат с вырезом, разрез сил горизонтален, 2-вертикален
  
  
  
  !OPEN (10,FILE='E.dat')
  !do nj=100,1000,100
  !  ds=pi/nj
  !  ds2=ds*h/R_1
  !end do
  !close(10)
  
  pause
  !call main_3(2)
  !call main_3_cheking
end
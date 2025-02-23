subroutine pg_init_boundline_geom(ibndl,npanel,x,y)
!dec$ attributes dllexport:: pg_init_boundline_geom
!инициализировать геометрию для участка границы
use gen_mod
integer(4) ibndl !номер участка границы
integer(4) npanel !количество панелей в текущем участке
real(8) x(npanel+1) !x - концов панелей
real(8) y(npanel+1) !y - концов панелей
call init_boundline_geom(ibndl,npanel,x,y,.false.)
call init_geom_detale1(ibndl,1,npanel+1)
end

subroutine init_boundline_geom(ibndl,npanel,x,y,invert)
!инициализировать геометрию для участка границы
use pgmod
integer(4) ibndl !номер участка границы
integer(4) npanel !количество панелей в текущем участке
real(8) x(npanel+1) !x - концов панелей
real(8) y(npanel+1) !y - концов панелей
type(TBoundline), pointer :: b
logical invert
integer(4) i,j
call deallocate_boundline(gsarea%i,gsbnd%i,ibndl)
call allocate_boundline_geom(ibndl,npanel)
b=>gsbnd%line(ibndl)
if (invert) then
  do i=1,npanel+1
    j=npanel-i+2
    b%x(i)=x(j)
    b%y(i)=y(j)
  enddo
else
  b%x=x
  b%y=y
endif
end

subroutine init_geom_detale_bndl(ibndl)
use pgmod
type(TBoundline), pointer :: b
type(TBoundline_geomdetale), pointer :: gd
integer(4) ibndl
b=>gsbnd%line(ibndl)
gd=>gsbnd%geom_detale0(gsbnd%ngeom_detale0)
if (b%igd_begin==0) b%igd_begin=gsbnd%ngeom_detale0
b%igd_end=gsbnd%ngeom_detale0
end

subroutine init_geom_detale1(ibndl,i_begin,i_end)
use pgmod
type(TBoundline_geomdetale), pointer :: gd
integer(4) ibndl,i_begin,i_end
call allocate_geom_detale_for_next
gd=>gsbnd%geom_detale0(gsbnd%ngeom_detale0)
gd%mode=1
gd%ibndl=ibndl
gd%i_begin=i_begin
gd%i_end=i_end
call init_geom_detale_bndl(ibndl)
end

subroutine pg_add_boundline_geom(ibndl,npanel,x,y)
!dec$ attributes dllexport:: pg_add_boundline_geom
!инициализировать геометрию для участка границы
!добавить к существующей
use gen_mod
integer(4) ibndl !номер участка границы
integer(4) npanel !количество панелей в текущем участке
real(8) x(npanel+1) !x - концов панелей
real(8) y(npanel+1) !y - концов панелей
integer(4) n
call add_boundline_geom(ibndl,npanel,x,y,n)
call init_geom_detale1(ibndl,n,npanel+n)
end

subroutine add_boundline_geom(ibndl,npanel,x,y,n)
!инициализировать геометрию для участка границы
!добавить к существующей
use pgmod
integer(4) ibndl !номер участка границы
integer(4) npanel !количество панелей в текущем участке
real(8) x(npanel+1) !x - концов панелей
real(8) y(npanel+1) !y - концов панелей
real(8), allocatable :: x1(:),y1(:)
integer(4) npanel1,n
type(TBoundline), pointer :: b
b=>gsbnd%line(ibndl)
n=b%npanel+1
npanel1=b%npanel+npanel
allocate(x1(npanel1+1),y1(npanel1+1))
x1(1:b%npanel+1)=b%x
y1(1:b%npanel+1)=b%y
x1(b%npanel+1:npanel1+1)=x
y1(b%npanel+1:npanel1+1)=y
call init_boundline_geom(ibndl,npanel1,x1,y1,.false.)
deallocate(x1,y1)
end

subroutine pg_init_boundline_geomline(ibndl,npanel,x1,y1,x2,y2)
!dec$ attributes dllexport:: pg_init_boundline_geomline
!инициализировать геометрию для участка границы в виде прямой линии
use gen_mod
integer(4) ibndl !номер участка границы
integer(4) npanel !количество панелей в текущем участке
real(8) x1,y1,x2,y2  !координаты концов отрезка
real(8) x(npanel+1) !x - концов панелей
real(8) y(npanel+1) !y - концов панелей
call init_boundline_geomline_array(npanel,x1,y1,x2,y2,x,y)
call init_boundline_geom(ibndl,npanel,x,y,.false.)
call init_geom_detale2(ibndl,1,npanel+1,x1,y1,x2,y2)
end

subroutine init_geom_detale2(ibndl,i_begin,i_end,x1,y1,x2,y2)
use pgmod
type(TBoundline_geomdetale), pointer :: gd
integer(4) ibndl,i_begin,i_end
real(8) x1,y1,x2,y2
call allocate_geom_detale_for_next
gd=>gsbnd%geom_detale0(gsbnd%ngeom_detale0)
gd%mode=2
gd%ibndl=ibndl
gd%i_begin=i_begin
gd%i_end=i_end
gd%x=x1
gd%y=y1
gd%x2=x2
gd%y2=y2
call init_geom_detale_bndl(ibndl)
end

subroutine pg_add_boundline_geomline(ibndl,npanel,x1,y1,x2,y2)
!dec$ attributes dllexport:: pg_add_boundline_geomline
!инициализировать геометрию для участка границы в виде прямой линии
!добавить к существующей
use gen_mod
integer(4) ibndl !номер участка границы
integer(4) npanel !количество панелей в текущем участке
real(8) x1,y1,x2,y2  !координаты концов отрезка
real(8) x(npanel+1) !x - концов панелей
real(8) y(npanel+1) !y - концов панелей
integer(4) n
call init_boundline_geomline_array(npanel,x1,y1,x2,y2,x,y)
call add_boundline_geom(ibndl,npanel,x,y,n)
call init_geom_detale2(ibndl,n,npanel+n,x1,y1,x2,y2)
end

subroutine init_boundline_geomline_array(npanel,x1,y1,x2,y2,x,y)
!инициализировать геометрию для участка границы в виде прямой линии
use gen_mod
integer(4) npanel !количество панелей в текущем участке
real(8) x1,y1,x2,y2  !координаты концов отрезка
real(8) x(npanel+1) !x - концов панелей
real(8) y(npanel+1) !y - концов панелей
integer(4) i
real(8) t,sred
do i=1,npanel+1
  t=(i-d1)/npanel
  x(i)=sred(d0,d1,x1,x2,t)
  y(i)=sred(d0,d1,y1,y2,t)
enddo
end

subroutine init_npanel(ds,x1,y1,x2,y2,npanel)
use pgmod
integer(4) npanel
real(8) ds          !примерный шаг
real(8) x1,y1,x2,y2  !координаты концов отрезка
integer(4) min_np
real(8) s,DistanceToPoint,eps
min_np=gs%const%di*2+4
eps=1.0d-6
s=DistanceToPoint(x1,y1,x2,y2)+eps
npanel=s/ds
if (npanel<min_np) npanel=min_np
end

subroutine pg_init_boundline_geomlineds(ibndl,ds,x1,y1,x2,y2)
!dec$ attributes dllexport:: pg_init_boundline_geomlineds
!инициализировать геометрию для участка границы в виде прямой линии с фиксированным ds
use pgmod
integer(4) ibndl !номер участка границы
integer(4) npanel !количество панелей в текущем участке
real(8) ds          !примерный шаг
real(8) x1,y1,x2,y2  !координаты концов отрезка
call init_npanel(ds,x1,y1,x2,y2,npanel)
call pg_init_boundline_geomline(ibndl,npanel,x1,y1,x2,y2)
end

subroutine pg_add_boundline_geomlineds(ibndl,ds,x1,y1,x2,y2)
!dec$ attributes dllexport:: pg_add_boundline_geomlineds
!инициализировать геометрию для участка границы в виде прямой линии с фиксированным ds
!добавить к существующей
use pgmod
integer(4) ibndl !номер участка границы
integer(4) npanel !количество панелей в текущем участке
real(8) ds          !примерный шаг
real(8) x1,y1,x2,y2  !координаты концов отрезка
call init_npanel(ds,x1,y1,x2,y2,npanel)
call pg_add_boundline_geomline(ibndl,npanel,x1,y1,x2,y2)
end

subroutine pg_init_boundline_geompolylineds(ibndl,ds,x,y,n)
!dec$ attributes dllexport:: pg_init_boundline_geompolylineds
!инициализировать геометрию для участка границы в виде полигона с фиксированным ds
use pgmod
integer(4) ibndl !номер участка границы
integer(4) npanel !количество панелей в текущем участке
real(8) ds          !примерный шаг
integer(4) n       !колчество отрезков в полилинии
real(8) x(n+1),y(n+1)  !координаты концов отрезков полилинии
integer(4) min_np
integer(4) np(n),i,k,k1
real(8), allocatable :: xx(:),yy(:)
real(8) s,DistanceToPoint,eps
min_np=gs%const%di*2+4
eps=1.0d-6
npanel=0
do i=1,n
  s=DistanceToPoint(x(i),y(i),x(i+1),y(i+1))+eps
  np(i)=s/ds
  if (np(i)<min_np) np(i)=min_np
  npanel=npanel+np(i)
enddo
allocate(xx(npanel+1),yy(npanel+1))
k=1
do i=1,n
  k1=k+np(i)
  call init_boundline_geomline_array(np(i),x(i),y(i),x(i+1),y(i+1),xx(k:k1),yy(k:k1))
  call init_geom_detale2(ibndl,k,k1,x(i),y(i),x(i+1),y(i+1))
  k=k1
enddo
call init_boundline_geom(ibndl,npanel,xx,yy,.false.)
deallocate(xx,yy)
end

subroutine pg_init_boundline_geomlineds2(ibndl,ds1,ds2,x1,y1,x2,y2)
!dec$ attributes dllexport:: pg_init_boundline_geomlineds2
!инициализировать геометрию для участка границы в виде прямой линии с переменным ds
!меняющимся от ds1 до ds2
use pgmod
integer(4) ibndl !номер участка границы
integer(4) npanel !количество панелей в текущем участке
real(8) x1,y1,x2,y2  !координаты концов отрезка
real(8) ds1          !шаг в начале отрезка
real(8) ds2          !шаг в конце отрезка
type(TAreaValue) x,y
integer(4) min_np
min_np=gs%const%di*2+4
call init_boundline_geomlineds2_array(ds1,ds2,x1,y1,x2,y2,x,y,npanel,min_np)
call init_boundline_geom(ibndl,npanel,x%v,y%v,.false.)
deallocate(x%v,y%v)
call init_geom_detale2(ibndl,1,npanel+1,x1,y1,x2,y2)
end

subroutine pg_add_boundline_geomlineds2(ibndl,ds1,ds2,x1,y1,x2,y2)
!dec$ attributes dllexport:: pg_add_boundline_geomlineds2
!инициализировать геометрию для участка границы в виде прямой линии с переменным ds
!меняющимся от ds1 до ds2
!добавить к существующей
use pgmod
integer(4) ibndl !номер участка границы
integer(4) npanel !количество панелей в текущем участке
real(8) x1,y1,x2,y2  !координаты концов отрезка
real(8) ds1          !шаг в начале отрезка
real(8) ds2          !шаг в конце отрезка
type(TAreaValue) x,y
integer(4) n
integer(4) min_np
min_np=gs%const%di*2+4
call init_boundline_geomlineds2_array(ds1,ds2,x1,y1,x2,y2,x,y,npanel,min_np)
call add_boundline_geom(ibndl,npanel,x%v,y%v,n)
deallocate(x%v,y%v)
call init_geom_detale2(ibndl,n,npanel+n,x1,y1,x2,y2)
end

subroutine init_boundline_geomlineds2_array(ds1,ds2,x1,y1,x2,y2,vx,vy,npanel,min_np)
!инициализировать геометрию для участка границы в виде прямой линии с переменным ds
!меняющимся от ds1 до ds2
use pgmod
integer(4) npanel !количество панелей в текущем участке
real(8) x1,y1,x2,y2  !координаты концов отрезка
real(8) ds1          !шаг в начале отрезка
real(8) ds2          !шаг в конце отрезка
integer(4) min_np    !минимальное число панелей (если 0 - не учитывается)
logical use_min_np
real(8) ds1_,ds2_
integer(4) i,get_ss_npanel
real(8) smax,sred
real(8), allocatable :: ss(:)
real(8), pointer :: x(:),y(:)
type(TAreaValue), target :: vx,vy
real(8) DistanceToPoint
logical ds1gtds2
smax=DistanceToPoint(x1,y1,x2,y2)
ds1gtds2=ds1>ds2
if (ds1gtds2) then
  ds1_=ds2
  ds2_=ds1
else
  ds1_=ds1
  ds2_=ds2
endif
npanel=get_ss_npanel(ds1_,ds2_,smax)
use_min_np=npanel<min_np
if (use_min_np) npanel=min_np
allocate(ss(npanel+1),vx%v(npanel+1),vy%v(npanel+1))
x=>vx%v
y=>vy%v
if (use_min_np) then
  call sred_aray(d0,smax,ss,npanel+1)
else
  call get_ss_geom(ss,npanel,ds1_,ds2_,smax)
  if (ds1gtds2) then
    do i=1,npanel+1
      x(i)=smax-ss(i)
    enddo
    do i=1,npanel+1
      ss(i)=x(npanel+2-i)
    enddo
  endif
endif
do i=1,npanel+1
  x(i)=sred(d0,smax,x1,x2,ss(i))
  y(i)=sred(d0,smax,y1,y2,ss(i))
enddo
deallocate(ss)
end

subroutine pg_init_boundline_geomlineds3(ibndl,ds1,ds2,kk,mode,x1,y1,x2,y2,min_np)
!dec$ attributes dllexport:: pg_init_boundline_geomlineds3
!инициализировать геометрию для участка границы в виде прямой линии с переменным ds
!меняющимся от ds1 до ds2 в геометрической прогрессии с коэффициентом kk, затем с равномерным шагом ds2
use pgmod
integer(4) ibndl !номер участка границы
integer(4) npanel !количество панелей в текущем участке
integer(4) npanel_k  !количество панелей на сгущающемся участке
real(8) x1,y1,x2,y2  !координаты концов отрезка
real(8) ds1          !шаг в начале отрезка
real(8) ds2          !шаг в конце (или в середине) отрезка
real(8) kk           !коэффициент изменения длины (>1!!!)
integer(4) mode      !тип сгущения
					 !1 - -- --- ---- ---- ----
					 !2 ---- ---- ---- --- -- - 
					 !3 - -- --- --- --- --- -- -
integer(4) min_np    !минимальное число панелей
					 !если получается число панелей,меньше минимального, 
					 !то делается равномерное разбиение с минимальным числом панелей
type(TAreaValue) x,y
call init_boundline_geomlineds3_array(ds1,ds2,kk,mode,x1,y1,x2,y2,min_np,x,y,npanel,npanel_k)
call init_boundline_geom(ibndl,npanel,x%v,y%v,.false.)
deallocate(x%v,y%v)
call init_geom_detale2(ibndl,1,npanel+1,x1,y1,x2,y2)
end

subroutine pg_add_boundline_geomlineds3(ibndl,ds1,ds2,kk,mode,x1,y1,x2,y2,min_np)
!dec$ attributes dllexport:: pg_add_boundline_geomlineds3
!инициализировать геометрию для участка границы в виде прямой линии с переменным ds
!меняющимся от ds1 до ds2 в геометрической прогрессии с коэффициентом kk, затем с равномерным шагом ds2
!добавить к существующей
use pgmod
integer(4) ibndl !номер участка границы
integer(4) npanel !количество панелей в текущем участке
integer(4) npanel_k  !количество панелей на сгущающемся участке
real(8) x1,y1,x2,y2  !координаты концов отрезка
real(8) ds1          !шаг в начале отрезка
real(8) ds2          !шаг в конце (или в середине) отрезка
real(8) kk           !коэффициент изменения длины (>1!!!)
integer(4) mode      !тип сгущения
					 !1 - -- --- ---- ---- ----
					 !2 ---- ---- ---- --- -- - 
					 !3 - -- --- --- --- --- -- -
integer(4) min_np    !минимальное число панелей
					 !если получается число панелей,меньше минимального, 
					 !то делается равномерное разбиение с минимальным числом панелей
type(TAreaValue) x,y
integer(4) n
call init_boundline_geomlineds3_array(ds1,ds2,kk,mode,x1,y1,x2,y2,min_np,x,y,npanel,npanel_k)
call add_boundline_geom(ibndl,npanel,x%v,y%v,n)
deallocate(x%v,y%v)
call init_geom_detale2(ibndl,n,npanel+n,x1,y1,x2,y2)
end

subroutine init_boundline_geomlineds3_array(ds1,ds2,kk,mode,x1,y1,x2,y2,min_np,x,y,npanel,npanel_k)
use pgmod
integer(4) npanel !количество панелей в текущем участке
integer(4) npanel_k  !количество панелей на сгущающемся участке
real(8) x1,y1,x2,y2  !координаты концов отрезка
real(8) ds1          !шаг в начале отрезка
real(8) ds2          !шаг в конце (или в середине) отрезка
real(8) kk           !коэффициент изменения длины (>1!!!)
integer(4) mode      !тип сгущения
					 !1 - -- --- ---- ---- ----
					 !2 ---- ---- ---- --- -- - 
					 !3 - -- --- --- --- --- -- - (получается гарантированное четное число панелей)
integer(4) min_np    !минимальное число панелей
					 !если получается число панелей,меньше минимального, 
					 !то делается равномерное разбиение с минимальным числом панелей
type(TAreaValue) x,y
if (mode==2) then
  call init_boundline_geomlineds3_array_mode1(ds2,ds1,kk,1,x2,y2,x1,y1,min_np,x,y,npanel,npanel_k)
  call swap_array(x%v,npanel+1)
  call swap_array(y%v,npanel+1)
else
  call init_boundline_geomlineds3_array_mode1(ds1,ds2,kk,mode,x1,y1,x2,y2,min_np,x,y,npanel,npanel_k)
endif
end 

subroutine init_boundline_geomlineds3_array_mode1(ds1,ds2,kk,mode,x1,y1,x2,y2,min_np,vx,vy,npanel,npanel_k)
use pgmod
integer(4) npanel !количество панелей в текущем участке
integer(4) npanel_k  !количество панелей на сгущающемся участке
real(8) x1,y1,x2,y2  !координаты концов отрезка
real(8) ds1          !шаг в начале отрезка
real(8) ds2          !шаг в конце (или в середине) отрезка
real(8) kk           !коэффициент изменения длины (>1!!!)
integer(4) mode      !тип сгущения
					 !1 - -- --- ---- ---- ----
					 !3 - -- --- --- --- --- -- - (получается гарантированное четное число панелей)
integer(4) min_np    !минимальное число панелей
					 !если получается число панелей,меньше минимального, 
					 !то делается равномерное разбиение с минимальным числом панелей
real(8) dsr          !шаг на равномерном участке
real(8) s,s1,kk2,kk_,kk1
integer(4) i,n,n1,mnp
real(8) smax,ds,sred,dsmin
type(TAreaValue), target :: vx,vy
real(8), allocatable :: ss(:)
real(8), pointer :: x(:),y(:)
real(8) DistanceToPoint
logical use_min_np
kk1=d1
kk2=d1
smax=DistanceToPoint(x1,y1,x2,y2)
dsmin=d0
if (min_np>0) dsmin=smax/min_np
kk_=kk
if (ds1>ds2) kk_=d1/kk
dsr=ds2
if (mode==3) smax=smax*d5
n=(dlog(ds2/ds1)/dlog(kk_)+d5) !количество панелей на сгущающемся участке
n1=0
s=ds1*(kk_**n-d1)/(kk_-d1) !длина сгущающегося учатка
if (smax>s) then
  n1=((smax-s)/dsr+d5) !количество панелей на равномерном участке
else
  n=(dlog(smax*(kk_-d1)/ds1+d1)/dlog(kk_)+d5)
  s=ds1*(kk_**n-d1)/(kk_-d1) !длина сгущающегося учатка
endif
npanel=n+n1
mnp=min_np
if (mode==3) mnp=(min_np+1)/2 !для нечетного +1
do while (npanel<mnp.and.n>0)
  if (kk_>d1) then
    n=n-1
  else
    n=n+1
  endif
  s=d0
  if (n>0) s=ds1*(kk_**n-d1)/(kk_-d1) !длина сгущающегося участка
  dsr=ds1*kk_**n
  if (dsr<dsmin) then
    npanel=0
    exit
  endif
  n1=((smax-s)/dsr+d5) !количество панелей на равномерном участке
  if (n1<0) n1=0
  npanel=n+n1
enddo
npanel_k=n
use_min_np=npanel<mnp
if (use_min_np) npanel=mnp
if (mode==3) npanel=npanel*2
allocate(ss(npanel+1))
allocate(vx%v(npanel+1),vy%v(npanel+1))
x=>vx%v
y=>vy%v
ss=d0
if (use_min_np) then
  ss=smax/npanel
  ss(1)=d0
else
  s=d0
  if (n>0) s=ds1*(kk_**n-d1)/(kk_-d1) !длина сгущающегося участка
  s1=dsr*n1   !длина равномерного участка
  if (s>d0) then
    kk2=(smax-s1)/s   !коэффициент растяжения/сжатия сгущающегося участка
  elseif (s1>d0) then
    kk1=smax/s1   !коэффициент растяжения/сжатия равномерного участка
  else 
    n=0
	  n1=0
  endif
  ds=ds1*kk2
  do i=1,n
    ss(i+1)=ds
    ds=ds*kk_
  enddo
  do i=1,n1
    ss(i+n+1)=dsr*kk1
  enddo
  if (mode==3) then
    n=npanel/2
    do i=1,n
      ss(i+n+1)=ss(n-i+2)
    enddo
    smax=smax*d2
  endif
endif
do i=1,npanel+1
  if (i>1) ss(i)=ss(i-1)+ss(i)
  x(i)=sred(d0,smax,x1,x2,ss(i))
  y(i)=sred(d0,smax,y1,y2,ss(i))
enddo
deallocate(ss)
end

subroutine pg_init_boundline_geomcircle(ibndl,npanel,x0,y0,r,gam1,gam2)
!dec$ attributes dllexport:: pg_init_boundline_geomcircle
!инициализировать геометрию для участка границы в виде дуги окружности
use gen_mod
integer(4) ibndl !номер участка границы
integer(4) npanel !количество панелей в текущем участке
real(8) x0,y0       !координаты центра круга
real(8) r           !радиус
real(8) gam1,gam2   !углы начала и конца дуги
call init_boundline_geom_ellipse(ibndl,npanel,x0,y0,r,r,d0,gam1,gam2)
call init_geom_detale3(ibndl,1,npanel+1,x0,y0,gam1,gam2,r)
end

subroutine init_geom_detale3(ibndl,i_begin,i_end,x,y,gam1,gam2,r)
use pgmod
type(TBoundline_geomdetale), pointer :: gd
integer(4) ibndl,i_begin,i_end
real(8) x,y,gam1,gam2,r
call allocate_geom_detale_for_next
gd=>gsbnd%geom_detale0(gsbnd%ngeom_detale0)
gd%mode=3
gd%ibndl=ibndl
gd%i_begin=i_begin
gd%i_end=i_end
gd%x=x
gd%y=y
gd%gam1=gam1
gd%gam2=gam2
gd%r=r
call init_geom_detale_bndl(ibndl)
end

subroutine pg_init_boundline_geom_ellipse(ibndl,npanel,x0,y0,rx,ry,g0,gam1,gam2)
!dec$ attributes dllexport:: pg_init_boundline_geom_ellipse
!инициализировать геометрию для участка границы в виде дуги окружности
use gen_mod
integer(4) ibndl !номер участка границы
integer(4) npanel !количество панелей в текущем участке
real(8) x0,y0       !координаты центра круга
real(8) rx          !полуось по оси x
real(8) ry          !полуось по оси y
real(8) g0          !направление первой полуоси (x)
real(8) gam1,gam2   !углы начала и конца дуги (отсчитываемые от g0)
call init_boundline_geom_ellipse(ibndl,npanel,x0,y0,rx,ry,g0,gam1,gam2)
if (rx==ry) then
  call init_geom_detale3(ibndl,1,npanel+1,x0,y0,gam1,gam2,rx)
else
  call init_geom_detale1(ibndl,1,npanel+1)
endif
end

subroutine init_boundline_geom_ellipse(ibndl,npanel,x0,y0,rx,ry,g0,gam1,gam2)
!инициализировать геометрию для участка границы в виде дуги окружности
use gen_mod
integer(4) ibndl !номер участка границы
integer(4) npanel !количество панелей в текущем участке
real(8) x0,y0       !координаты центра круга
real(8) rx          !полуось по оси x
real(8) ry          !полуось по оси y
real(8) g0          !направление первой полуоси (x)
real(8) gam1,gam2   !углы начала и конца дуги (отсчитываемые от g0)
real(8) x(npanel+1) !x - концов панелей
real(8) y(npanel+1) !y - концов панелей
integer(4) i
real(8) g,dg
dg=(gam2-gam1)/npanel
do i=1,npanel+1
  g=gam1+(i-d1)*dg+g0
  x(i)=x0+rx*dcos(g)
  y(i)=y0+ry*dsin(g)
enddo
call init_boundline_geom(ibndl,npanel,x,y,.false.)
end

subroutine pg_init_boundline_geom_ellipse_ds(ibndl,ds,x0,y0,rx,ry,g0,gam1,gam2)
!dec$ attributes dllexport:: pg_init_boundline_geom_ellipse_ds
!инициализировать геометрию для участка границы в виде дуги окружности с фиксированным ds
use gen_mod
integer(4) ibndl !номер участка границы
integer(4) npanel !количество панелей в текущем участке
real(8) ds          !примерный шаг
real(8) x0,y0       !координаты центра круга
real(8) rx          !полуось по оси x
real(8) ry          !полуось по оси y
real(8) g0          !направление первой полуоси (x)
real(8) gam1,gam2   !углы начала и конца дуги (отсчитываемые от g0)
real(8) s,eps,r
eps=1.0d-6
r=rx
if (rx<ry) r=ry
s=dabs(gam1-gam2)*r+eps
npanel=s/ds
call pg_init_boundline_geom_ellipse(ibndl,npanel,x0,y0,rx,ry,g0,gam1,gam2)
end

subroutine pg_init_boundline_geomcircleds(ibndl,ds,x0,y0,r,gam1,gam2)
!dec$ attributes dllexport:: pg_init_boundline_geomcircleds
!инициализировать геометрию для участка границы в виде дуги окружности с фиксированным ds
use gen_mod
integer(4) ibndl !номер участка границы
integer(4) npanel !количество панелей в текущем участке
real(8) ds          !примерный шаг
real(8) x0,y0       !координаты центра круга
real(8) r           !радиус
real(8) gam1,gam2   !углы начала и конца дуги
real(8) s,eps
eps=1.0d-6
s=dabs(gam1-gam2)*r+eps
npanel=s/ds
call pg_init_boundline_geomcircle(ibndl,npanel,x0,y0,r,gam1,gam2)
end

subroutine pg_copy_boundline(ibndl,ip2,ia2,ibnd2,ibndl2)
!dec$ attributes dllexport:: pg_copy_boundline
!скопировать участок с границы другой области
!с инверсией обхода - случай общей границы между двумя областями
use pgmod
integer(4) ibndl !номер участка границы
integer(4) ip2,ia2,ibnd2,ibndl2  !номер задачи, области, границы и участка, откуда будет произведено копирование
call copy_boundline_dxdy_rotate(ibndl,ip2,ia2,ibnd2,ibndl2,d0,d0,d0,d0,d0,.true.)
end

subroutine pg_copy_boundline_dxdy(ibndl,ip2,ia2,ibnd2,ibndl2,dx,dy)
!dec$ attributes dllexport:: pg_copy_boundline_dxdy
!скопировать участок с границы другой области
!со смещением
use pgmod
integer(4) ibndl !номер участка границы
integer(4) ip2,ia2,ibnd2,ibndl2  !номер задачи, области, границы и участка, откуда будет произведено копирование
real(8) dx,dy !смещение
call copy_boundline_dxdy_rotate(ibndl,ip2,ia2,ibnd2,ibndl2,dx,dy,d0,d0,d0,.false.)
end

subroutine pg_copy_boundline_dxdy_rotate(ibndl,ip2,ia2,ibnd2,ibndl2,dx,dy,angle,xr,yr)
!dec$ attributes dllexport:: pg_copy_boundline_dxdy_rotate
!скопировать участок с границы другой области
!со смещением
use pgmod
integer(4) ibndl !номер участка границы
integer(4) ip2,ia2,ibnd2,ibndl2  !номер задачи, области, границы и участка, откуда будет произведено копирование
real(8) dx,dy !смещение
real(8) angle  !угол поворота против часовой стрелки (в радианах)
real(8) xr,yr  !точка, вокруг которой осуществляется поворот
call copy_boundline_dxdy_rotate(ibndl,ip2,ia2,ibnd2,ibndl2,dx,dy,angle,xr,yr,.false.)
end

subroutine pg_copy_bound_dxdy(ip2,ia2,ibnd2,dx,dy)
!dec$ attributes dllexport:: pg_copy_bound_dxdy
!скопировать границу другой области
!со смещением
use pgmod
integer(4) ip2,ia2,ibnd2  !номер задачи, области и границы, откуда будет произведено копирование
real(8) dx,dy !смещение
call pg_copy_bound_dxdy_rotate(ip2,ia2,ibnd2,dx,dy,d0,d0,d0)
end

subroutine pg_copy_bound_dxdy_rotate(ip2,ia2,ibnd2,dx,dy,angle,xr,yr)
!dec$ attributes dllexport:: pg_copy_bound_dxdy_rotate
!скопировать границу другой области
!со смещением
use pgmod
integer(4) ip2,ia2,ibnd2  !номер задачи, области и границы, откуда будет произведено копирование
real(8) dx,dy !смещение
real(8) angle  !угол поворота против часовой стрелки (в радианах)
real(8) xr,yr  !точка, вокруг которой осуществляется поворот
integer(4) i
type(TBound), pointer :: b
b=>gsMain%ggs(ip2)%a(ia2)%bnd(ibnd2)
call pg_allocate_boundlines(b%nline)
do i=1,b%nline
  call copy_boundline_dxdy_rotate(i,ip2,ia2,ibnd2,i,dx,dy,angle,xr,yr,.false.)
enddo
end

subroutine pg_copy_allbound_dxdy(ip2,ia2,dx,dy)
!dec$ attributes dllexport:: pg_copy_allbound_dxdy
!скопировать все границы другой области
!со смещением
use pgmod
integer(4) ip2,ia2  !номер задачи, области, откуда будет произведено копирование
real(8) dx,dy !смещение
call pg_copy_allbound_dxdy_rotate(ip2,ia2,dx,dy,d0,d0,d0)
end

subroutine pg_copy_allbound_dxdy_rotate(ip2,ia2,dx,dy,angle,xr,yr)
!dec$ attributes dllexport:: pg_copy_allbound_dxdy_rotate
!скопировать все границы другой области
!со смещением
use pgmod
integer(4) ip2,ia2  !номер задачи, области, откуда будет произведено копирование
real(8) dx,dy !смещение
real(8) angle  !угол поворота против часовой стрелки (в радианах)
real(8) xr,yr  !точка, вокруг которой осуществляется поворот
integer(4) i
type(TArea), pointer :: a
a=>gsMain%ggs(ip2)%a(ia2)
call pg_allocate_bounds(a%nb)
do i=1,a%nb
  call pg_bind_bound(i)
  call pg_copy_bound_dxdy_rotate(ip2,ia2,i,dx,dy,angle,xr,yr)
enddo
end

subroutine copy_boundline_dxdy_rotate(ibndl,ip2,ia2,ibnd2,ibndl2,dx,dy,angle,xr,yr,invert)
!скопировать участок с границы другой области
!со смещением и инверсией
use pgmod
integer(4) ibndl !номер участка границы
integer(4) ip2,ia2,ibnd2,ibndl2  !номер задачи, области, границы и участка, откуда будет произведено копирование
real(8) dx,dy !смещение
logical invert !инвертировать обход
real(8) angle  !угол поворота против часовой стрелки (в радианах)
real(8) xr,yr  !точка, вокруг которой осуществляется поворот
integer(4) i,ib,ie
integer(4) i1,i2,istep
type(TBound), pointer :: b
type(TBoundline), pointer :: bl
type(TBoundline_geomdetale), pointer :: gd
real(8), allocatable :: xx(:), yy(:)
complex(8) z,z0,ett,z2,dz
b=>gsMain%ggs(ip2)%a(ia2)%bnd(ibnd2)
if (b%npanel==0) call gs_print_stop("Error copy_boundline_dxdy_rotate. No geom_postprocessor")
bl=>b%line(ibndl2)
if (dx==d0.and.dy==d0.and.angle==d0) then
  call init_boundline_geom(ibndl,bl%npanel,bl%x,bl%y,invert)
else
  allocate(xx(bl%npanel+1),yy(bl%npanel+1))
  if (angle==d0) then
    xx=bl%x+dx
    yy=bl%y+dy
  else
    z0=dcmplx(xr,yr)
    ett=cdexp(ii*angle)
    dz=dcmplx(dx,dy)
    do i=1,bl%npanel+1
      z=dcmplx(bl%x(i),bl%y(i))
      z=(z-z0)*ett+z0+dz
      xx(i)=dreal(z)
      yy(i)=dimag(z)
    enddo
  endif
  call init_boundline_geom(ibndl,bl%npanel,xx,yy,invert)
  deallocate(xx,yy)
endif
i1=bl%igd_begin
i2=bl%igd_end
istep=1
if (invert) then
  call swap_int(i1,i2)
  istep=-1
endif
do i=i1,i2,istep
  gd=>b%geom_detale0(i)
  ib=gd%i_begin-bl%i_begin+1
  ie=gd%i_end-bl%i_begin+1
  if (invert) then
    call swap_int(ib,ie)
    ib=bl%npanel-ib+2
    ie=bl%npanel-ie+2
  endif
  if (gd%mode==2.or.gd%mode==3) then
    if (angle==d0) then
      z=dcmplx(gd%x+dx,gd%y+dy)
      if (gd%mode==2) z2=dcmplx(gd%x2+dx,gd%y2+dy)
    else
      z=dcmplx(gd%x,gd%y)
      z=(z-z0)*ett+z0+dz
      if (gd%mode==2) then
        z2=dcmplx(gd%x2,gd%y2)
        z2=(z2-z0)*ett+z0+dz
      endif
    endif
  endif
  select case(gd%mode)
  case(1)
    call init_geom_detale1(ibndl,ib,ie)
  case(2)
    if (invert) then
      call init_geom_detale2(ibndl,ib,ie,dreal(z2),dimag(z2),dreal(z),dimag(z))
    else
      call init_geom_detale2(ibndl,ib,ie,dreal(z),dimag(z),dreal(z2),dimag(z2))
    endif
  case(3)
    if (invert) then
      call init_geom_detale3(ibndl,ib,ie,dreal(z),dimag(z),gd%gam2+angle,gd%gam1+angle,gd%r)
    else
      call init_geom_detale3(ibndl,ib,ie,dreal(z),dimag(z),gd%gam1+angle,gd%gam2+angle,gd%r)
    endif
  end select
enddo
call pg_set_boundline_gu_mode(ibndl,bl%gu_mode)
end

subroutine pg_copy_subdomain_dxdy_rotate(ia1,ia2,dx,dy,angle,xr,yr)
!dec$ attributes dllexport:: pg_copy_subdomain_dxdy_rotate
!создать область как копию другой подобласти, включая тип уравнения, границу, сетку в области
use pgmod
integer(4) ia1  !номер области, откуда будет произведено копирование
integer(4) ia2  !номер области, куда будет произведено копирование
real(8) dx,dy !смещение
real(8) angle  !угол поворота против часовой стрелки (в радианах)
real(8) xr,yr  !точка, вокруг которой осуществляется поворот
type(TArea), pointer :: a1
a1=>gs%a(ia1)
call pg_bind_domain(ia2)
select case (a1%type_eq(1))
case (8,11)
  call pg_set_domain_equation_syst(a1%type_eq(1),a1%type_eq(2))
case (22:26)
  call pg_set_domain_equation(a1%type_eq(1))  
  call pg_set_eq_var(a1%eq_var,a1%n_eq_var)
case default
  call pg_set_domain_equation(a1%type_eq(1))  
end select
call pg_set_areaconst(1,a1%const%k_helm)
call pg_set_areaconst(2,a1%const%k_oss)
call pg_set_areaconst(3,a1%const%Re)
call pg_copy_geom_domain_dxdy_rotate(gs%i,ia1,dx,dy,angle,xr,yr)
if (a1%a%npart>0) call pg_allocate_area_gu
end

subroutine pg_copy_geom_domain_dxdy_rotate(ip,ia,dx,dy,angle,xr,yr)
!dec$ attributes dllexport:: pg_copy_geom_domain_dxdy_rotate
!скопировать геометрию другой области (границу, сетку в области)
!задача и область должны быть текущими
!тип уравнения в текущей области должен быть задан!
use pgmod
integer(4) ia  !номер области, откуда будет произведено копирование
integer(4) ip  !номер задачи, откуда будет произведено копирование
real(8) dx,dy !смещение
real(8) angle  !угол поворота против часовой стрелки (в радианах)
real(8) xr,yr  !точка, вокруг которой осуществляется поворот
real(8), parameter :: eps=1.0d-8
type(TArea), pointer :: a1
integer(4) i
logical eq_eps
a1=>gsMain%ggs(ip)%a(ia)
call pg_copy_allbound_dxdy_rotate(ip,ia,dx,dy,angle,xr,yr)
call pg_geom_postprocessor
if (a1%a%npart>0) then
  call pg_allocate_area(a1%a%npart)
  do i=1,a1%a%npart
    call pg_bind_areapart(i)
    call pg_init_areapart_geom_ref_dxdy_rotate(ip,ia,i,dx,dy,angle,xr,yr)
  enddo
  call pg_areageom_postprocessor
endif
gsarea%dx=dx
gsarea%dy=dy
gsarea%cosa=dcos(angle)
gsarea%sina=dsin(angle)
gsarea%a_ref=>a1
gsarea%type_rotate=4
if (eq_eps(gsarea%sina,d0,eps)) then
  if (gsarea%cosa>eps) then
    gsarea%type_rotate=0
  else
    gsarea%type_rotate=1
  endif
elseif (eq_eps(gsarea%cosa,d0,eps)) then
  if (gsarea%sina>eps) then
    gsarea%type_rotate=2
  else
    gsarea%type_rotate=3
  endif
endif
if (gsarea%type_eq(1)==20) gsarea%oss_resolve=gsarea%type_rotate>0.or.gsarea%dy/=d0
end

subroutine pg_geom_postprocessor
!dec$ attributes dllexport:: pg_geom_postprocessor
!инициализация центров панелей и доп массивов
use pgmod
integer(4) i,k,p,m
complex(8) dz
type(TBoundline), pointer :: bl
type(TBoundline2), pointer :: bl2
type(TBoundline_Collocate), pointer :: gc
type(TBound), pointer :: b
real(8) g
do k=1,gsarea%nb
  b=>gsarea%bnd(k)
  select case (b%boundLineType)
  case (1)
    m=0
    do p=1,b%nline
      bl=>b%line(p)
      m=m+bl%npanel
    enddo
    call allocate_bound_geomgu(b,m)
    m=0
    do p=1,b%nline
      bl=>b%line(p)
	    bl%i_begin=m+1
	    do i=1,bl%npanel
	      m=m+1
	      b%x(m)=bl%x(i)
    	  b%y(m)=bl%y(i)
      enddo
      bl%i_end=m
    enddo
    m=m+1
    !b%x(m)=b%x(1)
    !b%y(m)=b%y(1)
    i=bl%npanel+1  !для незамкнутой границы
    b%x(m)=bl%x(i)
    b%y(m)=bl%y(i)
    do i=1,b%npanel+1
      b%z(i)=dcmplx(b%x(i),b%y(i))
    enddo
    b%s(1)=d0
    do i=1,b%npanel
      b%xc(i)=(b%x(i)+b%x(i+1))*d5
      b%yc(i)=(b%y(i)+b%y(i+1))*d5
      b%zc(i)=dcmplx(b%xc(i),b%yc(i))
      dz=b%z(i+1)-b%z(i)
      b%l(i)=cdabs(dz)
      b%ett(i)=dz/b%l(i)
      b%s(i+1)=b%s(i)+b%l(i)
      b%sc(i)=(b%s(i+1)+b%s(i))*d5
    enddo
    call geom_detale_postprocessor(k)
  case (2)
    do p=1,b%nline
      bl2=>b%line2(p)
	    if (.not.allocated(bl2%gc)) cycle !!!не нужно, если это решение, считанное из файла
      do i=1,gsarea%nu
        gc=>bl2%gc(i)
        do m=1,gc%n
		      gc%g(m)=pi2*gc%s_gu(m)
          g=bl2%g0+gc%g(m)*bl2%dir
          gc%ztc(m)=-bl2%r*cdexp(ii*g)
          gc%z(m)=bl2%zc-gc%ztc(m)
        enddo
      enddo
    enddo
  end select
enddo
end

subroutine geom_detale_postprocessor(ib)
use pgmod
integer(4), parameter :: nmode=3
integer(4) ib,i,j,k
real(8) v,vect_p2,s,scal_p2
type(TBound), pointer :: b
type(TBoundline), pointer :: bl
type(TBoundline_geomdetale), pointer :: gd,gd2
type(TBoundline_geomdetale), allocatable,target :: gdtemp(:)
integer(4) mode_sort(nmode)
data mode_sort /2,3,1/
b=>gsarea%bnd(ib)
if (b%ngeom_detale0==0) then
  gd=>b%geom_detale0(1)
  gd%mode=1
  gd%i_begin=1
  gd%i_end=b%npanel+1
  b%ngeom_detale0=1
  k=1
else
  do i=1,b%ngeom_detale0
    gd=>b%geom_detale0(i)
    bl=>b%line(gd%ibndl)
    k=bl%i_begin-1
    gd%i_begin=gd%i_begin+k
    gd%i_end=gd%i_end+k
  enddo
endif
allocate(gdtemp(b%ngeom_detale0))
gd=>gdtemp(1)
gd=b%geom_detale0(1)
k=1
do i=2,b%ngeom_detale0
  gd2=>b%geom_detale0(i)
  if (gd%mode==gd2%mode) then
    select case (gd%mode)
    case (1)
      gd%i_end=gd2%i_end
      cycle
    case (2)
      v=vect_p2(gd%x2-gd%x,gd%y2-gd%y,gd2%x2-gd2%x,gd2%y2-gd2%y)
      s=scal_p2(gd%x2-gd%x,gd%y2-gd%y,gd2%x2-gd2%x,gd2%y2-gd2%y)
      if (dabs(v)<1.0d-8.and.s>d0) then
        gd%i_end=gd2%i_end
        gd%x2=gd2%x2
        gd%y2=gd2%y2
        cycle
      endif
    endselect
  endif
  k=k+1
  gd=>gdtemp(k)
  gd=gd2
enddo
b%ngeom_detale=k
allocate(b%geom_detale(k))
k=0
do j=1,nmode
  do i=1,b%ngeom_detale
    gd=>gdtemp(i)
    if (gd%mode.ne.mode_sort(j)) cycle
    k=k+1
    b%geom_detale(k)=gd
  enddo
enddo
deallocate(gdtemp)
do i=1,b%ngeom_detale
  gd=>b%geom_detale(i)
  gd%s1=b%s(gd%i_begin)
  gd%s2=b%s(gd%i_end)
  gd%panelLmax=maxval(b%l(gd%i_begin:gd%i_end-1))
enddo
end

subroutine geom_collocate_points
use pgmod
integer(4) i,k,p,j,i0
integer(4) ncp,ncpl
type(TArea), pointer :: a
type(TBoundline), pointer :: bl
type(TBoundline2), pointer :: bl2
type(TBound_Collocate_points), pointer :: cp
type(TBound), pointer :: b
!на границе
do i0=1,gs%na
  a=>gs%a(i0)
  cp=>a%cpp
  if (cp%inited) cycle
  ncp=0
  do k=1,a%nb
    b=>a%bnd(k)
    select case (b%boundLineType)
    case (1)
  	  do p=1,b%nline
        bl=>b%line(p)
  		  do i=1,a%nu  
          ncpl=0
          if (bl%cp_line(i)) then
  		      if (allocated(bl%use_gu)) then
  		        do j=1,bl%npanel
    	          if (bl%use_gu(j,i)) ncpl=ncpl+1
              enddo
            else
  		        ncpl=bl%npanel
            endif
  		    endif
          ncp=ncp+ncpl
        enddo
  	  enddo
    case (2)
  	  do p=1,b%nline
        bl2=>b%line2(p)
  	    do i=1,a%nu
  	      ncp=ncp+bl2%gc(i)%n
  	    enddo
  	  enddo
    end select
  enddo
  !точки коллокации
  call allocate_bound_collocate_points(i0,ncp)
  ncpl=0 !индекс текущей точки коллокации
  do i=1,a%nu
    do k=1,a%nb
      b=>a%bnd(k)
      select case (b%boundLineType)
      case (1)
        do p=1,b%nline
          bl=>b%line(p)
          if (.not.bl%cp_line(i)) cycle
          do j=1,bl%npanel
            if (allocated(bl%use_gu)) then
              if (.not.bl%use_gu(j,i)) cycle
            endif
            ncpl=ncpl+1
            cp%i(ncpl)=j+bl%i_begin-1
            cp%ibnd(ncpl)=k
  	        cp%iu(ncpl)=i
          enddo
        enddo
      case (2)
        do p=1,b%nline
          bl2=>b%line2(p)
  	      do j=1,bl2%gc(i)%n
            ncpl=ncpl+1
            cp%i(ncpl)=j
            cp%ibndl(ncpl)=p
            cp%ibnd(ncpl)=k
            cp%iu(ncpl)=i
          enddo
        enddo
      end select
    enddo
  enddo
  cp%inited=.true.
enddo
end

subroutine geom_collocate_points_area
use pgmod
integer(4) i,k,j,i0
type(TArea), pointer :: a
type(areatype), pointer :: aa
!в области
do i0=1,gs%na
  a=>gs%a(i0)
  aa=>a%a
  if (a%haveAreaEq.and.(.not.aa%cppa%inited)) then
    !call allocate_area_collocate_points(i0,ncp) перенесли в get_ind_area!
    k=0
    do i=1,aa%ntr
      do j=1,a%umaxtr
        if (aa%areaind(i,j)==0) cycle
        k=k+1
  	    aa%cppa%iu(k)=j
  	    aa%cppa%itr(k)=i
  	  enddo
    enddo
    aa%cppa%inited=.true.
  endif
enddo
end

subroutine pg_areageom_postprocessor
!dec$ attributes dllexport:: pg_areageom_postprocessor
!инициализация центров треугольников
use pgmod
integer(4) i,j,k,n_all,ntr_all,npe_all
complex(8) zmz
type(areapart), pointer :: a
type(areatype), pointer :: aa
aa=>gsarea%a
n_all=0
ntr_all=0
npe_all=0
do i=1,aa%npart
  a=>aa%part(i)
  n_all=n_all+a%n
  ntr_all=ntr_all+a%ntr
  if (npe_all<a%npe) npe_all=a%npe
enddo
call allocate_area_geom(n_all,ntr_all,npe_all)
j=0 !n
k=0 !ntr
do i=1,aa%npart
  a=>aa%part(i)
  a%n_begin=j+1
  a%n_end=j+a%n
  a%ntr_begin=k+1
  a%ntr_end=k+a%ntr
  if (associated(a%apref)) then
    aa%zm(a%n_begin:a%n_end)=a%a_ref*a%apref%zm(:)+a%b_ref*dconjg(a%apref%zm(:))+a%c_ref
    aa%trm(1:a%npe,a%ntr_begin:a%ntr_end)=a%apref%trm(:,:)+j
    aa%npe_ar(a%ntr_begin:a%ntr_end)=a%apref%npe
    call test_invert_trm(i)
  else
    aa%zm(a%n_begin:a%n_end)=a%zm(:)
    aa%trm(1:a%npe,a%ntr_begin:a%ntr_end)=a%trm(:,:)+j
    aa%npe_ar(a%ntr_begin:a%ntr_end)=a%npe
  endif
  !call deallocate_areapart(gsarea%i,i)
  j=a%n_end
  k=a%ntr_end
enddo
aa%geom_inited=.true.
do i=1,aa%ntr
  do k=1,aa%npe_ar(i)
    aa%zmc(i)=aa%zmc(i)+zmz(k,i,gsarea%i)
  enddo
  aa%zmc(i)=aa%zmc(i)/aa%npe_ar(i)
enddo
end

subroutine test_invert_trm(ipart)
use pgmod
integer(4) ipart,i,j,k,t
complex(8) z1,z2,z3,v1,v2
real(8) s,st,vect_p
type(areapart), pointer :: a
type(areatype), pointer :: aa
aa=>gsarea%a
a=>aa%part(ipart)
s=d0
z1=aa%zm(aa%trm(1,a%ntr_begin))
do i=1,a%npe-2
  z2=aa%zm(aa%trm(i+1,a%ntr_begin))
  z3=aa%zm(aa%trm(i+2,a%ntr_begin))
  v1=z2-z1
  v2=z3-z1
  st=vect_p(v1,v2)
  s=s+st
enddo
if (s<d0) then
  do j=a%ntr_begin,a%ntr_end
    do i=2,(a%npe-1)/2+1
      k=a%npe-i+2
      t=aa%trm(i,j)
      aa%trm(i,j)=aa%trm(k,j)
      aa%trm(k,j)=t
    enddo
  enddo
endif
end

function get_ss_npanel(ds1,ds2,smax)
use gen_mod
integer(4) get_ss_npanel
real(8) ds1,ds2,smax
get_ss_npanel=d2*smax/(ds2+ds1)
end

subroutine get_ss_geom(ss,npanel,ds1,ds2,smax)
use gen_mod
real(8) ds1,ds2,k,smax
integer(4) i,npanel
real(8) ss(npanel+1) 
!шаг по радиусу растет от ds1 до ds2 (учитывается в get_ss_npanel)
!первый шаг в точночти равен ds1
k=(d2*smax/npanel-d2*ds1)/(npanel-d1)
ss(1)=d0
do i=1,npanel
  ss(i+1)=ds1+k*(i-d1)
enddo
if (npanel>2) then
  k=(smax-ds2-ds1)/(smax-ss(2)-ss(npanel+1))
  ss(2)=ds1
  ss(npanel+1)=ds2
  do i=3,npanel
    ss(i)=ss(i)*k
  enddo
endif 
do i=2,npanel+1
  ss(i)=ss(i-1)+ss(i)
enddo
end

function zmz(i,j,knd)
!координата вершины треугольника
use pgmod
complex(8) zmz
integer(4) i   !номер вершины (1,2,3,<4>)
integer(4) j   !номер треугольника
integer(4) knd !номер области
type(areatype), pointer :: a
a=>gs%a(knd)%a
zmz=a%zm(a%trm(i,j))
end

function zmz2(i,j)
!координата вершины треугольника
use pgmod
complex(8) zmz2
integer(4) i   !номер вершины (1,2,3,<4>)
integer(4) j   !номер треугольника
zmz2=gsarea%a%zm(gsarea%a%trm(i,j))
end

function xmc(i,knd)
!x центра треугольника
use pgmod
real(8) xmc
integer(4) i   !номер треугольника
integer(4) knd !номер области
type(areatype), pointer :: a
a=>gs%a(knd)%a
xmc=dreal(a%zmc(i))
end

function ymc(i,knd)
!y центра треугольника
use pgmod
real(8) ymc
integer(4) i   !номер треугольника
integer(4) knd !номер области
type(areatype), pointer :: a
a=>gs%a(knd)%a
ymc=dimag(a%zmc(i))
end

function elemarea(i)
!площадь элемента
use pgmod
integer(4) i   !номер треугольника
real(8) elemarea,s,vect_p
integer(4) j
complex(8) v1,v2
s=d0
do j=3,gsarea%a%npe_ar(i)
  v1=gsarea%a%zm(gsarea%a%trm(j-1,i))-gsarea%a%zm(gsarea%a%trm(1,i))
  v2=gsarea%a%zm(gsarea%a%trm(j,i))-gsarea%a%zm(gsarea%a%trm(1,i))
  s=s+vect_p(v1,v2)
enddo
elemarea=s*d5
end

subroutine pg_init_cylider_particle(ibndl,x0,y0,r,gam0,dir,type_ga)
!dec$ attributes dllexport:: pg_init_cylider_particle
!инициализировать стандартную геометрию для участка границы в круга для boundLineType=2
use pgmod
integer(4) ibndl !номер участка границы
real(8) x0,y0       !координаты центра круга
real(8) r           !радиус
real(8) gam0        !смещение точки с s=0 от оси x против часовой стрелки
integer(4) type_ga
integer(4) dir
type(TBoundline2), pointer :: b
b=>gsbnd%line2(ibndl)
b%i=ibndl
b%zc=dcmplx(x0,y0)
b%r=r
b%g0=gam0
b%dir=dir
b%type_ga=type_ga
call allocate_boundline_geom2(ibndl)
end

subroutine pg_init_areapart_geom_ref(ip,ia,ipart,ax,bx,cx,ay,by,cy)
!dec$ attributes dllexport:: pg_init_areapart_geom_ref
!инициализировать участок ссылкой на другой участок и преобразованием координат
!x1=ax*x+bx*y+cx
!y1=ay*x+by*y+cy
!z1=A*z+B*conj(z)+C
use pgmod
integer(4) ip !номер задачи
integer(4) ia !номер области
integer(4) ipart !номер участка
real(8) ax,bx,cx,ay,by,cy
complex(8) aa,bb,cc
type(areapart), pointer :: aref
aa=dcmplx((ax+by)*d5,(ay-bx)*d5)
bb=dcmplx((ax-by)*d5,(ay+bx)*d5)
cc=dcmplx(cx,cy)
aref=>gsMain%ggs(ip)%a(ia)%a%part(ipart)
call init_areapart_geom_ref(aref,aa,bb,cc)
end

recursive subroutine init_areapart_geom_ref(aref,aa,bb,cc)
!инициализировать участок ссылкой на другой участок и преобразованием координат
!z1=A*z+B*conj(z)+C
use pgmod
complex(8) aa,bb,cc
complex(8) aa2,bb2,cc2
complex(8) aa1,bb1,cc1
type(areapart), target :: aref
if (associated(aref%apref)) then
  aa1=aref%a_ref
  bb1=aref%b_ref
  cc1=aref%c_ref
  aa2=aa*aa1+bb*dconjg(bb1)
  bb2=aa*bb1+bb*dconjg(aa1)
  cc2=aa*cc1+bb*dconjg(cc1)+cc
  call init_areapart_geom_ref(aref%apref,aa2,bb2,cc2)
else
  gsareapart%apref=>aref
  gsareapart%a_ref=aa
  gsareapart%b_ref=bb
  gsareapart%c_ref=cc
  gsareapart%ntr=gsareapart%apref%ntr
  gsareapart%n=gsareapart%apref%n
  gsareapart%npe=gsareapart%apref%npe
endif
end

subroutine pg_init_areapart_geom_ref_dxdy_rotate(ip,ia,ipart,dx,dy,angle,xr,yr)
!dec$ attributes dllexport:: pg_init_areapart_geom_ref_dxdy_rotate
!инициализировать участок ссылкой на другой участок с поворотом и смещением
use pgmod
integer(4) ip !номер задачи
integer(4) ia !номер области
integer(4) ipart !номер участка
real(8) dx,dy !смещение
real(8) angle  !угол поворота против часовой стрелки (в радианах)
real(8) xr,yr  !точка, вокруг которой осуществляется поворот
real(8) ax,bx,cx,ay,by,cy,sina,cosa
sina=dsin(angle)
cosa=dcos(angle)
ax=cosa
bx=-sina
cx=-xr*cosa+yr*sina+xr+dx
ay=sina
by=cosa
cy=-xr*sina-yr*cosa+yr+dy
call pg_init_areapart_geom_ref(ip,ia,ipart,ax,bx,cx,ay,by,cy)
end

subroutine pg_init_subdmain_connectivity(test_point_coordinates)
!dec$ attributes dllexport:: pg_init_subdmain_connectivity
!определение внутренних границ между сабдоменами
!для boundLineType=1
use pgmod
logical test_point_coordinates !проверить помимо geomdetale совпадение всех точек
integer(4) i1,i2,j1,j2,k1,k2,igd1,igd2,ip1,ip2
type(TArea), pointer :: a1,a2
type(TBound), pointer :: b1,b2
type(TBoundline), pointer :: bl1,bl2
type(TBoundline_geomdetale), pointer :: gd1,gd2
real(8), parameter :: eps=1.0d-8
logical identical_boundline,eq_eps,test_gd1
do i1=1,gs%na
  a1=>gs%a(i1)
  do j1=1,a1%nb
    b1=>a1%bnd(j1)
    if (b1%boundLineType/=1) cycle
    do k1=1,b1%nline
      bl1=>b1%line(k1)
      if (bl1%ci%is_internal) cycle
      do i2=i1+1,gs%na
        a2=>gs%a(i2)
        do j2=1,a2%nb
          b2=>a2%bnd(j2)
          if (b2%boundLineType/=1) cycle
          do k2=1,b2%nline
            bl2=>b2%line(k2)
            if (bl2%ci%is_internal) cycle
            if (bl1%igd_end-bl1%igd_begin/=bl2%igd_end-bl2%igd_begin) cycle
            identical_boundline=.false.
            igd2=bl2%igd_end
            do igd1=bl1%igd_begin,bl1%igd_end
              gd1=>b1%geom_detale0(igd1)
              gd2=>b2%geom_detale0(igd2)
              !проверка тип линий
              identical_boundline=gd1%mode==gd2%mode 
              if (.not.identical_boundline) exit
              !DEC$ IF DEFINED (DEBUG)
              if (gd1%mode==0) call gs_print_stop('Error pg_init_subdmain_connectivity. gd%mode=0!!!') !сюда не должны попадать
              !DEC$ ENDIF
              !проверка количества панелей
              identical_boundline=gd1%i_end-gd1%i_begin==gd2%i_end-gd2%i_begin
              if (.not.identical_boundline) exit
              !проверка совпадения координат начала первого с концом второго
              identical_boundline=eq_eps(b1%x(gd1%i_begin),b2%x(gd2%i_end),eps).and.eq_eps(b1%y(gd1%i_begin),b2%y(gd2%i_end),eps)
              if (.not.identical_boundline) exit
              !проверка совпадения координат начала второго с концом первого
              identical_boundline=eq_eps(b1%x(gd1%i_end),b2%x(gd2%i_begin),eps).and.eq_eps(b1%y(gd1%i_end),b2%y(gd2%i_begin),eps)
              if (.not.identical_boundline) exit
              test_gd1=.false.
              select case (gd1%mode)
              case (1)
                test_gd1=.true.
              case (3)
                identical_boundline=eq_eps(gd1%x,gd2%x,eps).and.eq_eps(gd1%y,gd2%y,eps)
              end select
              if (.not.identical_boundline) exit
              !проверяем все координаты точек
              if (test_point_coordinates.or.test_gd1) then
                ip2=gd2%i_end-1
                do ip1=gd1%i_begin+1,gd1%i_end-1
                  identical_boundline=eq_eps(b1%x(ip1),b2%x(ip2),eps).and.eq_eps(b1%y(ip1),b2%y(ip2),eps)
                  if (.not.identical_boundline) exit
                  ip2=ip2-1
                enddo
              endif
              if (.not.identical_boundline) exit
              igd2=igd2-1
            enddo
            if (identical_boundline) then
              bl1%ci%is_internal=.true.
              bl1%ci%ia=i2
              bl1%ci%ibnd=j2
              bl1%ci%ibndl=k2
              bl2%ci%is_internal=.true.
              bl2%ci%ia=i1
              bl2%ci%ibnd=j1
              bl2%ci%ibndl=k1
            endif
          enddo
        enddo
      enddo
    enddo
  enddo
enddo
end

subroutine pg_init_subdmain_block_order(get_first_bound)
!dec$ attributes dllexport:: pg_init_subdmain_graph
!инициализировать порядок блоков в блочно-диагональной матрице
use pgmod
logical get_first_bound   !функция (без аргументов) возвращающая .true. для внешних границ первого ряда для текущего gsbndl
external get_first_bound
logical is_first_bound
integer(4) i,j,k,n,n_old,i0
type(TArea), pointer :: a,a2
type(TBound), pointer :: b
type(TBoundline), pointer :: bl
!первый проход - определение первых областей с использованием внешней функции get_first_bound
n=gs%na
do i=1,gs%na
  call pg_bind_domain(i)
  do j=1,gsarea%nb
    call pg_bind_bound(j)
    do k=1,gsbnd%nline
      call pg_bind_boundline(k)
      is_first_bound=get_first_bound()
      if (is_first_bound) then
        gsarea%iorder=1
        n=n-1
        exit
      endif
    enddo
    if (is_first_bound) exit
  enddo
enddo
if (n==gs%na) call gs_print_stop("Error pg_init_subdmain_block_order! (no first bound)")
i0=1
do while(n>0)
  n_old=n
  do i=1,gs%na
    a=>gs%a(i)
    if (a%iorder/=i0) cycle
    do j=1,a%nb
      b=>a%bnd(j)
      if (b%boundLineType/=1) cycle
      do k=1,b%nline
        bl=>b%line(k)
        if (.not.bl%ci%is_internal) cycle
        a2=>gs%a(bl%ci%ia)
        if (a2%iorder==0) then
          a2%iorder=i0+1
          n=n-1
        endif
      enddo
    enddo
  enddo
  i0=i0+1
  if (n_old==n) exit
enddo
if (n>0) call gs_print_stop("Error pg_init_subdmain_block_order! (n>0)")
end

subroutine pg_rebuild_boundline(ibndl,n0)
!dec$ attributes dllexport:: pg_rebuild_boundline_ds
!перестроить участок границы с равномерным шагом ds
!!!только для участков границы, которые имеют один geom_detale (не использовался add при создании)
use pgmod
integer(4) ibndl !номер участка границы
integer(4) n0    !желаемое количество панелей
call rebuild_boundline(ibndl,d0,n0,2)
end

subroutine pg_rebuild_boundline_ds(ibndl,ds0)
!dec$ attributes dllexport:: pg_rebuild_boundline_ds
!перестроить участок границы с равномерным шагом ds
!!!только для участков границы, которые имеют один geom_detale (не использовался add при создании)
use pgmod
integer(4) ibndl !номер участка границы
real(8) ds0       !желаемый шаг
call rebuild_boundline(ibndl,ds0,0,1)
end

subroutine rebuild_boundline(ibndl,ds0,n0,mode)
!перестроить участок границы с равномерным шагом ds
!!!только для участков границы, которые имеют один geom_detale (не использовался add при создании)
use pgmod
integer(4) ibndl !номер участка границы
real(8) ds0      !желаемый шаг
integer(4) n0    !желаемое количество панелей
integer(4) mode  !1 - по ds
                 !2 - по n
real(8), allocatable :: s(:),x(:),y(:),bb(:),cc(:,:)
type(TBoundline), pointer :: bl
type(TBoundline_geomdetale), pointer :: gd
integer(4) i,n
real(8) t,l,ds
bl=>gsbnd%line(ibndl)
if (bl%igd_begin/=bl%igd_end) call gs_print_stop("Error pg_rebuild_boundline_ds! bl%igd_begin/=bl%igd_end")
allocate(s(bl%npanel+1))
s(1)=d0
do i=1,bl%npanel
  ds=dsqrt((bl%x(i+1)-bl%x(i))**2+(bl%y(i+1)-bl%y(i))**2)
  s(i+1)=s(i)+ds
enddo
allocate(bb(bl%npanel+1),cc(4,bl%npanel+1))
l=s(bl%npanel+1)
n=n0
if (mode==1) n=idnint(l/ds0)
allocate(x(n+1),y(n+1))
call dcsakm(bl%npanel+1,s,bl%x,bb,cc)
x(1)=bl%x(1)
do i=2,n
  t=(i-d1)/n*l
  x(i)=dcsval(t,n,bb,cc)
enddo
x(n+1)=bl%x(bl%npanel+1)
call dcsakm(bl%npanel+1,s,bl%y,bb,cc)
y(1)=bl%y(1)
do i=2,n
  t=(i-d1)/n*l
  y(i)=dcsval(t,n,bb,cc)
enddo
y(n+1)=bl%y(bl%npanel+1)
call init_boundline_geom(ibndl,n,x,y,.false.)
gd=>gsbnd%geom_detale0(bl%igd_begin)
if (gd%mode==1) then
  gd%i_begin=1
  gd%i_end=n+1
endif
deallocate(s,bb,cc,x,y)
end

subroutine pg_getboundxy_by_s(s,x,y,i)
!dec$ attributes dllexport:: pg_getboundxy_by_s
!получить координаты x,y точки текущей границы по дуговой абсциссе s
use pgmod
real(8) s !дуговая абсцисса
real(8) x,y !выходные координаты x,y
integer(4) i ! номер панели, с которой начнется поиск
real(8) s1,l,sred
l=gsbnd%s(gsbnd%npanel+1)
s1=s-FLOOR(s/l)*l
if (s1<=d0) then
  i=1
elseif (s1>=l) then
  i=gsbnd%npanel
else
  do while (gsbnd%s(i)>s1.or.gsbnd%s(i+1)<s1)
    i=i+1
    if (i>gsbnd%npanel) i=1
  enddo
endif
x=sred(gsbnd%s(i),gsbnd%s(i+1),gsbnd%x(i),gsbnd%x(i+1),s1)
y=sred(gsbnd%s(i),gsbnd%s(i+1),gsbnd%y(i),gsbnd%y(i+1),s1)
end
subroutine init_geom2mode(mode,square_cell,nng,par)
!инициализация геометрии внешнейобласти по mode
use pgmod2
integer(4) nng,mode,nb
logical square_cell
real(8) par(5),r_ellipse,r_circle,r_virus,xb(nmax3),yb(nmax3)
external r_ellipse,r_circle,r_virus
nng=nj_
if (mode==0) then
  !circle
  par(1)=h2
  call init_geom2_23(square_cell,r_circle,par,nng,mode)
elseif (mode==1) then
  !ellipse
  par(1)=d1
  par(2)=d5
  call init_geom2_23(square_cell,r_ellipse,par,nng,mode)
elseif (mode==2) then
  !virus
  nng=151
  par(1)=d1
  par(2)=0.1d0
  par(3)=12.0d0
  call init_geom2_23(square_cell,r_virus,par,nng,mode)
elseif (mode==3) then
  !square
  xb(1)=-h2
  yb(1)=d0
  xb(2)=-h2
  yb(2)=h2
  xb(3)=h2
  yb(3)=h2
  xb(4)=h2
  yb(4)=d0
  nb=4
  nng=151
  call init_rpoly(nb,xb,yb)
  call init_geom2_4(square_cell,xb,yb,nb,nng)
elseif (mode==4) then
  !triangle
  xb(1)=-d1
  yb(1)=d0
  xb(2)=d1
  yb(2)=d1
  xb(3)=d1
  yb(3)=d0
  nb=3
  nng=151
  call init_rpoly(nb,xb,yb)
  call init_geom2_4(square_cell,xb,yb,nb,nng)
endif
end

subroutine init_geom2(square_cell)
!инициализация геометрии ячейки кувабара 
use pgmod2
real(8) par(5),r_circle
external r_circle
logical square_cell
par(1)=d1
call init_geom2_23(square_cell,r_circle,par,nj_,0)
end

subroutine init_geom2_vnesh(square_cell,rmin,rmin2,ds1,ds2)
!инициализация геометрии кроговой или прямоугольной ячейки кувабара с панелями равной длины
!без геометрии внутренней области
!на внутреннем круге панели меньшего размера для уточнения решения
use pgmod2
logical square_cell
real(8) ds1,ds2,rmin,rmin2
integer(4) j,pg_get_int

!прав линия симметрии
call pg_init_boundline_geomlineds2(1,ds1,ds2,rmin,d0,h,d0)
call pg_bind_boundline(1)
njr_=pg_get_int(2,1)
if (square_cell) then
  !внеш граница правая сторона
  call pg_init_boundline_geomlineds(2,ds2,h,d0,h,h)
  !внеш граница верх
  call pg_init_boundline_geomlineds(3,ds2,h,h,-h,h)
  !внеш граница левая сторона
  call pg_init_boundline_geomlineds(4,ds2,-h,h,-h,d0)
  j=5
else
  !внеш круг
  call pg_init_boundline_geomcircleds(2,ds2,d0,d0,h,d0,pi)
  j=3
endif
!лев линия симметрии
call pg_init_boundline_geomlineds2(j,ds2,ds1,-h,d0,-rmin2,d0)
end

subroutine init_geom2_23(square_cell,fr,par,nng,mode)
!инициализация геометрии круговой или прямоугольной ячейки с панелями равной длины
!с обобщением на внутреннюю геометрию произвольной формы
!на внутреннем круге панели меньшего размера для уточнения решения
use pgmod2
logical square_cell
real(8) gg(nmax),dg,ds1,ds2,rmin,rmin2
external fr   !функция внутренней геометрии
real(8) fr,par(5) !параметры функции
integer(4) i,j,nj,nj1,nng,mode
real(8) x(nmax), y(nmax)

if (square_cell) then
  j=6
else
  j=4
endif
call pg_allocate_boundlines(j)
ds1=pi/nng
ds2=0.2d0
if (h<=2.0d0) ds2=0.05d0
if (ds2<ds1) ds2=ds1
!ds1=ds1*2
!ds2=ds2*2
rmin=fr(d0,par)
rmin2=fr(pi,par)
!внешняя граница
call init_geom2_vnesh(square_cell,rmin,rmin2,ds1,ds2)
!внутр граница
if (mode==0) then
  call pg_init_boundline_geomcircleds(j,ds1,d0,d0,par(1),pi,d0)
  nj_=nng
else
  nj=pi/ds1
  nj_=nj
  nj1=nj+1
  dg=pi/nj
  do i=1,nj1
    gg(i)=pi-(i-d1)*dg
  enddo
  do i=1,nj1
    x(i)=fr(gg(i),par)*dcos(gg(i))
    y(i)=fr(gg(i),par)*dsin(gg(i))
  enddo
  call pg_init_boundline_geom(j,nj,x,y)
endif
end

subroutine init_geom3_dy(dy)
!инициализация геометрии для квадратной области со сдвигом по Y
use pgmod2
real(8) dy
nj_=10
call pg_allocate_boundlines(4)
call pg_init_boundline_geomline(1,nj_,d0,dy,d1,dy)
call pg_init_boundline_geomline(2,nj_,d1,dy,d1,d1+dy)
call pg_init_boundline_geomline(3,nj_,d1,d1+dy,d0,d1+dy)
call pg_init_boundline_geomline(4,nj_,d0,d1+dy,d0,dy)
end

subroutine init_geom3
!инициализация геометрии для квадратной области
use gen_mod
call init_geom3_dy(d0)
end

subroutine init_geom3_2
!инициализация геометрии для половины квадратной области
use pgmod2
real(8) ds,dx,dy
real(8), parameter :: eps=1.0d-8
dx=d1/nsx
dy=d1/nsy
ds=dx
if (dy<dx) ds=dy
nj_=10
ds=ds/nj_
nj_=(dx+eps)/ds
nj_2=(dy+eps)/ds
call pg_allocate_boundlines(4)
call pg_init_boundline_geomlineds(1,ds,d0,d0,dx,d0)
call pg_init_boundline_geomlineds(2,ds,dx,d0,dx,dy)
call pg_init_boundline_geomlineds(3,ds,dx,dy,d0,dy)
call pg_init_boundline_geomlineds(4,ds,d0,dy,d0,d0)
end

subroutine init_geom4
!инициализация геометрии внутренней области 
use pgmod2
type(TBound), pointer :: b
type(TBoundline), pointer :: bl
call pg_allocate_boundlines(2)
!полукруг
b=>gs%a(1)%bnd(1)
call pg_copy_boundline(gs%i,1,1,1,b%nline)
!ось симметрии
bl=>gsbnd%line(1)
call pg_init_boundline_geomlineds(2,pi/bl%npanel,bl%x(bl%npanel+1),d0,bl%x(1),d0)
end

subroutine get_rr(rr,dg,jrmax)
!инициализация массива rr с возрастанием, чтобы ячейки были близки к квадратным
use pgmod2
integer(4) jrmax
real(8) dg,rr(nmax)
call ga_get_rr(h2,hv,rr,dg,jrmax,nmax,5,0)
end

function r_circle(t,par)
use gen_mod
real(8) r_circle,t,par(5),t1
t1=t !чтобы не было предупреждения о неиспользуемой переменной
r_circle=par(1)
end

function r_ellipse(t,par)
use gen_mod
real(8) r_ellipse,t,par(5),a2,b2,t2
a2=par(1)**2
b2=par(2)**2
t2=dtan(t)**2
r_ellipse=dsqrt(a2*b2*(d1+t2)/(b2+a2*t2))
end

function r_virus(t,par)
use gen_mod
real(8) r_virus,t,par(5),r,g,mb
r=par(1)
g=par(2)
mb=par(3)
r_virus=r*(1+g*dcos(mb*t))
end

function r_poly(t,par)
!радиус-вектор для границы полилинии
use pgmod2
real(8) r_poly,t,sin_fi_tt,s,par(5),t1
integer(4) i,j
t1=par(1)  !чтобы не было предупреждения о неиспользуемой переменной
do i=2,rpoly_n
  j=i-1
  if (t<=rpoly_g(i)) exit
enddo
sin_fi_tt=dsin(t-rpoly_tt(j))
s=(rpoly_y(j)*dcos(t)-rpoly_x(j)*dsin(t))/sin_fi_tt
r_poly=dsqrt(rpoly_x(j)*rpoly_x(j)+rpoly_y(j)*rpoly_y(j)+s*s+d2*s*(rpoly_x(j)*rpoly_ksi(j)+rpoly_y(j)*rpoly_eta(j)))
end

subroutine init_rpoly(nb,xb,yb)
!инициализация массивов для функции r_poly
use pgmod2
real(8) xb(nmax3),yb(nmax3),s
integer(4) nb, i,j,k
rpoly_n=nb
j=nb
do i=1,nb
  rpoly_x(i)=xb(j)
  rpoly_y(i)=yb(j)
  rpoly_g(i)=datan2(yb(j),xb(j))
  if (i>1) then
    k=i-1
    rpoly_ksi(k)=rpoly_x(i)-rpoly_x(k)
	rpoly_eta(k)=rpoly_y(i)-rpoly_y(k)
	s=dsqrt(rpoly_ksi(k)*rpoly_ksi(k)+rpoly_eta(k)*rpoly_eta(k))
	rpoly_ksi(k)=rpoly_ksi(k)/s
	rpoly_eta(k)=rpoly_eta(k)/s
	rpoly_tt(k)=datan2(rpoly_eta(k),rpoly_ksi(k))
  endif
  j=j-1
enddo
rpoly_g(1)=d0
rpoly_g(nb)=pi
end

subroutine init_vnesh(nng,bndg,bndrv)
!инициализация внешней границы для прямоугольной области
use pgmod2
integer(4) nng
real(8) bndg(nmax),bndrv(nmax)
call ga_init_vneshg(nng,bndg,bndrv,h,h,3)
end

subroutine init_geom2_4(square_cell,xb,yb,nb,nng)
!инициализация геометрии круговой или прямоугольной ячейки с панелями равной длины
!с обобщением на внутреннюю геометрию полигональной формы
!на внутреннем круге панели меньшего размера для уточнения решения
use pgmod2
logical square_cell
real(8) ds1,ds2,dx,dy,rmin,rmin2,xb(nmax3),yb(nmax3)
real(8) ss,ds
integer(4) i,j,nng,nb,nx,k,nl
real(8) x(nmax), y(nmax)
if (square_cell) then
  nl=6
else
  nl=4
endif
call pg_allocate_boundlines(nl)
ds1=h2*pi/nng
ds2=d2*ds1 !0.2d0
!ds1=ds1*2
!ds2=ds2*2
rmin=xb(nb)
rmin2=-xb(1)
call init_geom2_vnesh(square_cell,rmin,rmin2,ds1,ds2)
!внутр граница
j=1
x(j)=xb(1)
y(j)=yb(1)
do i=2,nb
  dx=xb(i)-xb(i-1)
  dy=yb(i)-yb(i-1)
  ss=dsqrt(dx*dx+dy*dy)
  nx=ss/ds1
  ds=ss/nx
  if (nx<1) nx=1
  do k=1,nx
    j=j+1
	ss=(k-d0)/nx
    x(j)=xb(i-1)+ss*dx
	y(j)=yb(i-1)+ss*dy
  enddo
enddo
j=j-1
nj_=j
call pg_init_boundline_geom(nl,j,x,y)
end

subroutine init_geom4_1                           
!инициализация геометрии полукругом
use pgmod2
real(8) ds
call pg_allocate_boundlines(2)
ds=pi/nj_
!круг
call pg_init_boundline_geomcircleds(1,ds,d0,d0,d1,d0,pi)
!линия
call pg_init_boundline_geomlineds(2,ds,-d1,d0,d1,d0)
end

subroutine init_geom4_2                           
!инициализация геометрии полукругом
!с тремя участками границы
use pgmod2
real(8) ds
call pg_allocate_boundlines(3)
ds=pi/nj_
!линия
call pg_init_boundline_geomlineds(1,ds,-d0,d0,d1,d0)
!круг
call pg_init_boundline_geomcircleds(2,ds,d0,d0,d1,d0,pi)
!линия
call pg_init_boundline_geomlineds(3,ds,-d1,d0,d0,d0)
end

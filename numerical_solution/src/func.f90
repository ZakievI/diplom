function sred(a,b,fa,fb,c)
! среднее
real(8) sred,a,b,fa,fb,c
sred=fa+(fb-fa)/(b-a)*(c-a)
end

subroutine sred_aray(a,b,m,n)
!равномерное разбиение от а до b
integer(4) n
real(8) a,b,m(n),d
d=(b-a)/(n-1)
m=(/(a+d*i,i=0,n-1)/)
end

function csred(a,b,fa,fb,c)
! среднее
real(8) a,b,c
complex(8) csred,fa,fb
csred=fa+(fb-fa)/(b-a)*(c-a)
end

function scal_p(a,b)
!скалярное произведение
real(8) scal_p
complex(8) a,b
scal_p=dreal(dconjg(a)*b)
end

function scal_p2(ax,ay,bx,by)
!скалярное произведение
real(8) scal_p2
real(8) ax,ay,bx,by
scal_p2=ax*bx+ay*by
end

function vect_p(a,b)
!векторное произведение
real(8) vect_p
complex(8) a,b
vect_p=dimag(dconjg(a)*b)
end

function vect_p2(ax,ay,bx,by)
!векторное произведение
real(8) vect_p2
real(8) ax,ay,bx,by
vect_p2=ax*by-ay*bx
end

subroutine matr_vect_p(mm,v,r,n,m,lda)
!произведение матрицы mm(n,m) на вектор v(n)
!результат r(m)
use gen_mod
integer(4) n,m !число строк и столбцов
integer(4) lda !максимальное число элеменов по первому индексу
real(8) mm(n,m),v(m),r(n)
call dgemv('N', n, m, d1, mm, lda, v, 1, d0, r, 1)
end

function collenear3z(z1,z2,z3)
use gen_mod
complex(8) z1,z2,z3,a,b
real(8) vect_p,x
logical collenear3z
a=z2-z1
b=z3-z1
x=vect_p(a,b)/cdabs(a*b)
collenear3z=dabs(x)<1.0d-8
end

function get_s_section(z,z1,z2)
!вычисление дуговой абсциссы проекции точки z на отрезок (z1,z2)
!s=0 - проекция совпадает с z1
!s=1 - проекция совпадает с z2
complex(8) z,z1,z2,dz21,dz1
real(8) get_s_section
dz21=z2-z1
dz1=z-z1
get_s_section=dreal(dconjg(dz21)*dz1)/cdabs(dz21)**2
end

subroutine sort_mass(n,a)
!сотрировка массива
integer(4) i,j,n,jmin
real(8) a(n),tmin
do i=1,n-1
  tmin=a(i)
  jmin=i
  do j=i+1,n
    if (a(j)<tmin) then
	  tmin=a(j)
	  jmin=j
	endif
  enddo
  a(jmin)=a(i)
  a(i)=tmin
enddo
end

subroutine angles_razryv(nj1,alpha)
!удаление разрыва при скачках функций на 2\pi
use gen_mod
integer(4) i,nj1
real(8) dd,alpha(nj1),dpi
  dpi=pi !*1.5d0
  dd=d0
  do i=2,nj1
    alpha(i)=alpha(i)+dd
    if ((alpha(i)-alpha(i-1))>dpi) then 
      dd=dd-pi2
	  alpha(i)=alpha(i)-pi2
    endif
    if ((alpha(i)-alpha(i-1))<-dpi) then 
      dd=dd+pi2
	  alpha(i)=alpha(i)+pi2
    endif
  enddo
end

subroutine spline_razryv(xx,yy,bb,cc,nr,nmax,razryv)
!построение сплайна от функции с разрывной производной
use gen_mod
integer(4) i,j,k,nr,nmax,razryv(nmax)
real(8) xx(nmax),yy(nmax),bb(nmax),cc(4,nmax)
do i=1,nr-1
  j=razryv(i)
  k=razryv(i+1)
  call dcsakm(k-j+1,xx(j:k),yy(j:k),bb(j:k),cc(:,j:k))
enddo
end

subroutine sort_for_spline(nj1,xx,yy,bb,yy1)
!сортировка точек для сплайна
use gen_mod
integer(4) i,nj1,j,jmin
real(8) xx(nj1),yy(nj1),bb(nj1),yy1(nj1)
real(8) tmin,fmin
bb=xx
yy1=yy
do i=1,nj1-1
  tmin=bb(i)
  fmin=yy1(i)
  jmin=i
  do j=i+1,nj1
    if (bb(j)<tmin) then
	  tmin=bb(j)
	  fmin=yy1(j)
	  jmin=j
	endif
  enddo
  bb(jmin)=bb(i)
  bb(i)=tmin
  yy1(jmin)=yy1(i)
  yy1(i)=fmin
enddo
end

subroutine spline_linear(nj1,xx,yy,bb,cc)
!линейный сплайн
use gen_mod
integer(4) i,nj1
real(8) xx(nj1),yy(nj1),bb(nj1),cc(4,nj1),yy1(nj1)
call sort_for_spline(nj1,xx,yy,bb,yy1)
cc=d0
do i=1,nj1-1
  cc(1,i)=yy1(i) 
  cc(2,i)=(yy1(i+1)-yy1(i))/(bb(i+1)-bb(i))
enddo
cc(1,nj1)=yy1(nj1)
cc(2,nj1)=cc(2,nj1-1)
end


subroutine spline_const(nj1,xx,yy,bb,cc)
!линейный сплайн
use gen_mod
integer(4) nj1
real(8) xx(nj1),yy(nj1),bb(nj1),cc(4,nj1),yy1(nj1)
call sort_for_spline(nj1,xx,yy,bb,yy1)
cc=d0
cc(1,:)=yy1(:) 
end

subroutine findline(nt,tj,tt,x0,y0,lx,ly)
!определение формы линии по углу наклона касательной
use gen_mod
integer(4) nt       !число точек
real(8) tj(nt)      !дуговая абсцисса 
real(8) tt(nt)      !угол наклона касательной к кривой
real(8) x0,y0       !начальная точка
real(8) lx(nt),ly(nt) !координаты точек кривой
complex(8) f(nt)

f=c1
call findlinez(nt,tj,tt,f,x0,y0,lx,ly)
end

subroutine findlinez(nt,tj,tt,f,x0,y0,lx,ly)
!определение формы линии интегрированием уравнения dz = f(\zeta)d\zeta
use gen_mod
integer(4) nt       !число точек
real(8) tj(nt)      !дуговая абсцисса 
real(8) tt(nt)      !угол наклона касательной к кривой
complex(8) f(nt)    !функция f(\zeta)
real(8) x0,y0       !начальная точка
real(8) lx(nt),ly(nt) !координаты точек кривой
integer(4) i
complex(8) fc(nt)

do i=1,nt
  fc(i)=f(i)*cdexp(ii*tt(i))
enddo
call findlinez1(nt,tj,fc,x0,y0,lx,ly)
end

subroutine findlinez1(nt,tj,f,x0,y0,lx,ly)
!определение формы линии интегрированием уравнения dz = f(\zeta)d\sigma
!d\sigma=|d\zeta|
use gen_mod
integer(4) nt       !число точек
real(8) tj(nt)      !дуговая абсцисса 
complex(8) f(nt)    !функция f(\zeta)
real(8) x0,y0       !начальная точка
real(8) lx(nt),ly(nt) !координаты точек кривой
integer(4) i
real(8) ff(nt)
real(8) bb(nt),cc(4,nt)

ff=dreal(f)
call dcsakm(nt,tj,ff,bb,cc)
lx(1)=x0
do i=2,nt
  lx(i)=lx(i-1)+dcsitg(tj(i-1),tj(i),nt-1,bb,cc)
enddo
ff=dimag(f)
call dcsakm(nt,tj,ff,bb,cc)
ly(1)=y0
do i=2,nt
  ly(i)=ly(i-1)+dcsitg(tj(i-1),tj(i),nt-1,bb,cc)
enddo
end

function cint(nt,tj,f,t0,t1)
!вычисление комплексного интеграла \int_t0^t1 f(t) dt
use gen_mod
complex(8) cint
integer(4) nt       !число точек
real(8) tj(nt)      !дуговая абсцисса 
complex(8) f(nt)    !функция f(t)
real(8) t0,t1       !пределы интегрирования 
real(8) d_int
cint=c1*d_int(nt,tj,dreal(f),t0,t1)+ii*d_int(nt,tj,dimag(f),t0,t1)
end

function d_int(nt,tj,f,t0,t1)
!вычисление вещественного интеграла \int_t0^t1 f(t) dt
use gen_mod
real(8) d_int
integer(4) nt       !число точек
real(8) tj(nt)      !дуговая абсцисса 
real(8) f(nt)    !функция f(t)
real(8) t0,t1       !пределы интегрирования 
real(8) bb(nt),cc(4,nt)
call dcsakm(nt,tj,f,bb,cc)
d_int=dcsitg(t0,t1,nt-1,bb,cc)
end


function cintz(nt,t,f)
!вычисление комплексного интеграла \int_t0^t1 f(t) dt
use gen_mod
complex(8) cintz
integer(4) nt       !число точек
complex(8) t(nt)    !координаты линии интегрирования
complex(8) f(nt)    !функция f(t)
complex(8) cint
real(8) tj(nt)      !дуговая абсцисса 
complex(8) f1(nt)
real(8) xx(nt),yy(nt),bb(nt),cc(4,nt)
integer(4) i
tj(1)=d0
do i=2,nt
  tj(i)=tj(i-1)+cdabs(t(i)-t(i-1))
enddo
do i=1,nt
  xx(i)=dreal(t(i))
  yy(i)=dimag(t(i))
enddo
call dcsakm(nt,tj,xx,bb,cc)
do i=1,nt
  xx(i)=dcsder(1,tj(i),nt-1,bb,cc)
enddo
call dcsakm(nt,tj,yy,bb,cc)
do i=1,nt
  yy(i)=dcsder(1,tj(i),nt-1,bb,cc)
enddo
do i=1,nt
  f1(i)=f(i)*(c1*xx(i)+ii*yy(i))
enddo
cintz=cint(nt,tj,f1,tj(1),tj(nt))
end

function DistanceToPoint(x1,y1,x2,y2)
use gen_mod
real(8) x1,y1,x2,y2,pcx,pcy,DistanceToPoint
pcx=x1-x2
pcy=y1-y2
DistanceToPoint=dsqrt(pcx*pcx+pcy*pcy);
end


function DistanceToSegment(ax,ay,bx,by,cx,cy,s,px,py)
!расстояние до отрезка
!s - нормированная дуговая абсцисса проекции на отрезок =0 для первой точки =1 для второй точки (<0 или >1) - если точка вне отрезка
use gen_mod
real(8) ax,ay,bx,by,cx,cy,s,px,py,acx,acy,abx,aby,abl2,DistanceToSegment,DistanceToPoint
acx=cx-ax
acy=cy-ay
abx=bx-ax
aby=by-ay
abl2=abx*abx+aby*aby
if (abl2<1.0d-10) then
  s=0
  DistanceToSegment=dsqrt(acx*acx+acy*acy)
  return
endif
s=(acx*abx+acy*aby)/abl2
if (s<=d0) then
  px=ax
  py=ay
  DistanceToSegment=dsqrt(acx*acx+acy*acy)
  s=d0
elseif (s>=d1) then
  px=bx
  py=by
  DistanceToSegment=DistanceToPoint(cx,cy,bx,by)
  s=d1
else
  px=ax+s*abx
  py=ay+s*aby
  DistanceToSegment=DistanceToPoint(cx,cy,px,py)
endif
end

function DistanceToLine2(xx,yy,x,y,npanel,j,smin,px_min,py_min,first_i)
!расстояние до линии
use gen_mod
integer(4) npanel               !число отрезков
real(8) x(npanel+1),y(npanel+1) !координаты узлов
real(8) xx,yy                   !координаты точки
integer(4) j                    !индекс отрезка с ближайшей точкой
real(8) smin                    !локальная дуговая абсцисса на j-м отрезке (0..1)
real(8) px_min,py_min           !ближайшая точка на линии
integer(4) first_i              !индекс первого сегмента, с которого нужно начинать поиск (по умолчанию = 1)
real(8) DistanceToLine2          !расстояние до линии
integer(4) i
real(8) dist,minDist,s,px,py,DistanceToSegment
minDist=1.0d10
do i=first_i,npanel
  dist=DistanceToSegment(x(i),y(i),x(i+1),y(i+1),xx,yy,s,px,py)
  if (dist<minDist) then
    minDist=dist
    j=i
    smin=s
    px_min=px
    py_min=py
  endif
enddo
if (smin<d0) then
  smin=d0
elseif (smin>d1) then
  smin=d1
endif
DistanceToLine2=minDist
end

function DistanceToLine(xx,yy,x,y,npanel,j,smin)
!расстояние до линии
use gen_mod
integer(4) npanel               !число отрезков
real(8) x(npanel+1),y(npanel+1) !координаты узлов
real(8) xx,yy                   !координаты точки
integer(4) j                    !индекс отрезка с ближайшей точкой
real(8) smin                    !локальная дуговая абсцисса на j-м отрезке (0..1)
real(8) DistanceToLine          !расстояние до линии
real(8) px,py,DistanceToLine2
DistanceToLine=DistanceToLine2(xx,yy,x,y,npanel,j,smin,px,py,1)
end

function DistanceToCircleArc(xx,yy,x0,y0,r,gam1,gam2,gam)
!расстояние до дуги окружности
use gen_mod
real(8) xx,yy !координаты точки
real(8) x0,y0 !координаты центра круга
real(8) r !радиус
real(8) gam1,gam2 !угловые координаты начала и конца дуги
real(8) gam !углаовая коордианата ближайшей точки на дуге
real(8) DistanceToCircleArc,DistanceToPoint
real(8) rr,r1,r2
integer(4) i
call rr_tt(xx-x0,yy-y0,rr,gam)
if (gam2>gam1) then
  i=int((gam-gam1)/pi2)
  if (gam-gam1<d0) i=i-1
else
  i=int((gam-gam2)/pi2)
  if (gam-gam2<d0) i=i-1
endif
if (i.ne.0) gam=gam-i*pi2
if ((gam-gam1)*(gam-gam2)<=d0) then
  DistanceToCircleArc=dabs(r-rr)
else
  r1=DistanceToPoint(xx,yy,x0+r*dcos(gam1),x0+r*dsin(gam1))
  r2=DistanceToPoint(xx,yy,x0+r*dcos(gam2),x0+r*dsin(gam2))
  if (r1<r2) then
    DistanceToCircleArc=r1
    gam=gam1
  else
    DistanceToCircleArc=r2
    gam=gam2
  endif
endif
end

subroutine rr_tt(x,y,rr,tt)
use gen_mod
real(8) rr,tt
real(8) x,y
rr=dsqrt(x**2+y**2)
tt=datan2(y,x)
if (tt<d0) tt=tt+pi2
end

subroutine swap_array(a,n)
use gen_mod
integer(4) n,i,j
real(8) a(n),t
do i=1,n/2
  t=a(i)
  j=n-i+1
  a(i)=a(j)
  a(j)=t
enddo
end

subroutine swap_int(i1,i2)
integer(4) i1,i2,j
j=i1
i1=i2
i2=j
end

subroutine swap_double(i1,i2)
real(8) i1,i2,j
j=i1
i1=i2
i2=j
end

function ceq_eps(x,y,eps)
logical ceq_eps
complex(8) x,y
real(8) eps
ceq_eps=cdabs(x-y)<eps
end

function eq_eps(x,y,eps)
logical eq_eps
real(8) eps,x,y
eq_eps=dabs(x-y)<eps
end

subroutine line_array_ds(x1,y1,x2,y2,ds,x,y,s,n,min_np)
!разбиение отрезка с заданным шагом ds
use gen_mod
integer(4) n,j
real(8) x1,y1,x2,y2 !координаты концов отрезка
real(8) ds !шаг интегрирования
real(8), allocatable :: x(:),y(:),s(:)
integer(4) min_np !минимальное число точек
real(8) dx,dy,ll,ds1
dx=x2-x1
dy=y2-y1
ll=dsqrt(dx**2+dy**2)
n=ll/ds
if (n<min_np) n=min_np
dx=dx/n
dy=dy/n
ds1=ll/n
n=n+1
allocate(x(n),y(n),s(n))
x(1:n)=(/(x1+j*dx,j=0,n-1)/)
y(1:n)=(/(y1+j*dy,j=0,n-1)/)
s(1:n)=(/(j*ds1,j=0,n-1)/)
!!!освобождение массива снаружи
end

function test_array_ne_0(x,n)
!проверяем, что в массиве есть ненулевые элементы
use gen_mod
integer(4) n,k
real(8) x(n)
logical is_x_ne_0,test_array_ne_0
is_x_ne_0=.false.
do k=1,n
  is_x_ne_0=x(k)/=d0
  if (is_x_ne_0) exit
enddo
test_array_ne_0=is_x_ne_0
end

function itoa(i) result(res)
!integer to string
  character(:),allocatable :: res
  integer(4),intent(in) :: i
  character(range(i)+2) :: tmp
  write(tmp,'(i0)') i
  res = trim(tmp)
end

function ftoa(f) result(res)
!real(8) to string
  character(:),allocatable :: res
  real(8),intent(in) :: f
  character(30) :: tmp
  write(tmp,*) f
  res = trim(tmp)
end

subroutine LS_approx(n,m,x,y,cc)
!нахождение коэффициентов аппроксимационного полинома f(x)=\sum_i=0^n cc(i)*x^i
!методом наименьших квадратов
use gen_mod
integer(4) n !степень полинома
integer(4) m !число точек (x_j,y_j)
real(8) x(m),y(m) !координаты аппроксимируемых точек
real(8) cc(0:n)   !коэффициенты аппроксимационного полинома
real(8) mm(0:n,0:n),b(0:n)
integer(4) i,j,k
mm=d0
b=d0
do k=0,n
  do i=0,n
    do j=1,m
      mm(k,i)=mm(k,i)+x(j)**i*x(j)**k
    enddo
  enddo
  do j=1,m
    b(k)=b(k)+y(j)*x(j)**k
  enddo
enddo
call DLSLRG(n+1,mm,n+1,b,1,cc)
end

subroutine LS_approx2(n,m,x,y,w,cc,fcn_LS_approx)
!нахождение коэффициентов аппроксимационной функции F(x)=\sum_i=1^n cc(i)*f_i(x)
!методом наименьших квадратов
use gen_mod
integer(4) n !количество базисных функций
integer(4) m !число точек (x_j,y_j)
real(8) x(m),y(m) !координаты аппроксимируемых точек
real(8) w(m)      !вес каждой точки L=sum_(j=1)^m w(j)*(F(x_j)-y_j)^2
real(8) cc(n)   !коэффициенты аппроксимационной функции
real(8) fcn_LS_approx
external fcn_LS_approx !базисные функции fcn_LS_approx(i,x), i - номер функции, x - координата x 
real(8) mm(n,n),b(n),fk(m)
integer(4) i,j,k
mm=d0
b=d0
do k=1,n
  do j=1,m
    fk(j)=fcn_LS_approx(k,x(j))*w(j)
    b(k)=b(k)+y(j)*fk(j)
  enddo
  do i=1,n
    do j=1,m
      mm(k,i)=mm(k,i)+fcn_LS_approx(i,x(j))*fk(j)
    enddo
  enddo
enddo
call DLSLRG(n,mm,n,b,1,cc,fcn_LS_approx)
end

subroutine LS_approx3(n,m,x,y,cc,fcn_LS_approx)
!нахождение коэффициентов аппроксимационной функции общего вида F(cc,x)
!методом наименьших квадратов
use gen_mod
use func_mod
integer(4) n !количество искомых коэффициентов
integer(4) m !число точек (x_j,y_j)
real(8), target :: x(m),y(m) !координаты аппроксимируемых точек
real(8) cc(n)   !коэффициенты аппроксимационной функции
real(8) fcn_LS_approx
external fcn_LS_approx !аппроксимирующая функция fcn_LS_approx(n,cc,x) x - координата одной точки x 
real(8) w(m)
w=d1
call LS_approx4(n,m,x,y,w,cc,fcn_LS_approx)
end

SUBROUTINE func_FCN_LS_approx3(M, N, X, F)
use func_mod
INTEGER(4) N,M,j
REAL(8) X(N), F(M)
do j=1,M
  F(j)=func_fcn_LS_approx(n,x,func_fcn_LS_approx_x(j))-func_fcn_LS_approx_y(j)
  F(j)=F(j)*func_fcn_LS_approx_w(j)
enddo
end

subroutine LS_approx4(n,m,x,y,w,cc,fcn_LS_approx)
!нахождение коэффициентов аппроксимационной функции общего вида F(cc,x)
!методом наименьших квадратов
!с добавлением возможности учета веса каждой точки
use gen_mod
use func_mod
integer(4) n !количество искомых коэффициентов
integer(4) m !число точек (x_j,y_j)
real(8), target :: x(m),y(m) !координаты аппроксимируемых точек
real(8), target :: w(m)      !вес каждой точки L=sum_(j=1)^m w(j)*(F(cc,x_j)-y_j)^2
real(8) cc(n)   !коэффициенты аппроксимационной функции
real(8) fcn_LS_approx
external fcn_LS_approx !аппроксимирующая функция fcn_LS_approx(n,cc,x,j,w) x - координата одной точки x; w - возвращаемый вес
external func_FCN_LS_approx3
real(8) x0(n),FSCALE(m),XSCALE(n),RPARAM(7),FVEC(m)
integer(4) IPARAM(7)
real(8), allocatable :: FJAC(:,:)
func_fcn_LS_approx=>fcn_LS_approx
func_fcn_LS_approx_x=>x
func_fcn_LS_approx_y=>y
func_fcn_LS_approx_w=>w
FSCALE=d1
XSCALE=d1
CALL dU4LSF (IPARAM, RPARAM)
allocate(FJAC(m,n))
x0=cc
CALL dUNLSF (func_FCN_LS_approx3, m, n, x0, XSCALE, FSCALE, IPARAM, RPARAM, cc, FVEC, FJAC, m)
deallocate(FJAC)
end

function GetFourierSeriesCoeff(n,s,f,smin,smax,k,mode) result(res)
!вычисление коэффициентов ряда Фурье f(g)=a0+sum_1^\infty(ak*cos(kg)+bk*sin(kg))
![smin,smax] -> [0,2Pi]
!!!на много точнее, чем через быстрое преобразование Фурье (FFT), но на порядок дольше
use gen_mod
integer(4) n !число точек
real(8) s(n) !разбиение [smin,smax]
real(8) f(n) !функция на разбиении
real(8) smin,smax
integer(4) k !интекс коэффициента
integer(4) mode !0 - ak, 1 - bk
real(8) res,d
real(8), allocatable :: g(:) !разбиение [0,2Pi]
real(8), allocatable :: bb(:),cc(:,:),f1(:)
integer(4) i,n1
n1=n+1
allocate(bb(n1),cc(4,n1),g(n1),f1(n1))
d=pi2/(smax-smin)
forall (i=1:n) g(i)=(s(i)-smin)*d
g(n1)=g(1)+pi2
if (k==0) then
  f1(1:n)=f
else
  if (mode==0) then
    forall (i=1:n) f1(i)=f(i)*dcos(k*g(i))
  else
    forall (i=1:n) f1(i)=f(i)*dsin(k*g(i))
  endif
endif
f1(n1)=f1(1)
call dcsper(n1,g,f1,bb,cc)
res=dcsitg(g(1),g(n1),n,bb,cc)/pi
if (k==0) res=res*d5
deallocate(bb,cc,g,f1)
end

recursive function dbsin(n,x) result(res)
!функция Бесселя I_n(x)
use gen_mod
integer(4) n
real(8) x,res
if (n<0) then
  write(*,*) "Error dbsin. n<0!!!"
  stop
endif
select case (n)
case (0)
  res=dbsi0(x)
case (1)
  res=dbsi1(x)
case default
  res=dbsin(n-2,x)-2*(n-1)*dbsin(n-1,x)/x
endselect
end

recursive function dbskn(n,x) result(res)
!функция Бесселя K_n(x)
use gen_mod
integer(4) n
real(8) x,res
if (n<0) then
  write(*,*) "Error dbskn. n<0!!!"
  stop
endif
select case (n)
case (0)
  res=dbsk0(x)
case (1)
  res=dbsk1(x)
case default
  res=dbskn(n-2,x)+2*(n-1)*dbskn(n-1,x)/x
endselect
end

function real_test_inf_nan(x,x1) result (res)
use IEEE_ARITHMETIC
real(8) x,x1,res
if (ieee_is_nan(x).or.(.not.ieee_is_finite(x))) then
  res=x1
else
  res=x
endif
end

function FindZFromGrid(xx,yy,x0,y0,dx,dy,nx,ny,ff,emptyVal) result(res)
!билейная интерполяция
use gen_mod
real(8) xx,yy !координаты точки
real(8) x0,y0 !левый нижний угол сетки
real(8) dx,dy !шаг сетки
integer(4) nx,ny !число ячеек сетки
real(8) ff(nx+1,ny+1) !массив функции
real(8) emptyVal !пустое значение, в точках, где функция не вычислена
real(8) Ksi,Eta,res,x,y
real(8) x_1,z1,z2,z3,z4,y_1
integer(4) i,j
x=xx
y=yy
x=Max(x,x0)
x=Min(x,x0+nx*dx)
y=Max(y,y0)
y=Min(y,y0+ny*dy)
i=int((x-x0)/dx)
j=int((y-y0)/dy)
i=Max(i,0)
i=Min(i,nx-1)
j=Max(j,0)
j=Min(j,ny-1)
x_1=x0+i*dx
y_1=y0+j*dy
i=i+1
j=j+1
z1=ff(i,j)
z2=ff(i+1,j)
z3=ff(i+1,j+1)
z4=ff(i,j+1)
if (z1.ne.emptyVal.and.z2.ne.emptyVal.and.z3.ne.emptyVal.and.z4.ne.emptyVal) then
  ksi=(x-x_1)/dx*d2-d1
  eta=(y-y_1)/dy*d2-d1
  res=(z1*(d1-ksi)*(d1-eta)+z2*(d1+ksi)*(d1-eta)+z3*(d1+ksi)*(d1+eta)+z4*(d1-ksi)*(d1+eta))*0.25d0
else
  res=emptyVal
endif
end

function split_str(str,i,sep,word) result(res)
!выделение слов из строки, разделенных заданным набором разделителей
character(*) str !строка
character(*) word !(out) выделенное слово
character(*) sep !разделители, может быть несколько символов
integer(4) i  !позиция, начиная с которой ищется слово, на выходе обновляется, начальное значение = 1
logical res,is_sep  !выходное значение =true, если достигнут конец строки
integer(4) j,k,b
res=.false.
word=''
b=-1
do j=i,len(str)
  is_sep=.false.
  do k=1,len(sep)
    is_sep=str(j:j)==sep(k:k)
    if (is_sep) exit
  enddo
  if (is_sep) then
    if (b>0) exit
  else
    if (b<0) b=j
  endif
enddo
if (b>0.and.j>b) word=str(b:j-1)
res=j>len(str)
if (res) then
  i=j
else
  i=j+1
endif
end

function getfpoly(n,cc,x)
!вычисление полинома, заданного массивом коэффициентов в точке x
use gen_mod
integer(4) n !степень полинома
real(8) cc(0:n) !коэффициенты полинома
real(8) x !точка
real(8) s,getfpoly
integer(4) i
s=cc(0)
do i=1,n
  s=s+cc(i)*x**i
enddo
getfpoly=s
end

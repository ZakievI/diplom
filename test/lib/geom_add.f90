subroutine ga_init_vneshg(nng,bndg,bndrv,h,h1,min_ng)
!dec$ attributes dllexport:: ga_init_vneshg
!инициализация внешней границы для прямоугольной области
use gen_mod
integer(4) nng      !число отрезков в разбиении (число точек - 1)
real(8) bndg(nng+1) !разбиение по углу gamma (out)
real(8) bndrv(nng+1)!длина радиус-вектора до внешней границы
real(8) h           !высота прямоугольной области 
real(8) h1          !половина ширины прямоугольлной области
integer(4) min_ng   !минимальное число отрезков на вертикальной границе
integer(4) i,ng,nj2,j
real(8) alfh,dg
alfh=datan(h/h1)
bndg=d0
dg=pi/nng
ng=alfh/dg
if (ng<min_ng) ng=min_ng
do i=1,ng
  bndg(i)=(i-d1)/ng*alfh
  j=nng+1-i+1
  bndg(j)=pi-bndg(i)
  bndrv(i)=h1/dcos(bndg(i))
  bndrv(j)=bndrv(i)
enddo
nj2=nng-2*ng
do i=1,nj2+1
  j=i+ng
  bndg(j)=(i-d1)/nj2*(pi-d2*alfh)+alfh
  bndrv(j)=h/dsin(bndg(j))
enddo
end

subroutine ga_get_rr(rmin,rmax,rr,dg,jrmax,nmax,min_nr,max_nr)
!dec$ attributes dllexport:: ga_get_rr
!инициализация массива rr с возрастанием, чтобы ячейки были близки к квадратным
use gen_mod
real(8) rmin !мин r
real(8) rmax !макс r
integer(4) nmax  !размер буфера rr
real(8) rr(nmax) !массив разбиений от rmin до rmax
real(8) dg !шаг, с корого начинать разбиение (шаг на круге)
integer(4) jrmax !получившееся число точек в rr
integer(4) min_nr !минимальное число точек в rr (если = 0, то не учитывается)
integer(4) max_nr !максимальное число точек в rr (если = 0, то не учитывается)
integer(4) j
real(8) k
j=1
rr(1)=rmin
do
  if (rr(j)>rmax.and.(j>=min_nr.or.min_nr==0)) exit
  if (j>=max_nr.and.max_nr>0) exit
  rr(j+1)=rr(j)*(d1+dg)
  j=j+1
enddo
jrmax=j
k=(rmax-rr(1))/(rr(jrmax)-rr(1))
do j=2,jrmax
  rr(j)=(rr(j)-rr(1))*k+rr(1)
enddo
end

subroutine ga_init_mesh_rcell(rr,nx1,ny1,gg_last,use_symmetr)
!dec$ attributes dllexport:: ga_init_mesh_rcell
!инициализация узлов сетки в полукруглой ячейке для произвольного разбиения по радиусу
!с учетом заданного разбиения по углу для последней точки по радиусу
use gen_mod
integer(4) nx1  !число точек в массиве rr
real(8) rr(nx1) !разбиение по радиусу (может быть неравномерным, может начинаться с 0)
integer(4) ny1  !число точек в массиве gg_last
real(8) gg_last(ny1) !разбиение по углу на последней окружности (от 0 до pi включительно, может бюыть неравномерным)
logical use_symmetr !true - четное число отреззков по углу
integer(4) nx,i,j,k,ngmax,ng_last,ktr,n,ntr
real(8) ds1,ds2,dg,ds,kcos
integer(4), allocatable :: ng(:),ind(:,:),ind_tr(:,:)
real(8), allocatable :: gg(:,:)
real(8), allocatable :: x(:), y(:)
integer(4), allocatable :: tr(:)
real(8) gg_begin,gg_end,gg_delta
ng_last=ny1-1 !число отрезков по углу на последней окружности
nx=nx1-1 !число отрезков по лучу единичного радиуса
gg_begin=gg_last(1)
gg_end=gg_last(ny1)
gg_delta=gg_end-gg_begin
kcos=dcos(pi/6)
allocate(ng(nx1))
do i=1,nx1
  ds1=1.0d10
  ds2=1.0d10
  if (i>1) then
    ds1=rr(i)-rr(i-1)
  endif
  if (i<nx1) then
    ds2=rr(i+1)-rr(i)
  endif
  ds=min(ds1,ds2)/kcos
  ng(i)=gg_delta*rr(i)/ds
  if (use_symmetr.and.mod(ng(i),2)==1) ng(i)=ng(i)+1
enddo
ng(nx1)=ng_last
ngmax=0
do i=1,nx1
  ngmax=max(ngmax,ng(i))
enddo
ngmax=ngmax+1
!количество треугольников и узлов
ntr=ng(1)+ng(nx1)
n=ntr+2
do i=2,nx
  ntr=ntr+ng(i)*2
  n=n+ng(i)+1
enddo

call pg_allocate_areapart_geom(n,ntr,3)
allocate(x(n),y(n),tr(3*ntr))
allocate(gg(nx1,ngmax))
allocate(ind(nx1,ngmax))
allocate(ind_tr(ngmax*2,2))
gg=d0
ind=0
ind_tr=0
k=0
!координаты узлов сетки
do i=1,nx1
  if (rr(i)==d0) then
    dg=d0
  else
    dg=gg_delta/ng(i)
  endif
  do j=1,ng(i)+1
    k=k+1
	if (i==nx1) then
	  gg(i,j)=gg_last(j)
	else
	  gg(i,j)=(j-d1)*dg+gg_begin
	endif
	x(k)=rr(i)*dcos(gg(i,j))
	y(k)=rr(i)*dsin(gg(i,j))
	ind(i,j)=k
  enddo
enddo
!индексы вершин треугольников
k=0
do i=1,nx
  call find_ind_line(gg(i:i+1,1:ngmax),ind_tr,ngmax,ng(i)+1,ng(i+1)+1,ktr)
  do j=1,ktr-1
	  k=k+1
    tr(k)=ind(i,ind_tr(j,1))
	  k=k+1
    tr(k)=ind(i+1,ind_tr(j,2))
	  k=k+1
	  if (ind_tr(j,1)==ind_tr(j+1,1)) then
      tr(k)=ind(i+1,ind_tr(j+1,2))
	  else
	    tr(k)=ind(i,ind_tr(j+1,1))
	  endif
  enddo
enddo
call pg_init_areapart_geom_xy(x,y)
call pg_init_areapart_geom_tr(tr)
deallocate(x,y,tr)
deallocate(ng,gg,ind,ind_tr)
end

subroutine find_ind_line(gg,ind_tr,ngmax,ng1,ng2,k)
use gen_mod
integer(4) ngmax,ng1,ng2,ind_tr(ngmax*2,2)
real(8) gg(2,ngmax),s1,s2
integer(4) i,j,k
ind_tr=0
k=1
i=1
j=1
ind_tr(k,1)=i
ind_tr(k,2)=j
do while(i<ng1.or.j<ng2)
  if (i<ng1) then
    s1=dabs(gg(1,i+1)-gg(2,j))
  else
    s1=1.0d8
  endif
  if (j<ng2) then 
    s2=dabs(gg(1,i)-gg(2,j+1))
  else
    s2=1.0d8
  endif
  if (s1<s2) then
    i=i+1
  else 
    j=j+1
  endif
  k=k+1
  ind_tr(k,1)=i
  ind_tr(k,2)=j
enddo
end

subroutine ga_init_mesh_quadcell_ds(x0,x1,y0,y1,ds,npe)
!dec$ attributes dllexport:: ga_init_mesh_quadcell_ds
!инициализация узлов сетки в прямоугольной области с квадратными ячейками равными ds
use gen_mod
real(8) x0,x1,y0,y1 !координаты ограничивающего прямоугольника
real(8) ds  !шаг сетки
integer(4) npe !3 - треугольники, 4 - четырехугольники
external ftest_empty 
logical ftest_empty !функиция теста вида function ftest(x,y)
call init_mesh_quadcell_ds_test(x0,x1,y0,y1,ds,npe,ftest_empty,.false.)
end

subroutine ga_init_mesh_quadcell_ds_test(x0,x1,y0,y1,ds,npe,ftest)
!dec$ attributes dllexport:: ga_init_mesh_quadcell_ds_test
!инициализация узлов сетки в прямоугольной области с квадратными ячейками равными ds
use gen_mod
real(8) x0,x1,y0,y1 !координаты ограничивающего прямоугольника
real(8) ds  !шаг сетки
integer(4) npe !3 - треугольники, 4 - четырехугольники
external ftest 
logical ftest !функиция теста вида function ftest(x,y)
call init_mesh_quadcell_ds_test(x0,x1,y0,y1,ds,npe,ftest,.true.)
end

subroutine init_mesh_quadcell_ds_test(x0,x1,y0,y1,ds,npe,ftest,have_ftest)
!инициализация узлов сетки в прямоугольной области с квадратными ячейками равными ds
use gen_mod
real(8) x0,x1,y0,y1 !координаты ограничивающего прямоугольника
real(8) ds  !шаг сетки
integer(4) npe !3 - треугольники, 4 - четырехугольники
real(8), allocatable :: xx(:),yy(:)
external ftest
logical ftest,have_ftest
integer(4) nx1,ny1
nx1=(x1-x0)/ds+d5
if (nx1<1) nx1=1
nx1=nx1+1
ny1=(y1-y0)/ds+d5
if (ny1<1) ny1=1
ny1=ny1+1
allocate(xx(nx1),yy(ny1))
call sred_aray(x0,x1,xx,nx1)
call sred_aray(y0,y1,yy,ny1)
call init_mesh_quadcell_test(xx,yy,nx1,ny1,npe,ftest,have_ftest)
deallocate(xx,yy)
end

subroutine ga_init_mesh_quadcell(xx,yy,nx1,ny1,npe)
!dec$ attributes dllexport:: ga_init_mesh_quadcell
!инициализация узлов сетки в прямоугольной области
use gen_mod
integer(4) nx1 !число точек в массиве xx
integer(4) ny1 !число точек в массиве yy
real(8) xx(nx1) !координаты x точек на горизонтальных сторонах
real(8) yy(ny1) !координаты y точек на вертикальных сторонах
integer(4) npe !3 - треугольники, 4 - четырехугольники
external ftest_empty 
logical ftest_empty 
call init_mesh_quadcell_test(xx,yy,nx1,ny1,npe,ftest_empty,.false.)
end

function ftest_empty(x,y)
logical ftest_empty
real(8) x,y
ftest_empty=.true.
x=x
y=y
end

subroutine ga_init_mesh_quadcell_test(xx,yy,nx1,ny1,npe,ftest)
!dec$ attributes dllexport:: ga_init_mesh_quadcell_test
!инициализация узлов сетки в прямоугольной области
!с функцией теста для ячеек
use gen_mod
integer(4) nx1 !число точек в массиве xx
integer(4) ny1 !число точек в массиве yy
real(8) xx(nx1) !координаты x точек на горизонтальных сторонах
real(8) yy(ny1) !координаты y точек на вертикальных сторонах
integer(4) npe !3 - треугольники, 4 - четырехугольники
external ftest 
logical ftest !функиция теста вида function ftest(x,y)
call init_mesh_quadcell_test(xx,yy,nx1,ny1,npe,ftest,.true.)
end

subroutine init_mesh_quadcell_test(xx,yy,nx1,ny1,npe,ftest,have_ftest)
!инициализация узлов сетки в прямоугольной области
use gen_mod
integer(4) nx1 !число точек в массиве xx
integer(4) ny1 !число точек в массиве yy
real(8) xx(nx1) !координаты x точек на горизонтальных сторонах
real(8) yy(ny1) !координаты y точек на вертикальных сторонах
integer(4) i,j,k,k1,nx,ntr,ny,n
real(8), allocatable :: x(:), y(:)
integer(4), allocatable :: tr(:), ind(:)
integer(4) npe !3 - треугольники, 4 - четырехугольники
external ftest
logical ftest,have_ftest
real(8) xc,yc
nx=nx1-1
ny=ny1-1
n=nx1*ny1
ntr=nx*ny
if (npe==3) ntr=ntr*2
allocate(x(n),y(n),tr(npe*ntr))
!индексы вершин треугольников
if (have_ftest) then
  allocate(ind(n))
  ind=0
endif
k=0
do i=1,ny
  if (have_ftest) yc=(yy(i)+yy(i+1))*d5
  do j=1,nx
    if (have_ftest) then
      xc=(xx(j)+xx(j+1))*d5
      if (.not.(ftest(xc,yc))) cycle
    endif
    k1=(i-1)*nx1+j !индекс левого нижнего узла квадратной ячейки
    if (npe==3) then
      k=k+1
      tr(k)=k1
		  k=k+1
      tr(k)=k1+1
		  k=k+1
      tr(k)=k1+nx1
      k=k+1
      tr(k)=k1+1
		  k=k+1
      tr(k)=k1+nx1+1
		  k=k+1
      tr(k)=k1+nx1
    else
      k=k+1
      tr(k)=k1
		  k=k+1
      tr(k)=k1+1
      k=k+1
      tr(k)=k1+nx1+1
		  k=k+1
      tr(k)=k1+nx1
    endif
    if (have_ftest) then
      ind(k1)=-1
      ind(k1+1)=-1
      ind(k1+nx1)=-1
      ind(k1+nx1+1)=-1
    endif
  enddo
enddo    
if (have_ftest) then
  ntr=k/npe
  j=0
  do i=1,n
    if (ind(i)<0) then
      j=j+1
      ind(i)=j
    endif
  enddo
  forall (i=1:k) tr(i)=ind(tr(i))
  !координаты узлов сетки
  k1=0
  n=0
  do i=1,ny1
    do j=1,nx1
      k1=k1+1
      if (ind(k1)>0) then
        n=n+1
        x(n)=xx(j)
  		  y(n)=yy(i)
      endif
    enddo
  enddo
  deallocate(ind)
else 
  k=0
  !координаты узлов сетки
  do i=1,ny1
    do j=1,nx1
      k=k+1
      x(k)=xx(j)
  		y(k)=yy(i)
    enddo
  enddo
endif
call pg_allocate_areapart_geom(n,ntr,npe)
call pg_init_areapart_geom_xy(x,y)
call pg_init_areapart_geom_tr(tr)
deallocate(x,y,tr)
end

subroutine ga_init_mesh_rcell_quads(rr,gg,nr1,ng1)
!dec$ attributes dllexport:: ga_init_mesh_rcell_quads
!инициализация узлов сетки в круговой ячейке четырехугольниками
use gen_mod
integer(4) nr1 !число точек в массиве rr
integer(4) ng1 !число точек в массиве gg
real(8) rr(nr1) !координаты x точек по радиусу
real(8) gg(ng1) !координаты y точек по окружности
integer(4) i,j,k,k1,nr,ntr,ng,n
real(8), allocatable :: x(:), y(:)
integer(4), allocatable :: tr(:)
nr=nr1-1
ng=ng1-1
n=nr1*ng1
ntr=nr*ng
call pg_allocate_areapart_geom(n,ntr,4)
allocate(x(n),y(n),tr(4*ntr))
k=0
!координаты узлов сетки
do i=1,nr1
  do j=1,ng1
    k=k+1
    x(k)=rr(i)*dcos(gg(j))
		y(k)=rr(i)*dsin(gg(j))
  enddo
enddo
!индексы вершин треугольников
k=0
do i=1,nr
  do j=1,ng
    k=k+1
    k1=(i-1)*ng1+j !индекс левого нижнего узла квадратной ячейки
    tr(k)=k1
    k=k+1
    tr(k)=k1+ng1
    k=k+1
    tr(k)=k1+ng1+1
		k=k+1
    tr(k)=k1+1
  enddo
enddo    
call pg_init_areapart_geom_xy(x,y)
call pg_init_areapart_geom_tr(tr)
deallocate(x,y,tr)
end

subroutine ga_init_geom_square(x0,y0,mcellnj,a)
!dec$ attributes dllexport:: ga_init_geom_square
use pgmod
integer(4) mcellnj
integer(4) nj
real(8) x0,y0
real(8) a !половина стороны квадрата
nj=mcellnj/4
if (nj-gs%const%di*2<4) nj=gs%const%di*2+4
call pg_allocate_boundlines(4)
call pg_init_boundline_geomline(1,nj,x0-a,y0-a,x0-a,y0+a)
call pg_init_boundline_geomline(2,nj,x0-a,y0+a,x0+a,y0+a)
call pg_init_boundline_geomline(3,nj,x0+a,y0+a,x0+a,y0-a)
call pg_init_boundline_geomline(4,nj,x0+a,y0-a,x0-a,y0-a)
end

subroutine ga_init_geom_square_symmetr(x0,y0,mcellnj,a)
!dec$ attributes dllexport:: ga_init_geom_square_symmetr
use pgmod
integer(4) mcellnj
integer(4) nj
real(8) x0,y0
real(8) a !половина стороны квадрата
nj=mcellnj/4
if (nj-gs%const%di*2<4) nj=gs%const%di*2+4
nj=(nj+1)/2
call pg_allocate_boundlines(2)
call pg_init_boundline_geomline(1,nj,x0,y0-a,x0-a,y0-a)
call pg_add_boundline_geomline(1,2*nj,x0-a,y0-a,x0-a,y0+a)
call pg_add_boundline_geomline(1,nj,x0-a,y0+a,x0,y0+a)
call pg_init_boundline_geomline(2,nj,x0,y0+a,x0+a,y0+a)
call pg_add_boundline_geomline(2,2*nj,x0+a,y0+a,x0+a,y0-a)
call pg_add_boundline_geomline(2,nj,x0+a,y0-a,x0,y0-a)
end

subroutine ga_init_geom_rectangle(x0,y0,ds,a,b)
!dec$ attributes dllexport:: ga_init_geom_rectangle
use pgmod
real(8) x0,y0 !координаты левого нижнего угла
real(8) a,b !стороны прямоугольника
real(8) ds  !шаг
call pg_allocate_boundlines(4)
call pg_init_boundline_geomlineds(1,ds,x0,y0,x0+a,y0)
call pg_init_boundline_geomlineds(2,ds,x0+a,y0,x0+a,y0+b)
call pg_init_boundline_geomlineds(3,ds,x0+a,y0+b,x0,y0+b)
call pg_init_boundline_geomlineds(4,ds,x0,y0+b,x0,y0)
end

subroutine ga_init_geom_circle(x0,y0,r0,invert_obhod,nj)
!dec$ attributes dllexport:: ga_init_geom_circle
!инициализация цилиндрической геометрии
use pgmod
integer(4) nj
logical invert_obhod
real(8) x0,y0,r0
call ga_init_geom_circle2(x0,y0,r0,invert_obhod,nj,d0)
end

subroutine ga_init_geom_circle2(x0,y0,r0,invert_obhod,nj,gam0)
!dec$ attributes dllexport:: ga_init_geom_circle2
!инициализация цилиндрической геометрии с указанием начальной точки
use pgmod
integer(4) nj
logical invert_obhod
real(8) x0,y0,r0,gam1,gam2,gam0
if (invert_obhod) then
  gam1=pi2+gam0
  gam2=gam0
else
  gam1=gam0
  gam2=pi2+gam0
endif
call pg_allocate_boundlines(1)
call pg_init_boundline_geomcircle(1,nj,x0,y0,r0,gam1,gam2)
end

subroutine ga_init_geom_ellipse(x0,y0,rx,ry,g0,invert_obhod,nj)
!dec$ attributes dllexport:: ga_init_geom_ellipse
!инициализация цилиндрической геометрии
use pgmod
integer(4) nj
logical invert_obhod
real(8) x0,y0,rx,ry,g0,gam1,gam2
if (invert_obhod) then
  gam1=pi2
  gam2=d0
else
  gam1=d0
  gam2=pi2
endif
call pg_allocate_boundlines(1)
call pg_init_boundline_geom_ellipse(1,nj,x0,y0,rx,ry,g0,gam1,gam2)
end

subroutine ga_init_geom_circle_symmetr(x0,y0,r0,invert_obhod,nj)
!dec$ attributes dllexport:: ga_init_geom_circle_symmetr
!инициализация цилиндрической геометрии
use pgmod
integer(4) nj
logical invert_obhod
real(8) x0,y0,r0
call ga_init_geom_ellipce_symmetr(x0,y0,r0,r0,d0,invert_obhod,nj)
end

subroutine ga_init_geom_ellipce_symmetr(x0,y0,rx,ry,g0,invert_obhod,nj)
!dec$ attributes dllexport:: ga_init_geom_ellipce_symmetr
!инициализация цилиндрической геометрии
use pgmod
integer(4) nj,nj2
logical invert_obhod
real(8) x0,y0,rx,ry,g0,gam1,gam2,gam3
nj2=(nj+1)/2
if (invert_obhod) then
  gam1=pi*1.5d0
  gam2=pi5
  gam3=-pi5
else
  gam1=-pi5
  gam2=pi5
  gam3=pi*1.5d0
endif
call pg_allocate_boundlines(2)
call pg_init_boundline_geom_ellipse(1,nj2,x0,y0,rx,ry,g0,gam1,gam2)
call pg_init_boundline_geom_ellipse(2,nj2,x0,y0,rx,ry,g0,gam2,gam3)
end

subroutine ga_init_geom_cylider_particle(x0,y0,r0,type_ga)
!dec$ attributes dllexport:: ga_init_geom_cylider_particle
!инициализация цилиндрической геометрии с использованием boundLineType=2
!для задачи Стокса с обходом по часовой стрелке
use pgmod
real(8) x0,y0,r0
integer(4) type_ga !0 - 6 неизвестных, 1 - 4 неизвестных
call ga_init_geom_cylider_particle2(x0,y0,r0,type_ga,1,1)
end

subroutine ga_init_geom_cylider_particle2(x0,y0,r0,type_ga,nsn,n_shift)
!dec$ attributes dllexport:: ga_init_geom_cylider_particle2
!инициализация цилиндрической геометрии с использованием boundLineType=2
!с учетом сдвига точек коллокации
!для задачи Стокса с обходом по часовой стрелке
use pgmod
real(8) x0,y0,r0
integer(4) type_ga !0 - 6 неизвестных, 1 - 4 неизвестных
integer(4) nsn !количество вариантов точек коллокаций
integer(4) n_shift !номер массива точек коллокаций
integer(4) n
type(TBoundline2), pointer :: b
integer(4) i
call pg_allocate_boundlines2(1)
call pg_init_cylider_particle(1,x0,y0,r0,d0,-1,type_ga)
!***точки каллокации
n=3-type_ga
if (nsn==2) then
  call init_gsMain_s_ndns2(n)
else
  call init_gsMain_s_ndns(n)
endif
b=>gsbnd%line2(1)
do i=1,gsarea%nu
  call pg_allocate_gc2(1,i,n)
  b%gc(i)%s_gu=>gsMain%sn(n_shift)%s_ndns(i,:)
enddo
end

subroutine ga_GetFourierBounLineFuncApprox(psiom,n,ga)
!dec$ attributes dllexport:: ga_GetFourierBounLineFuncApprox
!для круговой частицы !получить первые три коэффициента разложения в ряд Фурье
use pgmod
type(TBounLineFuncApprox) ga
integer(4) n
real(8) COEF(n),psiom(n)
call ga_GetFourierBounLineFuncCoef(psiom,n,coef)
!восстановление функции по первым трем членам ряда
ga%kk=d1
ga%val(1)=COEF(1)*d5
ga%val(2)=-COEF(3)
ga%val(3)=COEF(2)
ga%valkk=ga%val
!вычисление полей поэффициентов и углов осей диполей
ga%cc=dsqrt(COEF(2)**2+COEF(3)**2)
ga%delta=datan2(-COEF(2),COEF(3))
end

subroutine ga_GetFourierBounLineFuncCoef(psiom,n,coef)
!dec$ attributes dllexport:: ga_GetFourierBounLineFuncCoef
!для круговой частицы !получить первые три коэффициента разложения в ряд Фурье
use pgmod
integer(4) j,n
real(8) ff(n),COEF(n),psiom(n)
forall (j=2:n) ff(j) = (psiom(j-1)+psiom(j))*d5
ff(1) = (psiom(n)+psiom(1))*d5
!разложение в ряд Фурье
CALL dFFTRF (n, ff, COEF)
forall (j=1:n) COEF(j) = d2*COEF(j)/n
end

subroutine ga_GetFourierBounLineFuncApprox_testhalf(psiom,n,ga,is_half)
!dec$ attributes dllexport:: ga_GetFourierBounLineFuncApprox_testhalf
!для круговой частицы получить первые три коэффициента разложения в ряд Фурье
!с половины круга при условии симметрии
use pgmod
type(TBounLineFuncApprox) ga
integer(4) j,n,n2
real(8) psiom(n)
real(8), allocatable :: psiom2(:)
logical is_half
if (is_half) then
  n2=n*2
  allocate(psiom2(n2))
  forall (j=1:n) psiom2(j)=psiom(j)
  forall (j=1:n) psiom2(j+n)=-psiom(n-j+1)
  call ga_GetFourierBounLineFuncApprox(psiom2,n2,ga)
  deallocate(psiom2)
else
  call ga_GetFourierBounLineFuncApprox(psiom,n,ga)
endif
end

subroutine ga_mesh_square(rmin,rmax,ng,rr,gg)
!dec$ attributes dllexport:: ga_mesh_square
!преобразование треугольной сетки, заполняющей круговую ячейку к треугольной, заполняющей прямоугольную ячейку
use pgmod
real(8) rmin !минимальный радиус исходной сетки
real(8) rmax !максимальный радиус исходной сетки
integer(4) ng !число отрезков (число элементов - 1) для массивов rr и gg
real(8) rr(ng+1) !массив r для внешней границы для каждого gg
real(8) gg(ng+1) !разбиение по gamma
integer(4) i,n,pg_get_int
real(8) bb(ng+1),cc(4,ng+1),g,rmax2,a1,b1,d,r2
real(8), allocatable :: x(:), y(:)
call spline_linear(ng+1,gg,rr,bb,cc)
n=pg_get_int(6,2)
allocate(x(n),y(n))
call pg_get_array_real(6,3,x,n)
call pg_get_array_real(6,4,y,n)
d=rmax-rmin
do i=1,n
  g=datan2(y(i),x(i))
  if (g<d0) g=g+pi2
  r2=dsqrt(x(i)**2+y(i)**2)
  rmax2=dcsval(g,ng,bb,cc)
  a1=(rmax2-rmin)/d
  b1=rmin*(d1-a1)
  r2=a1*r2+b1
  x(i)=r2*dcos(g)
  y(i)=r2*dsin(g)
enddo
call pg_init_areapart_geom_xy(x,y)
deallocate(x,y)
end

subroutine ga_drw_trmesh(ia)
!dec$ attributes dllexport:: ga_drw_trmesh
!вывод в файл сетки
use pgmod
integer(4) i,ia,k,i1,i2,i0,nb_all
type(TArea), pointer :: ar
type(TBound), pointer :: b
type(TBoundline2), pointer :: bl2
external ga_drw_trmesh_write_point
OPEN (1,FILE='trmesh.dat')
write(1,*) 'TITLE = "Triangle Mesh"'
write(1,*) 'VARIABLES = "X", "Y"'
call ga_drw_trmesh_func_zone2(ia,1,"area_e",ga_drw_trmesh_write_point)
i1=ia
i2=ia
if (ia==0) then
  i1=1
  i2=gs%na
endif
nb_all=0
do i0=i1,i2
  ar=>gs%a(i0)
  do k=1,ar%nb
    b=>ar%bnd(k)
    if (b%boundLineType==2) then
      nb_all=nb_all+1
    else
      nb_all=nb_all+b%npanel+1
    endif
  enddo
enddo 
if (nb_all==0) then
  write(1,"('ZONE T=""panels"", I=', i0, ', F=POINT')") 1
  call ga_drw_trmesh_write_point(1,1)
else
  write(1,"('ZONE T=""panels"", I=', i0, ', F=POINT')") nb_all
  do i0=i1,i2
    ar=>gs%a(i0)
    do k=1,ar%nb
      b=>ar%bnd(k)
      if (b%boundLineType==2) then
        bl2=>b%line2(1)
        WRITE(1,"(F9.5, ' ', F9.5)") dreal(bl2%zc), dimag(bl2%zc)
      else
        do i=1,b%npanel+1
          WRITE(1,"(F9.5, ' ', F9.5)") b%x(i), b%y(i)
        enddo
      endif
    enddo
  enddo
endif
close(1)
!trmesh.lay
end

subroutine ga_drw_trmesh_write_point(ifile,i)
use pgmod
integer(4) ifile,i
type(TBound), pointer :: b
if (i==0) then
  b=>gsarea%bnd(1)
  WRITE(ifile,"(F9.5, ' ', F9.5)") b%x(1), b%y(1)
else
  WRITE(ifile,"(F9.5, ' ', F9.5)") dreal(gsarea%a%zm(i)), dimag(gsarea%a%zm(i))
endif
end

subroutine ga_drw_trmesh_func(filename,vars,write_point)
!dec$ attributes dllexport:: ga_drw_trmesh_func
!вывод в файл сетки и функции текущей области
use pgmod
character(*) vars,filename
external write_point
call ga_drw_trmesh_func_start(1,filename,vars)
call ga_drw_trmesh_func_zone(1,"area_e",write_point)
close(1)
end

subroutine ga_drw_trmesh_func_start(ifile,filename,vars)
!dec$ attributes dllexport:: ga_drw_trmesh_func_start
!вывод в файл сетки и функции
!сохранение заголовка tecplot файла
use pgmod
integer(4) ifile
character(*) vars,filename
OPEN (ifile,FILE=filename)
write(ifile,*) 'TITLE = "Triangle Mesh"'
write(ifile,*) 'VARIABLES = ',vars
end

subroutine ga_drw_trmesh_func_zone(ifile,zonename,write_point)
!dec$ attributes dllexport:: ga_drw_trmesh_func_zone
!вывод в файл сетки и функции текущей области как одной зоны
use pgmod
integer(4) ifile
character(*) zonename
external write_point
call ga_drw_trmesh_func_zone2(gsarea%i,ifile,zonename,write_point)
end

subroutine ga_drw_trmesh_func_zone2(ia,ifile,zonename,write_point)
!dec$ attributes dllexport:: ga_drw_trmesh_func_zone
!вывод в файл сетки и функции всех областей текущей задачи как одной зоны
!!!функция меняет gsarea - выбирая текущую область
!write_point(ifile,i) i - номер узла сетки текущей области
use pgmod
integer(4) i,k,ifile
integer(4) ia !номер области или 0 - для вывода всех областей в одной зоне
character(*) zonename !имя зоны
character(200) formatstr
type(areatype), pointer :: a
integer(4) n_all,ntr_all,ns,ntrs,i1,i2,i0
type(TBound), pointer :: b
integer(4), allocatable :: nShift(:),ntrShift(:)
logical is_quad
if (ia==0) then
  n_all=0
  ntr_all=0
  allocate(nShift(gs%na),ntrShift(gs%na))
  is_quad=.false.
  do i=1,gs%na
    a=>gs%a(i)%a
    n_all=n_all+a%n
    ntr_all=ntr_all+a%ntr
    if (i==1) then
      nShift(1)=0
      ntrShift(1)=0
    else
      nShift(i)=nShift(i-1)+ns
      ntrShift(i)=ntrShift(i-1)+ntrs
    endif
    ns=a%n
    ntrs=a%ntr
    if (a%n>0.and.a%npe==4) is_quad=.true.
  enddo
  i1=1
  i2=gs%na
else
  a=>gs%a(ia)%a
  n_all=a%n
  ntr_all=a%ntr
  allocate(nShift(1),ntrShift(1))
  nShift(1)=0
  ntrShift(1)=0
  is_quad=a%npe==4
  i1=ia
  i2=ia
endif
if (n_all>0) then
  if (is_quad) then
    formatstr="('ZONE T="""//trim(zonename)//""" N=', i0, ', E=', i0, ', F=FEPOINT, ET=QUADRILATERAL')"
  else
    formatstr="('ZONE T="""//trim(zonename)//""" N=', i0, ', E=', i0, ', F=FEPOINT, ET=TRIANGLE')"
  endif
  write(ifile,trim(formatstr)) n_all, ntr_all
  do i0=i1,i2
    call pg_bind_domain(i0)
    do i=1,gsarea%a%n
      call write_point(ifile,i)
    enddo
  enddo
  do i0=i1,i2
    a=>gs%a(i0)%a
    if (is_quad.and.a%npe==3) then
      do i=1,a%ntr
        WRITE(1,"(4(' ',i0))") (a%trm(k,i)+nShift(i0),k=1,a%npe),a%trm(3,i)+nShift(i0)
      enddo
    else
      do i=1,a%ntr
        WRITE(1,"(4(' ',i0))") (a%trm(k,i)+nShift(i0),k=1,a%npe)
      enddo
    endif
  enddo
else
  b=>gs%a(i1)%bnd(1)
  write(1,"('ZONE T="""//trim(zonename)//""" N=', i0, ', E=', i0, ', F=FEPOINT, ET=TRIANGLE')") 3,1
  do i=1,3
    call write_point(ifile,0)
  enddo
  WRITE(1,"(4(' ',i0))") 1,2,3
endif
deallocate(nShift,ntrShift)
end

subroutine gaintf_prepare(nn,x1,y1,x2,y2,ds,xx,yy,mode,ifs)
use pgmod
integer(4) nn,j
real(8) x1,y1,x2,y2,ds
real(8),target :: xx(nn),yy(nn)
integer(4) mode !1 - заданы концы x1,y1,x2,y2,ds
                !2 - горизонтальная линия x,y1,n
                !3 - линия массивом координат x,y,n
type(intf_struct), target :: ifs
real(8) dx,dy
complex(8) zder
if (mode==1) then
  call line_array_ds(x1,y1,x2,y2,ds,ifs%x0,ifs%y0,ifs%s0,ifs%n,4)
  ifs%x=>ifs%x0
  ifs%y=>ifs%y0
  ifs%s=>ifs%s0
else
  ifs%x=>xx
  ifs%n=nn
  if (mode==2) then
    allocate(ifs%y0(nn))
    ifs%y=>ifs%y0
    ifs%y=y1
  else
    ifs%y=>yy
  endif
  allocate(ifs%s(nn))
endif
ifs%xder=ifs%x(ifs%n)-ifs%x(1)
ifs%yder=ifs%y(ifs%n)-ifs%y(1)
zder=dcmplx(ifs%xder,ifs%yder)
zder=zder/cdabs(zder)*ii
ifs%xder=dreal(zder)
ifs%yder=dimag(zder)
if (mode>1) then
  if (ifs%xder==d0) then
    if (ifs%x(ifs%n)>ifs%x(1)) then
      ifs%s=ifs%x-ifs%x(1)
    else
      ifs%s=-(ifs%x-ifs%x(1))
    endif
  elseif (ifs%yder==d0) then
    if (ifs%y(ifs%n)>ifs%y(1)) then
      ifs%s=ifs%y-ifs%y(1)
    else
      ifs%s=-(ifs%y-ifs%y(1))
    endif
  else
    ifs%s(1)=d0
    do j=2,ifs%n
      dx=ifs%x(j)-ifs%x(j-1)
      dy=ifs%y(j)-ifs%y(j-1)
      ifs%s(j)=ifs%s(j-1)+dsqrt(dx**2+dy**2)
    enddo
  endif
endif
end

subroutine gaintf_free(ifs)
use pgmod
type(intf_struct) ifs
if (allocated(ifs%x0)) deallocate(ifs%x0)
if (allocated(ifs%y0)) deallocate(ifs%y0)
if (allocated(ifs%s0)) deallocate(ifs%s0)
end

function ga_intf_oneline_base(ifs,nf)
!dec$ attributes dllexport:: ga_intf_oneline_base
!нахождение интеграла от функции вдоль прямой линии заданной массивами координат
use pgmod
type(intf_struct) ifs
integer(4) nf !номер функции
              !1 - psi
              !2 - dpsi/dn
              !3 - Delta psi
              !4 - d(Delta  psi)/dn
              !n - нормаль, получаемая поворотом касательной против часовой стрелки
integer(4) j
real(8) ga_intf_oneline_base,ff(ifs%n),pg_get_fun_xy,bb(ifs%n),cc(4,ifs%n)
ff(1)=pg_get_fun_xy(ifs%x(1),ifs%y(1),nf,ifs%xder,ifs%yder,3)
do j=2,ifs%n-1
  ff(j)=pg_get_fun_xy(ifs%x(j),ifs%y(j),nf,ifs%xder,ifs%yder,1)
enddo
ff(ifs%n)=pg_get_fun_xy(ifs%x(ifs%n),ifs%y(ifs%n),nf,ifs%xder,ifs%yder,3)
call dcsakm(ifs%n,ifs%s,ff,bb,cc)
ga_intf_oneline_base=dcsitg(ifs%s(1),ifs%s(ifs%n),ifs%n-1,bb,cc)
end

function ga_intf_oneline(x,y,n,nf)
!dec$ attributes dllexport:: ga_intf_oneline
!нахождение интеграла от функции вдоль прямой линии заданной массивами координат
use pgmod
integer(4) n !число точек в массиве линии
real(8) x(n),y(n) !массив координат линии (отрезок прямой)
integer(4) nf !номер функции
real(8) ga_intf_oneline,ga_intf_oneline_base
type(intf_struct) ifs
call gaintf_prepare(n,d0,d0,d0,d0,d0,x,y,3,ifs)
ga_intf_oneline=ga_intf_oneline_base(ifs,nf)
call gaintf_free(ifs)
end

function ga_intf_oneline_x(x,y,n,nf)
!dec$ attributes dllexport:: ga_intf_oneline_x
!нахождение интеграла от функции вдоль горизонтальной линии
use pgmod
integer(4) n !число точек в массиве линии
real(8) x(n) !массив x координат линии (отрезок прямой)
real(8) y !y координатf линии (отрезок прямой)
real(8) ga_intf_oneline_x,ga_intf_oneline_base
integer(4) nf !номер функции
type(intf_struct) ifs
call gaintf_prepare(n,d0,y,d0,d0,d0,x,y,2,ifs)
ga_intf_oneline_x=ga_intf_oneline_base(ifs,nf)
call gaintf_free(ifs)
end

function ga_intf_oneline_ds(x1,y1,x2,y2,ds,nf)
!dec$ attributes dllexport:: ga_intf_oneline_ds
!нахождение интеграла от функции по линии (заданы концы отрезка)
use pgmod
real(8) x1,y1,x2,y2 !координаты концов отрезка
real(8) ds !шаг интегрирования
integer(4) nf !номер функции
real(8) ga_intf_oneline_ds,ga_intf_oneline_base
type(intf_struct) ifs
call gaintf_prepare(0,x1,y1,x2,y2,ds,x1,y1,1,ifs)
ga_intf_oneline_ds=ga_intf_oneline_base(ifs,nf)
call gaintf_free(ifs)
end

function gagradp_oneline_base(ifs) result(dp)
use pgmod
real(8) ga_intf_oneline_base,dp,dp1,dp2
type(intf_struct) ifs
select case(gsarea%type_eq(1))
case(1,3,15)
case default
  call gs_print_stop("Error ga_intf_oneline!")
endselect
dp=0
if (gsarea%type_eq(1)/=1) then
  dp1=-ga_intf_oneline_base(ifs,4)
  if (gsarea%type_eq(1)==15) dp1=dp1*gsarea%const%mub
  dp=dp1
endif
if (gsarea%type_eq(1)/=3) then
  dp2=ga_intf_oneline_base(ifs,2)/gsarea%const%k
  dp=dp+dp2
endif
call gaintf_free(ifs)
end

function ga_gradp_oneline(x,y,n)
!dec$ attributes dllexport:: ga_gradp_oneline
!нахождение перепада давления вдоль прямой линии заданной массивами координат
use pgmod
integer(4) n !число точек в массиве линии
real(8) x(n),y(n) !массив координат линии (отрезок прямой)
real(8) ga_gradp_oneline,gagradp_oneline_base
type(intf_struct) ifs
call gaintf_prepare(n,d0,d0,d0,d0,d0,x,y,3,ifs)
ga_gradp_oneline=gagradp_oneline_base(ifs)
end

function ga_gradp_oneline_x(x,y,n)
!dec$ attributes dllexport:: ga_gradp_oneline_x
!нахождение перепада давления по горизонтальной линии
use pgmod
integer(4) n !число точек в массиве линии
real(8) x(n) !массив x координат линии (отрезок прямой)
real(8) y !y координатf линии (отрезок прямой)
real(8) ga_gradp_oneline_x,gagradp_oneline_base
type(intf_struct) ifs
call gaintf_prepare(n,d0,y,d0,d0,d0,x,y,2,ifs)
ga_gradp_oneline_x=gagradp_oneline_base(ifs)
end

function ga_gradp_oneline_ds(x1,y1,x2,y2,ds)
!dec$ attributes dllexport:: ga_gradp_oneline_ds
!нахождение перепада давления по линии (заданы концы отрезка)
use pgmod
real(8) x1,y1,x2,y2 !координаты концов отрезка
real(8) ds !шаг интегрирования
real(8) ga_gradp_oneline_ds,gagradp_oneline_base
type(intf_struct) ifs
call gaintf_prepare(0,x1,y1,x2,y2,ds,x1,y1,1,ifs)
ga_gradp_oneline_ds=gagradp_oneline_base(ifs)
end

subroutine ga_get_force_bl2(ib,fx,fy)
!dec$ attributes dllexport:: ga_get_force_bl2
!Вычисление сил, действующих на круговую частицу BoundLineType=2
use pgmod
integer(4) ib !номер границы
real(8) fx,fy !компоненты вектора силы
type(TBound), pointer :: b
type(TBoundline2), pointer :: bl2
type(TBounLineFuncApprox), pointer :: ga,ga1
complex(8) f
b=>gsarea%bnd(ib)
if (b%boundLineType/=2) call gs_print_stop('Error ga_get_force_bl2! boundLineType/=2')
bl2=>b%line2(1)
ga=>bl2%ga(3)
ga1=>bl2%ga(4)
f=-pi*bl2%r*(ga%val(2)+bl2%r*ga1%val(2)+ii*(ga%val(3)+bl2%r*ga1%val(3)))
if (bl2%g0/=d0) f=f*cdexp(ii*bl2%g0)
fx=dreal(f)
fy=dimag(f)
end

subroutine ga_get_force_circle(ib,fx,fy)
!dec$ attributes dllexport:: ga_get_force_circle
!Вычисление сил, действующих на круговую частицу BoundLineType=1
use pgmod
integer(4) ib !номер границы
real(8) fx,fy !компоненты вектора силы
call ga_get_force_circle_bl2(ib,1,gsarea%bnd(ib)%nline,fx,fy)
end

subroutine ga_get_force_circle_bl(ib,ibl,fx,fy)
!dec$ attributes dllexport:: ga_get_force_circle_bl
!Вычисление сил, действующих на круговой участок границы BoundLineType=1 (целый или половинку)
!для половинки круга вычисляется сила, действующая на половину
use pgmod
integer(4) ib !номер границы
integer(4) ibl !номер участка границы
real(8) fx,fy !компоненты вектора силы
call ga_get_force_circle_bl2(ib,ibl,ibl,fx,fy)
end

subroutine ga_get_force_circle_bl2(ib,ibl1,ibl2,fx,fy)
!dec$ attributes dllexport:: ga_get_force_circle_bl
!Вычисление сил, действующих на круговой группу участоков границы BoundLineType=1 (целый или половинку)
!для половинки круга вычисляется сила, действующая на половину
use pgmod
integer(4) ib !номер границы
integer(4) ibl1,ibl2 !номер первого и последнего участка границы
real(8) fx,fy !компоненты вектора силы
type(TBoundline_geomdetale), pointer :: gd
type(TBound), pointer :: b
type(TBoundline), pointer :: bl1,bl2
type(TBounLineFuncApprox) ga,ga1,ga2
complex(8) f
real(8) g0,dgam
logical is_inner_cont !контур внутренний
logical is_half !половина круга при условии симметрии по диаметру
real(8) f_half
integer(4) i_half,npanel
logical test_boundlines_circle
b=>gsarea%bnd(ib)
if (b%boundLineType/=1) call gs_print_stop('Error ga_get_force_circle_bl! boundLineType/=1')
if (.not.test_boundlines_circle(ib,ibl1,ibl2,dgam)) call gs_print_stop('Error ga_get_force_circle_bl! .not.test_boundlines_circle')
bl1=>b%line(ibl1)
bl2=>b%line(ibl2)
gd=>b%geom_detale0(bl1%igd_begin)
is_inner_cont=dgam<d0
f_half=dabs(dgam)/pi
i_half=nint(f_half)
if (dabs(f_half-i_half)>1d-6) i_half=-1
select case (i_half)
case (1)
  is_half=.true.
case (2)
  is_half=.false.
case default
  call gs_print_stop('Error ga_get_force_circle_bl! Test half')
endselect
f=c0
npanel=bl2%i_end-bl1%i_begin+1
select case (gsarea%type_eq(1))
case (3,15)
  call pg_allocate_and_set_ga2_base(ga,3)
  call pg_allocate_and_set_ga2_base(ga1,3)
  call ga_GetFourierBounLineFuncApprox_testhalf(b%psiom(bl1%i_begin:bl2%i_end,3),npanel,ga,is_half)
  call ga_GetFourierBounLineFuncApprox_testhalf(b%psiom(bl1%i_begin:bl2%i_end,4),npanel,ga1,is_half)
  if (is_inner_cont) then
    ga%val(2)=-ga%val(2)
    ga1%val(3)=-ga1%val(3)
  endif
  f=ga%val(2)-gd%r*ga1%val(2)-ii*(ga%val(3)-gd%r*ga1%val(3))
endselect
select case (gsarea%type_eq(1))
case (1,15)
  call pg_allocate_and_set_ga2_base(ga2,3)
  call ga_GetFourierBounLineFuncApprox_testhalf(b%psiom(bl1%i_begin:bl2%i_end,2),npanel,ga2,is_half)
  if (is_inner_cont) ga2%val(3)=-ga2%val(3)
  if (gsarea%type_eq(1)==15) f=f*gsarea%const%mub
  f=f+gd%r*(ga2%val(2)-ii*ga2%val(3))/gsarea%const%k
endselect
f=f*pi*gd%r
g0=gd%gam1 !!!!!!!!!!!! дописать вычисление
if (g0/=d0) f=f*cdexp(ii*g0)
if (is_half) f=f*d5
fx=dreal(f)
fy=dimag(f)
select case (gsarea%type_eq(1))
case (3,15)
  call pg_deallocate_ga_base(ga)
  call pg_deallocate_ga_base(ga1)
endselect
select case (gsarea%type_eq(1))
case (1,15)
  call pg_deallocate_ga_base(ga2)
endselect
end

function test_boundlines_circle(ib,ibl1,ibl2,dgam) result (is_circ)
use pgmod
integer(4) ib,ibl1,ibl2,i,j
logical is_circ
real(8), parameter :: eps=1.0d-8
type(TBoundline_geomdetale), pointer :: gd,gd_prev
type(TBound), pointer :: b
type(TBoundline), pointer :: bl
real(8) dg,dgam
logical eq_eps,ceq_eps
b=>gsarea%bnd(ib)
is_circ=.false.
gd_prev=>null()
dgam=d0
do i=ibl1,ibl2
  bl=>b%line(i)
  do j=bl%igd_begin,bl%igd_end
    gd=>b%geom_detale0(j)
    if (gd%mode/=3) return
    if (associated(gd_prev)) then
      if (.not.eq_eps(gd%x,gd_prev%x,eps)) return
      if (.not.eq_eps(gd%y,gd_prev%y,eps)) return
      if (.not.eq_eps(gd%r,gd_prev%r,eps)) return
      if (.not.ceq_eps(cdexp(ii*gd%gam1),cdexp(ii*gd_prev%gam2),eps)) return
    endif
    dg=gd%gam2-gd%gam1
    if (dgam/=d0) then
      if (dg*dgam<d0) return
    endif
    dgam=dgam+dg
    gd_prev=>gd
  enddo
enddo
is_circ=.true.
end

subroutine ga_get_force_circle_bound(ib,fx,fy)
!dec$ attributes dllexport:: ga_get_force_circle_bound
!Вычисление сил, действующих на круговую границу с проверкой типа границы
use pgmod
integer(4) ib !номер границы
real(8) fx,fy !компоненты вектора силы
select case (gsarea%bnd(ib)%boundLineType)
case (1)
  call ga_get_force_circle(ib,fx,fy)
case (2)
  call ga_get_force_bl2(ib,fx,fy)
case default 
  call gs_print_stop('Error ga_get_force_circle_bound!')
endselect
end

function ga_get_psi_particle(ib) result(psi1)
!dec$ attributes dllexport:: ga_get_psi_particle
!вернуть значение функции тока на границе-частице
use pgmod
integer(4) ib !номер границы
real(8) psi1
type(TBounLineFuncApprox), pointer :: ga_psi
type(TBound), pointer :: b
type(TBoundline2), pointer :: bl2
b=>gsarea%bnd(ib)
select case (b%boundLineType) 
case(1)
  psi1=b%psiom(1,1)
case(2)
  bl2=>b%line2(1)
  ga_psi=>bl2%ga(1)
  psi1=ga_psi%valkk(1)
endselect
end

function ga_get_om_particle(ib) result(om1)
!dec$ attributes dllexport:: ga_get_psi_particle
!вернуть значение функции eta на границе-частице (eta=-om)
use pgmod
integer(4) ib !номер границы
real(8) om1
type(TBounLineFuncApprox), pointer :: ga_om
type(TBound), pointer :: b
type(TBoundline2), pointer :: bl2
integer(4) i
b=>gsarea%bnd(ib)
select case (b%boundLineType) 
case(1)
  om1=d0
  do i=1,b%npanel
    om1=om1+b%psiom(i,3)*b%l(i)
  enddo
  om1=om1/b%s(b%npanel+1)
case(2)
  bl2=>b%line2(1)
  ga_om=>bl2%ga(3)
  om1=ga_om%valkk(1)
endselect
end

function ga_gradp_boundline(ib,ibl) result (f)
!dec$ attributes dllexport:: ga_gradp_boundline
!вычисление перепада давления по участку границы
!с учетом знака: f=p_end-p_begin
!для течения слева на право
use pgmod
integer(4) ib  !номер границы
integer(4) ibl !номер участка границы
real(8) f
integer(4) i
type(TBound), pointer :: b
type(TBoundline), pointer :: bl
b=>gsarea%bnd(ib)
bl=>b%line(ibl)
f=d0
select case (gsarea%type_eq(1))
case (3)
  do i=bl%i_begin,bl%i_end
    f=f-b%psiom(i,4)*b%l(i)
  enddo
case default
  call gs_print_stop('Error ga_gradp_boundline! type_equation/=3')
endselect
end

subroutine ga_readMeshMSh(fname)
!dec$ attributes dllexport:: ga_readMeshMSh
!считывает в текущий участок треугольную сетку из файла фармата *.msh
use pgmod
character(*) fname !имя файла
character(200) line
integer(4) mode,i,j,n,k
real(8) x,y,z
integer(4), allocatable :: tr(:,:)
integer(4) tmp(8)
mode=0
OPEN (1,FILE=fname)
do while (.not.eof(1))
  read(1,*) line
  if (trim(line)=='$Nodes') then
    mode=1
  elseif (trim(line)=='$Elements') then
    mode=2
  endif
  selectcase(mode)
  case (1) !nodes
    read(1,*) j
    call pg_allocate_areapart_geom_n(j)
    do i=1,j
      read(1,*) j,x,y,z
      gsareapart%zm(i)=dcmplx(x,y)
    enddo
    mode=0
  case (2) !elements
    read(1,*) j
    allocate(tr(3,j))
    tr=0
    n=0
    do i=1,j
      read(1,*) (tmp(k),k=1,2),(tmp(k),k=3,merge(8,0,tmp(2)==2))
      if (tmp(2)==2) then
        n=n+1
        forall (j=1:3) tr(j,n)=tmp(j+5)
      endif
    enddo
    call pg_allocate_areapart_geom_tr(n,3)
    do i=1,n
      forall (j=1:3) gsareapart%trm(j,i)=tr(j,i)
    enddo
    mode=0
  endselect
enddo
close(1)
deallocate(tr)
end

subroutine ga_readMeshVTK(fname)
!dec$ attributes dllexport:: ga_readMeshVTK
!считывает в текущий участок треугольную сетку из файла фармата *.vtk
!!!читает только треугольную сетку
use pgmod
character(*) fname !имя файла
character(200) line,word
integer(4) i,n,cl(4),j
logical split_str,res,point_read,cells_read
real(8) x,y,z
point_read=.false.
cells_read=.false.
OPEN (1,FILE=fname)
read(1,*) line !# vtk DataFile Version 3.0
read(1,*) line !vtk output
read(1,*) line
if (trim(line)/='ASCII') call gs_print_stop("Error ga_readMeshVTK! not ASCII")
do while (.not.eof(1))
  read(1,'(A)') line
  i=1
  res=split_str(line,i," ",word)
  if (trim(word)=='DATASET') then
    res=split_str(line,i," ",word)
    if (trim(word)/='UNSTRUCTURED_GRID') call gs_print_stop("Error ga_readMeshVTK! not UNSTRUCTURED_GRID") 
  elseif (trim(word)=='POINTS') then
    res=split_str(line,i," ",word)
    read(word,*) n
    call pg_allocate_areapart_geom_n(n)
    do i=1,n
      read(1,*) x,y,z
      gsareapart%zm(i)=dcmplx(x,y)
    enddo
    point_read=.true.
  elseif (trim(word)=='CELLS') then
    res=split_str(line,i," ",word)
    read(word,*) n
    call pg_allocate_areapart_geom_tr(n,3)
    do i=1,n
      read(1,*) (cl(j),j=1,4)
      if (cl(1)/=3) call gs_print_stop("Error ga_readMeshVTK! not triangle cell") 
      forall (j=1:3) gsareapart%trm(j,i)=cl(j+1)+1
    enddo
    cells_read=.true.
  endif
  if (point_read.and.cells_read) exit
  enddo
close(1)
end
subroutine init_mesh3(xmax,ymax,nj,ny)
!инициализация узлов сетки в прямоугольной области
use pgmod2
integer(4) nj,ny,nj1,ny1
real(8) xmax,ymax
nj1=nj
ny1=ny
call init_mesh3_1(d0,xmax,d0,ymax,nj1,ny1)
end

subroutine init_mesh3_1(xmin,xmax,ymin,ymax,nj,ny)
!инициализация узлов сетки в прямоугольной области со смещением
use pgmod2
integer(4) i
real(8) xx(nmax),yy(nmax),xmax,ymax,xmin,ymin,dx,dy
integer(4) nj,ny,nj1,ny1
nj1=nj
ny1=ny
dx=(xmax-xmin)/nj1
dy=(ymax-ymin)/ny1
do i=1,ny1+1
  yy(i)=(i-d1)*dy+ymin
enddo
if (need_thickening) call add_thickening(yy,ny1,n_thickening,0)
do i=1,nj1+1
  xx(i)=(i-d1)*dx+xmin
enddo
call ga_init_mesh_quadcell(xx,yy,nj1+1,ny1+1,npe_)
end

subroutine add_thickening(xx,nj,n,mode)
use pgmod2
real(8), parameter :: lam=1.5d0 
integer(4) nj,n,nj1,i
real(8) xx(nmax),xx1(nmax),sred,x
integer(4) mode !0-внчале, 1-вконце
nj1=nj+n
if (mode==0) then
  xx1(1)=xx(1)
  do i=1,n
    x=((i-d0)/(n+d1))**lam
    xx1(i+1)=sred(d0,d1,xx(1),xx(2),x)
  enddo
  xx1(2+n:nj+1+n)=xx(2:nj+1)
else
  xx1(nj1+1)=xx(nj+1)
  do i=1,n
    x=((i-d0)/(n+d1))**lam
    xx1(nj1-i+1)=sred(d0,d1,xx(nj+1),xx(nj),x)
  enddo
  xx1(1:nj)=xx(1:nj)
endif
nj=nj1
xx(1:nj+1)=xx1(1:nj+1)
end

subroutine init_meshval3
!значения в треугольниках для задачи пуассона или Гельмгольца для квадратной области
use pgmod
integer(4) i,ku,j
real(8) f_puasson,f_helmgolz,f_syst,f_v
real(8), allocatable :: g(:), x(:), y(:)
integer(4) ntr,pg_get_int
ntr=pg_get_int(4,1)
ku=gsarea%type_eq(1)
allocate(g(ntr),x(ntr),y(ntr))
call pg_get_array_real(4,1,x,ntr)
call pg_get_array_real(4,2,y,ntr)
if ((ku.eq.4).or.(ku.eq.6).or.(ku.eq.17).or.(ku.eq.18).or.(ku.eq.21)) then
  do i=1,ntr 
	  g(i)=f_puasson(x(i),y(i))
  enddo
  call pg_init_area_gu(g,1)
endif
if ((ku.eq.5).or.(ku.eq.6)) then
  do i=1,ntr 
    g(i)=f_helmgolz(x(i),y(i))
  enddo
  call pg_init_area_gu(g,2)
endif
if ((ku.eq.8).or.(ku.eq.11)) then
  do i=1,ntr 
    g(i)=f_syst(x(i),y(i))
  enddo
  call pg_init_area_gu(g,3)
endif
if ((ku.eq.16).or.(ku.eq.17)) then
  do i=1,ntr 
    g(i)=f_v(x(i),y(i),1)
  enddo
  call pg_init_area_gu(g,4)
  do i=1,ntr 
    g(i)=f_v(x(i),y(i),2)
  enddo
  call pg_init_area_gu(g,5)
endif
if (ku.eq.18) then
  do i=1,ntr 
    g(i)=f_v(x(i),y(i),2)
  enddo
  call pg_init_area_gu(g,4)
endif
if ((ku.eq.8).or.(ku.eq.11)) then
  ku=gsarea%type_eq(2)
  if (ku.eq.7) then
    do i=1,ntr 
      g(i)=f_puasson(x(i),y(i))
    enddo
    call pg_init_area_gu(g,1)
  endif
  if (ku.eq.14) then
    do i=1,ntr 
      g(i)=f_helmgolz(x(i),y(i))
    enddo
    call pg_init_area_gu(g,2)
  endif
endif
if (ku==22.or.ku==26) then
  do j=1,gsarea%n_eq_var
    if (gsarea%eq_var(j)) then
      do i=1,ntr 
        g(i)=f_v(x(i),y(i),j)
      enddo
      call pg_init_area_gu(g,j)
    endif
  enddo
endif
deallocate(g,x,y)
end

subroutine init_meshval_bri
!значения в треугольниках для уранения Бринкмана
use pgmod2
integer(4) ku2
ku2=gsarea%type_eq(2)
if (ku2==14) call pg_init_area_gu_const(-s_darci2,2)
call pg_init_area_gu_const(d1,3) !уравнения составляем для eta=-omega, поэтому тут знака "-" нет
end

subroutine init_mesh4
!инициализация узлов сетки в полукруге
use pgmod2
real(8) ds
integer(4) nx,i,nx1
real(8), allocatable :: rr(:)
ds=pi/nj_
nx=d1/ds !число отрезков по лучу единичного радиуса
ds=d1/nx
nx1=nx+1
allocate(rr(nx1))
do i=1,nx1
  rr(i)=(i-d1)*ds
enddo
call init_mesh4r(rr,nx1,nj_)
deallocate(rr)
end

subroutine init_mesh4_1
!инициализация узлов сетки в полукруге
!тестирование локального сгущения
use pgmod2
real(8) ds
integer(4) nx,i,nx1,j,k
real(8), allocatable :: rr(:)
real(8) dr,rr1
dr=1.0d-2
ds=pi/nj_
nx=d1/ds !число отрезков по лучу единичного радиуса
j=nx/2
ds=d1/nx
nx1=nx+1
allocate(rr(nx1+1))
k=0
do i=1,nx1
  rr1=(i-d1)*ds
  if (i==j) then
    rr(i)=rr1-dr
	rr(i+1)=rr1+dr
	k=1
	cycle
  endif
  rr(i+k)=rr1
enddo
call init_mesh4r(rr,nx1+1,nj_)
deallocate(rr)
end

subroutine init_mesh4_2(need_cylinder,updatetosquare)
!инициализация узлов сетки во всей ячейке как в единой области
use pgmod2
real(8) ds
integer(4) nx,i,nx1,j,k,i1,pg_get_int,n_add_r,n_add_g
logical need_cylinder,updatetosquare
real(8), allocatable :: rr(:), gg(:),rr2(:)
real(8) x(nmax),y(nmax)
call pg_bind_bound(1)
call pg_get_array_real(1,1,x,nmax)
call pg_get_array_real(1,2,y,nmax)
nx=0
if (need_cylinder) then
  ds=pi/nj_
  nx=d1/ds !число отрезков по лучу единичного радиуса
  ds=d1/nx
endif
j=nx
call pg_bind_boundline(1)
i=pg_get_int(2,1)
nx=nx+i
nx1=nx+1
n_add_r=0
n_add_g=0
if (need_thickening) then
  !n_add_r=n_thickening
  n_add_g=n_thickening*2
endif
allocate(rr(nx1+n_add_r))
if (need_cylinder) then 
  do i=1,j
    rr(i)=(i-d1)*ds
  enddo
endif
do i=j+1,nx1
  rr(i)=x(i-j)
enddo
call pg_bind_boundline(2)
k=pg_get_int(2,2)
if (updatetosquare) call pg_bind_boundline(4)
j=pg_get_int(2,3)-k+1
allocate(gg(j+1+n_add_g),rr2(j+1))
do i=1,j+1
  i1=i-1+k
  gg(i)=datan2(y(i1),x(i1))
  if (gg(i)<d0) gg(i)=gg(i)+pi2
  rr2(i)=dsqrt(x(i1)**2+y(i1)**2)
enddo
if (npe_==4.and.need_thickening) then
  call add_thickening(gg,j,n_thickening,0)
  call add_thickening(gg,j,n_thickening,1)
  !call add_thickening(rr,nx,n_thickening,0)
  !nx1=nx+1
endif
if (npe_==3) then
  call ga_init_mesh_rcell(rr,nx1,j+1,gg,.true.)
else
  call ga_init_mesh_rcell_quads(rr,gg,nx1,j+1)
endif
if (updatetosquare) then
  call ga_mesh_square(rr(1),rr(nx1),j,rr2,gg)
endif
deallocate(rr,gg,rr2) 
end

subroutine init_mesh4r(rr,nx1,ng_last)
!инициализация узлов сетки в полукруге для произвольного разбиения по радиусу
!равномерное разбиение на последнем радиусе
use gen_mod
integer(4) nx1,i,ng_last
real(8) dg
real(8), allocatable :: gg(:)
real(8) rr(nx1)
dg=pi/ng_last
allocate(gg(ng_last+1))
do i=1,ng_last+1
  gg(i)=(i-d1)*dg
enddo
call ga_init_mesh_rcell(rr,nx1,ng_last+1,gg,.true.)
deallocate(gg)
end

subroutine init_meshval_per
!значения в треугольниках для задачи осаждения взвеси на непроницаемом цилиндре (уравнение переноса)
use pgmod2
integer(4) i,ntr,pg_get_int
real(8), allocatable :: xx(:),yy(:),g(:,:)
real(8) pg_get_fun_1
ntr=pg_get_int(4,1)
allocate(xx(ntr),yy(ntr),g(ntr,2))
call pg_get_array_real(4,1,xx,ntr)
call pg_get_array_real(4,2,yy,ntr)
call pg_bind_problem(1)
call pg_bind_domain(1)
do i=1,ntr
  g(i,1)=pg_get_fun_1(xx(i),yy(i),2)*pe
  g(i,2)=-pg_get_fun_1(xx(i),yy(i),1)*pe
enddo
call pg_bind_problem(2)
call pg_bind_domain(1)
call pg_init_area_gu(g(:,1),4)
call pg_init_area_gu(g(:,2),5)
deallocate(xx,yy,g)
end

subroutine test_mesh_ref
!тестирование создания треугольной сетки с переменным шагом
use pgmod
call pg_allocate_problems(1)
call pg_bind_problem(1)
call pg_allocate_domains(1)
call pg_bind_domain(1)
call pg_allocate_area(2)
call pg_bind_areapart(1)
call init_mesh3(d1,d1,1,1)
call pg_bind_areapart(2)
!call pg_init_areapart_geom_ref(1,1,1,d1,d0,d1,d0,d1,d0) !смещение
call pg_init_areapart_geom_ref(1,1,1,-d1,d0,d0,d0,d1,d0) !отражение относительно оси Y
call pg_areageom_postprocessor
call ga_drw_trmesh(1)
end

subroutine per_mesh_aninodal
use pgmod2
use ani_mod
type(Tani_mesh) am
integer(4) i,nvfix
real(8) pg_get_fun_xy2,Quality
real(8), allocatable :: u(:),Metric(:,:)
!integer(4), allocatable :: fixedV(:)
integer(4) fixedV(1)
integer(4) in_bound
call pgani_meshfromarea(am,gsarea%a%ntr*10)
allocate(u(am%nv),Metric(3,am%nv))
!allocate(fixedV(am%nv))
nvfix=0
do i=1,am%nv
  u(i)=pg_get_fun_xy2(am%vrt(1,i),am%vrt(2,i),1,d0,d0,0,in_bound)
!  if (in_bound) then
!    nvfix=nvfix+1
!    fixedV(nvfix)=i
!  endif
enddo
!call ani_lmrNodal2MetricZZ(am,u,Metric)
call ani_lmrNodal2MetricVAR(am,u,Metric)  !эта процедура улавливает погран слой
Quality=0.5d0
call ani_mbaNodal2(am,Metric,nvfix,fixedV,am%nt,Quality)
call deallocate_area(1)
call pg_allocate_area(1)
call pg_bind_areapart(1)
call pgani_addareapart(am,1)
call pg_areageom_postprocessor
deallocate(u,Metric)
!deallocate(fixedV)
end
subroutine main_konform
!обтекание цилиндра течением стокса с использованием конформного отображени€
use pgmod2
integer(4) ku1,ku2
nmain=16
ku1=8  !8,11
ku2=2 

call pg_allocate_problems(1)
call pg_bind_problem(1)
call init_kuv(h)
call pg_allocate_domains(1)
call pg_bind_domain(1)
call pg_set_domain_equation_syst(ku1,ku2)
call pg_allocate_bounds(1)
call pg_bind_bound(1)
call init_geom3_konf3
call pg_geom_postprocessor
call init_gu2(.false.,.false.,.true.,.true.,.false.)
call pg_allocate_area(1)
call pg_bind_areapart(1)
!call init_mesh3(dlog(h),pi,njr_,nj_)
call init_mesh3(dlog(h),pi,nj_,njr_) !nj и ny перепутаны местами специально дл€ сгущени€ по оси x
call pg_areageom_postprocessor
call pg_allocate_area_gu
call init_meshval_konf
call pg_get_matrix
call pg_solve
call pg_get_psioml
call drw_solv2circ(nj_)
end


subroutine init_geom3_konf3
!инициализаци€ геометрии дл€ пр€моугольной области (пересчет из координат €чейки)
use pgmod2
integer(4) i,j,n,pg_get_int
real(8) g,x(nmax),y(nmax)
call init_geom2(.false.)
do j=1,4
  call pg_bind_boundline(j)
  n=pg_get_int(2,1)
  call pg_get_array_real(2,1,x,n+1)
  call pg_get_array_real(2,2,y,n+1)
  do i=1,n+1
    g=datan2(y(i),x(i))
    if (g<d0) g=g+pi2
    x(i)=d5*dlog(x(i)**2+y(i)**2)
    y(i)=g
  enddo
  call pg_init_boundline_geom(j,n,x,y)
enddo
end

subroutine init_meshval_konf
!значени€ в треугольниках дл€ задачи пуассона или √ельмгольца
use pgmod
real(8), allocatable :: g(:), x(:)
integer(4) ntr,pg_get_int
ntr=pg_get_int(4,1)
allocate(g(ntr),x(ntr))
call pg_get_array_real(4,1,x,ntr)
g=dexp(d2*x)   !!!помен€л тут знак, чтобы втора€ функци€ =-om
call pg_init_area_gu(g,3)
deallocate(g,x)
end

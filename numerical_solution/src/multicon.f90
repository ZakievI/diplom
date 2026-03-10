subroutine main_mc1
!течение стокса в €чейке кувабара через двусв€зную область
use pgmod2
integer(4) ku
nmain=17
call init_kuv(h)
ku=3
call pg_allocate_problems(1)
call pg_bind_problem(1)
call pg_allocate_domains(1)
call pg_bind_domain(1)
call pg_set_domain_equation(ku)
call pg_allocate_bounds(2)
call pg_bind_bound(1)
call ga_init_geom_circle(d0,d0,h,.false.,120)
call pg_bind_bound(2)
call ga_init_geom_circle(d0,d0,d1,.true.,60)
call pg_geom_postprocessor
call init_gu_mc1
call pg_get_matrix
call pg_solve
call pg_get_psioml
call drw_solv2circ(nj_)
end

subroutine main_mc2
real(8) ex
ex=0.5d0
call main_mc2_(ex,.true.)
end

subroutine main_mc3
use pgmod
integer(4) i
real(8) ex
OPEN (1,FILE='mc2e.dat')
do i=0,9
  ex=i*0.1d0
  call main_mc2_(ex,.false.)
  write(*,"('e=',f4.1,' psi=',f7.4)")ex, gs%a(1)%bnd(1)%psiom(1,1)
enddo
close(1)
end

subroutine main_mc2_(ex,need_drw)
!течение стокса в жидкой муфте
use pgmod
real(8) r_1,r_2,ex,om,c,e
integer(4) nnn,ku
logical need_drw
nnn=80
r_1=d5
r_2=d1
c=r_2-r_1
e=ex*c
om=d1

ku=3
call pg_allocate_problems(1)
call pg_bind_problem(1)
call pg_allocate_domains(1)
call pg_bind_domain(1)
call pg_set_domain_equation(ku)
call pg_allocate_bounds(2)
call pg_bind_bound(1)
call ga_init_geom_circle(-e,d0,r_1,.true.,nnn)
call pg_bind_bound(2)
call ga_init_geom_circle(d0,d0,r_2,.false.,nnn)
call pg_geom_postprocessor
call init_gu_mc2(om*r_1)
call pg_get_matrix
!call closing_matrix_mc2
call pg_closing_psiconst(1)
call pg_solve
if (need_drw) then
  call pg_get_psioml
  call drw_solv_mc2(r_1,r_2,e,om)
endif
end

subroutine init_gu_mc1
!инициализаци€ граничных условий дл€ задачи течени€ стокса в €чейке кувабара через двусв€зную область
use pgmod2
real(8) yc(nmax)
!внеш круг
call pg_bind_bound(1)
call pg_allocate_bound_gu
call pg_bind_boundline(1)
call pg_get_array_real(2,4,yc,nmax)
call pg_init_boundline_gu_val(1,1,yc)
call pg_init_boundline_gu_val_const(2,3,d0)
!внутр круг
call pg_bind_bound(2)
call pg_allocate_bound_gu
call pg_bind_boundline(1)
call pg_init_boundline_gu_val_const(1,1,d0)
call pg_init_boundline_gu_val_const(2,2,d0)
end

subroutine init_gu_mc2(omr)
!инициализаци€ граничных условий дл€ задачи течени€ стокса в €чейке кувабара через двусв€зную область
use pgmod2
real(8) omr
call pg_allocate_constvalind(1)
!внеш круг
call pg_bind_bound(1)
call pg_allocate_bound_gu
call pg_bind_boundline(1)
call pg_init_boundline_gu_constval(1,1,1)
call pg_init_boundline_gu_val_const(2,2,omr)
!внутр круг
call pg_bind_bound(2)
call pg_allocate_bound_gu
call pg_bind_boundline(1)
call pg_init_boundline_gu_val_const(1,1,d0)
call pg_init_boundline_gu_val_const(2,2,d0)
end

subroutine closing_matrix_mc2
!замыкаем —Ћј” задачи дл€ жидкой муфты
use pgmod
integer(4) i,np1,pg_get_int,nx
integer(4), allocatable :: b1psiind4(:)
real(8), allocatable :: eq(:)
real(8) b
type(col_info), target :: ci
integer(4), allocatable :: buff_c(:)
logical, allocatable :: ind(:)
allocate(ind(gs%m%nx),buff_c(gs%m%nx))
ind=.false.
call pg_bind_bound(1)
np1=pg_get_int(1,1)
allocate(b1psiind4(np1))
nx=pg_get_int(3,1)
allocate(eq(nx))
call pg_get_array_int2(1, 1, 4, b1psiind4, np1)
eq=d0
do i=1,np1
  eq(b1psiind4(i))=gsbnd%l(i)
enddo
do i=1,gsbnd%nline
  call init_ind_psiom(gsbnd,gsbnd%line(i),ind,4)
enddo
call init_col_info(gsarea,ci,ind,buff_c,gs%m%nx,.false.)
gsarea%m%ci=>ci
b=d0
call pg_add_closing_eq_to_area(gsarea%i,eq,b)
deallocate(b1psiind4,eq)
deallocate(ind,buff_c)
call deallocate_col_info(ci)
end

subroutine drw_solv_mc2(r_1,r_2,e,om)
!вывод в файл решени€ бигармонического уравнени€ (дл€ решени€ в многосв€зной области)
use pgmod2
integer(4) i,j,nng,nnr
real(8) g,xx,yy,f,pg_get_fun,fom,pg_get_fun2,pg_f_psioml,s,pg_get_fun_1
real(8) om,r_1,r_2,e,rr,rr1
real(8) fu,fv,fv2
complex(8) uv 
character(75) str

nnr=30
nng=60
OPEN (1,FILE='mc2.dat')
write(1,*) 'TITLE = "velosity"'
write(1,*) 'VARIABLES = "X", "Y", "psi", "om", "u", "v", "v2"'
write(1,"('ZONE T=""area_e"", I=', i0, ', J=', i0, ', F=POINT')") nnr+1,nng+1
call init_gs_time
do i=1,nng+1
  write(str,"('Write file area_e. ', I0, ' (from ', I0, ')')") i, nng+1
  call write_by_time(str)
  g=(i-d1)/nng*pi2
  rr=e*dcos(g)
  rr=rr+dsqrt(rr**2+r_2**2-e**2)
  do j=1,nnr+1
    rr1=r_1+(j-d1)/nnr*(rr-r_1)
    xx=-e+rr1*dcos(g)
 	yy=rr1*dsin(g)
    if (j==1) then
	  s=(pi2-g)*r_1
	  f=gsarea%bnd(1)%psiom(1,1)
	  fom=pg_f_psioml(1,s,3,0)
	  uv=ii*cdexp(ii*g)*r_1*om
	  fu=dreal(uv)
	  fv=dimag(uv)
	  fv2=cdabs(uv)
    elseif (j==nnr+1) then
	  f=d0
	  fom=d0
	  fu=d0
	  fv=d0
	  fv2=d0
	else
	  f=pg_get_fun(xx,yy)
	  fom=-pg_get_fun2(xx,yy)
	  fu=pg_get_fun_1(xx,yy,2)
	  fv=-pg_get_fun_1(xx,yy,1)
	  fv2=dsqrt(fu**2+fv**2)
	endif
	WRITE(1,"(F9.5, ' ', F9.5, ' ', F9.5, ' ', F11.5, ' ', F9.5, ' ', F9.5, ' ', F9.5, ' ', F9.5, ' ', F9.5, ' ', F9.5)") xx,yy,f,fom,fu,fv,fv2
  enddo
enddo
close(1)
end

subroutine main_mc4
!тестирование течени€ в квадратной €чейке с круговой частией типа bl2
use pgmod2
integer(4) ku,i,j,k
integer(4) mode ! 0 - одна частица, 1 - пориста€ среда
mode=1
nj_=10
mcell%x0=d0
mcell%y0=d0
mcell%dx=d1
mcell%dy=d1
mcell%nx=1
if (mode==1) mcell%nx=nj_
mcell%ny=mcell%nx
mcell%m=0.999d0
call save_mcell("mcell.dat",mcell)
ku=3
call dns_init_mcell1(mcell,pi)
call dns_init_mcell2(mcell)
call dns_init_mcell3(mcell,pi2,nj_)
call pg_allocate_problems(1)
call pg_bind_problem(1)
call pg_allocate_domains(1)
call pg_bind_domain(1)
call pg_set_domain_equation(ku)
call pg_allocate_bounds(1+mcell%nc)
call pg_bind_bound(1)
call init_geom3
k=1
do i=1,mcell%nx
  do j=1,mcell%ny
    k=k+1
    call pg_bind_bound(k)
    call ga_init_geom_cylider_particle(mcell%x(i),mcell%y(j),mcell%r,0)
  enddo
enddo
call pg_geom_postprocessor
call init_gu_mc4
call pg_get_matrix
call pg_solve
call pg_get_psioml
call drw_solv_mc4(mode)
end

subroutine init_gu_mc4
use pgmod2
call pg_bind_bound(1)
call pg_allocate_bound_gu
!низ
call pg_bind_boundline(1)
call pg_init_boundline_gu_val_const(1,1,d0)
call pg_init_boundline_gu_val_const(2,3,d0)
!право
call pg_bind_boundline(2)
call pg_init_boundline_gu_val_const(1,2,d0)
call pg_init_boundline_gu_val_const(2,4,d0)
!верх
call pg_bind_boundline(3)
call pg_init_boundline_gu_val_const(1,1,d1)
call pg_init_boundline_gu_val_const(2,3,d0)
!лево
call pg_bind_boundline(4)
call pg_init_boundline_gu_val_const(1,2,d0)
call pg_init_boundline_gu_val_const(2,4,d0)
!частица
call pg_init_bound_gu_cylider_particle_all
end

subroutine drw_solv_mc4(mode)
use pgmod2
real(8) rr(nmax),bndg(nmax),bndrv(nmax)
integer(4) nnr,nng,nng_pi,mode
external write_point_mc4
call gs_print('write 2.dat')
call init_gs_time
OPEN (1,FILE='2.dat')
write(1,*) 'TITLE = "velosity"'
write(1,*) 'VARIABLES = "X", "Y", "psi", "om", "u", "v"'
if (mode==0) then
  nng_pi=320
  call dns_prepare_drw_one_cell_g(mcell,nng_pi,nng,bndg,bndrv,nmax)
  call dns_prepare_drw_one_cell_r(mcell,pi/nng_pi,24,32,nnr,rr,nmax)
  call dns_drw_one_cell(1,mcell,1,1,nng,nnr,bndg,bndrv,rr,' ',.true.)
else
  call dns_write_mcell_in2(1,mcell,1,mcell%nx,1,mcell%ny,0,write_point_mc4)
endif
close(1)
call print_total_time
end

subroutine write_point_mc4(nfile,x,y,i1,j1,mode)
use pgmod2
integer(4) nfile,i1,j1,mode
real(8) x,y,psi1,om,u,v,ff(2),get_fun_xy_solve_bi
type(TBound_info) bi
call get_fun_xy_prepare_bi(x,y,3,bi)
psi1=get_fun_xy_solve_bi(x,y,1,d0,d0,bi)
om=-get_fun_xy_solve_bi(x,y,3,d0,d0,bi)
!u=get_fun_xy_solve_bi(x,y,2,d0,d1,bi)
!v=-get_fun_xy_solve_bi(x,y,2,d1,d0,bi)
call get_fun_xy_solve_gradient_bi(x,y,5,d0,d0,bi,ff)
u=ff(2)
v=-ff(1)
write(nfile,"(6(E15.6))") x,y,psi1,om,u,v
return
i1=i1
j1=j1
mode=mode
end
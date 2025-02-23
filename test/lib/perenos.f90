subroutine main_per1
!тестирование уравнения переноса в квадрате
use pgmod2
integer(4) ku
ku=17 !16,17
nmain=8
nfun=28
	   !ku=16
	   !26 - y*y+2*x      (vx=d1;vy=d0)
	   !27 - y*dexp(x)    (vx=d1;vy=d0)

	   !ku=17
       !28 - xx**2+yy**2  (vx=xx;vy=yy**2)
	   !29 - xx+yy        (vx=xx;vy=yy)
     !47 - x**2+y**2   (vx={0 : x<0.5; x : x>0.5} vy={y : x<0.5; 0 : x>0.5})

call pg_allocate_problems(1)
call pg_bind_problem(1)
call pg_allocate_domains(1)
call pg_bind_domain(1)
call pg_set_domain_equation(ku)
call pg_allocate_bounds(1)
call pg_bind_bound(1)
call init_geom3
call pg_geom_postprocessor
call init_gu3
call pg_allocate_area(1)
call pg_bind_areapart(1)
call init_mesh3(d1,d1,nj_,nj_)
call pg_areageom_postprocessor
call pg_allocate_area_gu
call init_meshval3
call pg_get_matrix
call pg_solve
call pg_get_psioml
call drw_solv3
end

subroutine main_per2
!осаждение взвеси в квадратном канале
use pgmod2
integer(4) ku
ku=16 
nmain=9
nfun=31
	   !30 - осаждение взвеси в канале    (vx=d1;vy=d0)
	   !31 - осаждение взвеси в канале    (vx=течение Пуазеля;vy=d0)
Pe=1.0d0
call pg_allocate_problems(1)
call pg_bind_problem(1)
call pg_allocate_domains(1)
call pg_bind_domain(1)
call pg_set_domain_equation(ku)
call pg_allocate_bounds(1)
call pg_bind_bound(1)
call init_geom3
call pg_geom_postprocessor
call init_gu3_per
call pg_allocate_area(1)
call pg_bind_areapart(1)
call init_mesh3(d1,d1,nj_,nj_)
call pg_areageom_postprocessor
call pg_allocate_area_gu
call init_meshval3
call pg_get_matrix
call pg_solve
call pg_get_psioml
call drw_solv3
end

subroutine main_per3
!осаждение взвеси при обтекании непроницаемого цилиндра
use pgmod2
real(8) par(5)
integer(4) nng

call pg_allocate_problems(2)

!задача обтекания
call pg_bind_problem(1)
call main2_(.false.)

!задача взвеси
Pe=10.0d0
nng=nj_
call pg_bind_problem(2)
par(1)=d1
call main_per3_2(par,nng)
call main_per3_3
call pg_get_psioml
!call save_dcdn
call drw_solv2mode(0,par,nng,.true.)
if (npe_==3) call main_per3_4
end

subroutine main_per3_2(par,nng)
!осаждение взвеси при обтекании непроницаемого цилиндра
!расчетная часть для взвеси
use pgmod2
integer(4) ku,nng
real(8) par(5)
ku=16 
call pg_allocate_domains(1)
call pg_bind_domain(1)
call pg_set_domain_equation(ku)
call pg_allocate_bounds(1)
call pg_bind_bound(1)
call init_geom2mode(0,.true.,nng,par)
call pg_geom_postprocessor
call init_gu2_per
call pg_allocate_area(1)
call pg_bind_areapart(1)
call init_mesh4_2(.false.,.true.)
call pg_areageom_postprocessor
call pg_allocate_area_gu
call ga_drw_trmesh(1)
end

subroutine main_per3_3
!осаждение взвеси при обтекании непроницаемого цилиндра
!расчетная часть для взвеси
call init_meshval_per
call pg_get_matrix
call pg_solve
end

subroutine main_per3_4
!перерасчет треугольной сетки по полю концентрации
use pgmod2
call per_mesh_aninodal
call pg_allocate_area_gu
call ga_drw_trmesh(1)
call main_per3_3
call pg_get_psioml
drw_e_tr=.true.
call drw_solv2circ(nj_)
!call drw_solv2mode(0,par,nng,.true.)
end

subroutine save_dcdn
use pgmod2
real(8) q(gs%a(1)%bnd(1)%npanel,2)
q(:,2)=gs%a(1)%bnd(1)%psiom(:,2)
q(:,1)=gs%a(1)%bnd(1)%s(gs%a(1)%bnd(1)%npanel+1)-gs%a(1)%bnd(1)%sc
end

subroutine main_per5
!решение уравнения Дарси с переменной проницаемостью в круговой ячейке
use pgmod2
integer(4) ku
ku=16 !16,17
nmain=10
nfun=32
s_darci=d2

call init_kuv_per_darci(s_darci)
call pg_allocate_problems(1)
call pg_bind_problem(1)
call pg_allocate_domains(1)
call pg_bind_domain(1)
call pg_set_domain_equation(ku)
call pg_allocate_bounds(1)
call pg_bind_bound(1)
call init_geom4_1                          
call pg_geom_postprocessor
call init_gu4_4
call pg_allocate_area(1)
call pg_bind_areapart(1)
call init_mesh4
call pg_areageom_postprocessor
call pg_allocate_area_gu
call init_meshval3
call pg_get_matrix
call pg_solve
call pg_get_psioml
call drw_solv_per_darci2
end

subroutine test_per_mesh
!тестирование создания треугольной сетки с переменным шагом
use pgmod2
integer(4) ku
ku=16 
call pg_allocate_problems(1)
call pg_bind_problem(1)
call pg_allocate_domains(1)
call pg_bind_domain(1)
call pg_set_domain_equation(ku)
call pg_allocate_bounds(1)
call pg_bind_bound(1)
call init_geom4_1                          
call pg_geom_postprocessor
call init_mesh4_1
call ga_drw_trmesh(1)
end
subroutine main_oss1
!тестирование осесимметричных уравнений в квадрате как частный случай уравнени€ переноса
use pgmod2
integer(4) ku
ku=20 !17,18,19,20
k_oss=d1
square_shift_dy=0.0d0
nmain=20
nfun=39
          !17,18
          !37 - y*y+x*x      (vx=0,vy=x*y)
		  !38 - x            (vx=0,vy=-1/y)

		  !19,20
		  !38 - x, k_oss=d1         уравнение Ћапласа
		  !39 - y**2/2, k_oss=-d1   уравнение дл€ функции тока
if (nfun==39) k_oss=-d1
need_thickening=.true.
call pg_allocate_problems(1)
call pg_bind_problem(1)
call pg_allocate_domains(1)
call pg_bind_domain(1)
call pg_set_domain_equation(ku)
if (ku==19.or.ku==20) call pg_set_areaconst(2,k_oss) 
call pg_allocate_bounds(1)
call pg_bind_bound(1)
call init_geom3_dy(square_shift_dy)
call pg_geom_postprocessor
call init_gu3
call pg_allocate_area(1)
call pg_bind_areapart(1)
call init_mesh3_1(d0,d1,square_shift_dy,square_shift_dy+d1,nj_,nj_)
call pg_areageom_postprocessor
call ga_drw_trmesh(1)
call pg_allocate_area_gu
if (ku==17.or.ku==18) then
  call init_meshval3
else
  call pg_init_area_gu_oss
endif
call pg_get_matrix
call pg_solve
call pg_get_psioml
call drw_solv3
end

subroutine main_oss2
!тестирование осесимметричных уравнений. обтекание сферы в €чейке кувабара
use pgmod2
integer(4) ku
ku=19 !19,20
nmain=21
k_oss=-d1 !дл€ функции тока
need_thickening=.true.
call pg_allocate_problems(1)
call pg_bind_problem(1)
call pg_allocate_domains(1)
call pg_bind_domain(1)
call pg_set_domain_equation(ku)
call pg_set_areaconst(2,k_oss) 
call pg_allocate_bounds(1)
call pg_bind_bound(1)
call init_geom2(.false.)
call pg_geom_postprocessor
call init_gu1
call pg_allocate_area(1)
call pg_bind_areapart(1)
call init_mesh4_2(.false.,.false.)
call pg_areageom_postprocessor
call ga_drw_trmesh(1)
call pg_allocate_area_gu
call pg_init_area_gu_oss
call pg_get_matrix
call pg_solve
call pg_get_psioml
  drw_e_tr=.false.
call drw_solv1
end


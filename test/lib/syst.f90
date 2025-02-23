subroutine main_s1
!система вида \Delta\psi=v(x,y)\omega  и уравнение на \omega
use pgmod2
integer(4) ku1,ku2
nmain=5 
nfun=21 
        !ku2=2
		!18 - get_fund=xx
        !19 - get_fund=xx**3

		!ku2=12
		!20 - get_fund=yy*dsin(k_helm*xx)

		!ku2=13
		!21 - get_fund=yy*dexp(k_helm*xx)
ku1=11
ku2=13 !2,12,13
k_helm=d2

call pg_allocate_problems(1)
call pg_bind_problem(1)
call pg_allocate_domains(1)
call pg_bind_domain(1)
call pg_set_domain_equation_syst(ku1,ku2)
if (ku2==12.or.ku2==13) call pg_set_areaconst(1,k_helm) 
call pg_allocate_bounds(1)
call pg_bind_bound(1)
call init_geom3
call pg_geom_postprocessor
call init_gu3_2
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

subroutine main_s2
!система вида \Delta\psi=v(x,y)\omega и уравнение на \omega (\omega в области неизвестна)
use pgmod2
integer(4) ku1,ku2
nmain=6 
nfun=24
        !ku2=2
        !19 - get_fund=xx**3  
		!23 - get_fund=xx**4  get_fund2=6*xx

		!ku2=7
		!22 - get_fund=xx**4  get_fund2=12*xx**2

		!ku2=12
		!20 - get_fund=yy*dsin(k_helm*xx)

		!ku2=13
		!21 - get_fund=yy*dexp(k_helm*xx)

		!ku2=14
		!24 - get_fund=dsin(xx*yy)  get_fund2=dsin(xx*yy)
		!25 - get_fund=yy*dexp(k_helm*xx) get_fund=k_helm**2*yy*dexp(k_helm*xx)

ku1=8
ku2=14 !2,7,12,13,14
k_helm=d2

call pg_allocate_problems(1)
call pg_bind_problem(1)
call pg_allocate_domains(1)
call pg_bind_domain(1)
call pg_set_domain_equation_syst(ku1,ku2)
if (ku2==12.or.ku2==13) call pg_set_areaconst(1,k_helm) 
call pg_allocate_bounds(1)
call pg_bind_bound(1)
call init_geom3
call pg_geom_postprocessor
call init_gu3_2
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
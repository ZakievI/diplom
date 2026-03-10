subroutine main_h
use pgmod2
integer(4) ku
ku=10 !9,10
nmain=15
nfun=9
       !ku=9
       !4 - get_funb=dsin(xx)
	   !6 - get_funb=dsin(k_helm*xx)
       !7 - get_funb=yy*dsin(xx)
       !8 - get_funb=yy*dsin(k_helm*xx)

	   !ku=10
	   !5 - get_funb=dexp(xx)
	   !9 - get_funb=yy*dexp(k_helm*xx) 

k_helm=d2
if (nfun==4.or.nfun==5.or.nfun==7) k_helm=d1

call pg_allocate_problems(1)
call pg_bind_problem(1)
call pg_allocate_domains(1)
call pg_bind_domain(1)
call pg_set_domain_equation(ku)
call pg_set_areaconst(1,k_helm) 
call pg_allocate_bounds(1)
call pg_bind_bound(1)
call init_geom3
call pg_geom_postprocessor
call init_gu3
call pg_get_matrix
call pg_solve
call pg_get_psioml
call drw_solv3
end

subroutine testBesselF
use pgmod2
integer(4) nf,i
real(8) f,calcintn,f1,f2,f3,e1,e2
complex(8) z
nf=18
do i=1,6
  call pg_allocate_problems(1)
  call pg_bind_problem(1)
  call pg_allocate_domains(1)
  call pg_bind_domain(1)
  call pg_set_domain_equation(5)
  k_helm=d2
  call pg_set_areaconst(1,k_helm) 
  call pg_allocate_bounds(1)
  call pg_bind_bound(1)
  call pg_allocate_boundlines(1)
  select case (i)
    case(1)
      call pg_init_boundline_geomline(1,1,d0,d0,d1,d0)
    case(2)
      call pg_init_boundline_geomline(1,1,d1,d0,d0,d0)
    case(3)
      call pg_init_boundline_geomline(1,1,d0,d0,d0,d1)
    case(4)
      call pg_init_boundline_geomline(1,1,d0,d1,d0,d0)
    case(5)
      call pg_init_boundline_geomline(1,1,d0,d0,d1,d1)
    case(6)
      call pg_init_boundline_geomline(1,1,d0,d0,d2,d1)
    end select
  call pg_geom_postprocessor
  call bind_AreaConst(1)
  gs_test_point_collenear_for_BesselG=.true.
  z=dcmplx(d2,d0)
  !z=dcmplx(d2,d2)
  f=calcintn(1,dreal(z),dimag(z),nf,1,1)
  f1=calcintn(1,dreal(z),dimag(z),nf+1,1,1)
  call calcintn_dual(1,dreal(z),dimag(z),nf,1,1,f2,f3)
  e1=f-f2
  e2=f1-f3
  call pg_deallocate_mem 
enddo
end
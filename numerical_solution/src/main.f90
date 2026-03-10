program main
use pgmod2

implicit none

call pg_start
call pg_start_write_log_file

!gs_test_areaval_eq_0=.true.
gs_use_cash=.true.
gs_cash_presolve=.true.
!gs_use_dual_integral_solve=.true.
!gs_test_point_collenear_for_BesselG=.true.
!gs_use_parallel_get_fun=0

!gs_write_matrix=.true.
only_analit=.false.

npe_=4
h=5.0d0
hv=h !размер внешней области для выгрузки в файл
h2=d1 !для квадрата - сторона квадрата, для круга - радиус
nj_=61
square_shift_dy=d0
drw_e_tr=.false.
need_thickening=.false.
n_thickening=5
slip_k=1.0d0
use_subdomain=.false.
solve_only_psi=.false.
nmain=0
nfun=0

!call testBesselF
!call test_qr
!call testDistanceToCircleArc
!call test_fast_cash_solve
!call test_cash_2_problem
!call testselectif
 
!call main2
!call main_p7
call main_h
!call main_s2
!call main_konform
!call main_mc4
!call main_kuvbri5
!call main_per1
!call main_oss1

!call test_per_mesh
!call test_mesh_ref

!call test_cp_cint2
!call test_intd

call deallocate_mem
end

subroutine deallocate_mem
use pgmod2
call pg_finish
if (allocated(parr)) deallocate(parr)
if (allocated(pknow)) deallocate(pknow)
end

subroutine main1
!потенциальное течение в ячейке кувабара
use pgmod2
integer(4) ku
nmain=11
ku=1
call pg_allocate_problems(1)
call pg_bind_problem(1)
call pg_allocate_domains(1)
call pg_bind_domain(1)
call pg_set_domain_equation(ku)
call pg_allocate_bounds(1)
call pg_bind_bound(1)
call init_geom2(.false.)
call pg_geom_postprocessor
call init_gu1
call pg_get_matrix
call pg_solve
call pg_get_psioml
call drw_solv1
end

subroutine main2
use pgmod
!течение Стокса в ячейке Кувабара
call pg_allocate_problems(1)
call pg_bind_problem(1)
call main2_(.true.)
end

subroutine main2_(need_save)
!течение Стокса в ячейке Кувабара
use pgmod2
logical square_cell,need_save
real(8) par(5),m
integer(4) nng,mode,ku,periodic
nmain=12
ku=3
square_cell=.true.
periodic=4
if (.not.square_cell) periodic=0
mode=0 !0 - circle
       !1 - ellipse
	   !2 - virus
	   !3 - square
	   !4 - triangle (poly)

if (mode==0.and.(.not.square_cell)) then
  call init_kuv(h)
  nfun=33
endif

if (.false..and.mode==3.and.square_cell) then
  m=0.99d0
  h2=h*dsqrt(1-m)
  h=d5*dsqrt(pi/(d1-m))
endif

if (.false..and.mode==0.and.square_cell) then
  m=0.9d0
  h=d5*dsqrt(pi/(d1-m))
  slip_k=slip_k*d2*h
endif


!call pg_allocate_problems(1)
!call pg_bind_problem(1)
call pg_allocate_domains(1)

!call pg_set_matrix_sparse_count_closing(2)

call pg_bind_domain(1)
call pg_set_domain_equation(ku)
call pg_allocate_bounds(1)
call pg_bind_bound(1)
call init_geom2mode(mode,square_cell,nng,par)
call pg_geom_postprocessor
!call init_gu2(square_cell,periodic>0,.true.,.false.,.false.)
call init_gu2(square_cell,periodic>0,.true.,.false.,.false.)
if (periodic==3.or.periodic==4) call init_gu_per(periodic)
call pg_get_matrix
if (periodic==1.or.periodic==2) call closing_matrix2_per(periodic)
call pg_solve
if (need_save) then
  call pg_get_psioml
  call drw_solv2mode(mode,par,nng,square_cell)
  !call save_vx
  !call calc_pmap
  !call calc_dp
endif
end

subroutine main3
!уравнение лапласа в квадрате
use pgmod2
integer(4) ku

!real(8) pg_f_psiomlxy,fff

nmain=1
nfun=2   !1 - xx
		 !2 - d1-xx-yy+xx*yy
ku=1
call pg_allocate_problems(1)
call pg_bind_problem(1)
call pg_allocate_domains(1)
call pg_bind_domain(1)
call pg_set_domain_equation(ku)
call pg_allocate_bounds(1)
call pg_bind_bound(1)
call init_geom3
!call ga_init_geom_circle(d0,d0,d1,.false.,nj_)
!call init_geom4_2
call pg_geom_postprocessor
call init_gu3
call pg_get_matrix
call pg_solve
call pg_get_psioml

!fff=pg_f_psiomlxy(d1,d0,1,0)

call drw_solv3
end

subroutine main_subdomain(main_first_area,modemesh)
!решение тестовой задачи в квадрате с разбиением на подобласти
use pgmod2
real(8) dx,dy,dxx,dyy
integer(4) i,j,k
integer(4) modemesh !0 - не нужна сетка в области
                    !1 - сетка и значения в узлах сетки                    
external main_first_area
logical use_rotate
use_subdomain=.true.
use_rotate=.true.
nsx=2
nsy=2
dx=d1/nsx
dy=d1/nsy
call pg_allocate_problems(1)
call pg_bind_problem(1)
call pg_allocate_domains(nsx*nsy)
gs%m%matrix_store_type=2
k=0
do i=1,nsx
  dxx=(i-d1)*dx
  do j=1,nsy
    dyy=(j-d1)*dy
    k=k+1
    if (k==1) then
      call pg_bind_domain(1)
      call main_first_area   !задание типа уравнения
      call pg_allocate_bounds(1)
      call pg_bind_bound(1)
      call init_geom3_2
      call pg_geom_postprocessor
      if (modemesh==1) then
        call pg_allocate_area(1)
        call pg_bind_areapart(1)
        call init_mesh3(dx,dy,nj_,nj_2)
        call pg_areageom_postprocessor
        call pg_allocate_area_gu
        call init_meshval3
      endif
    else
      if (k==2.and.nsx==nsy.and.use_rotate) then
        call pg_copy_subdomain_dxdy_rotate(1,k,dxx,dyy,pi5,dx*d5,dy*d5)
      else
        call pg_copy_subdomain_dxdy_rotate(1,k,dxx,dyy,d0,d0,d0)
      endif
      if (gs_use_cash) call pg_set_cash_ref(1,1)
      if (modemesh==1) call init_meshval3
    endif
  enddo
enddo
if (gsarea%nu==1) then
  call init_gu3_3
else
  call init_gu3_4
endif
call pg_get_matrix
!gs_block_matrix_solver_write_all_iteration=.true.
gs_block_matrix_solver_lambda=0.7d0
gs_block_matrix_solver_maxIter=1000
!gs_block_matrix_solver_maxIter_convergence=500
!gs_block_matrix_solver_eps=1.0d-9
call pg_solve
call pg_get_psioml
call drw_solv3
end

subroutine main3_2_first
use pgmod2
integer(4) ku
nmain=1
nfun=2   !1 - xx
		 !2 - d1-xx-yy+xx*yy
ku=1
call pg_set_domain_equation(ku)
end

subroutine main3_2
!уравнение лапласа в квадрате
!тестирование сабдоменов
use pgmod2
external main3_2_first
call main_subdomain(main3_2_first,0)
end

subroutine main4
!обтекание пористого цилиндра в пористой среде
use pgmod2
integer(4) ku,nng,mode,nbndl1
logical square_cell,gu_closing
real(8) par(5)

nmain=13
ku=1
kkr=1.0d0/2.0d0   !k_e/k_i

!kkr=3.9128d0
!kkr=0.25557d0

!kkr=d1/2.4563d0

gu_closing=.false.

mode=0 !0 - circle
	     !3 - square
square_cell=.false.

call pg_allocate_problems(1)
call pg_bind_problem(1)
call pg_allocate_domains(2)

!call pg_set_matrix_sparse_count_closing(2)

!area_e
call pg_bind_domain(1)
call pg_set_domain_equation(ku)
call pg_allocate_bounds(1)
call pg_bind_bound(1)
call init_geom2mode(mode,square_cell,nng,par)
call pg_geom_postprocessor
call init_gu4_1(square_cell)
nbndl1=gsbnd%nline
!area_i
call pg_bind_domain(2)
call pg_set_domain_equation(ku)
call pg_allocate_bounds(1)
call pg_bind_bound(1)
call init_geom4
call pg_geom_postprocessor
call init_gu4_2
if (gu_closing) call pg_init_boundline_gu_darcy_darcy(1,1,nbndl1,2,1,1,kkr)
call pg_get_matrix
if (.not.gu_closing) call pg_closing_matrix_darcidarci(1,1,nbndl1,2,1,1,kkr)
call pg_solve
call pg_get_psioml
call drw_solv4mode(mode,par,nng,square_cell)
end

subroutine get_alf_bet_bivers(m,kap,mode)
use pgmod2
real(8) m,kap,k,sqrt_mu
integer(4) mode !0-darci, 1-bri
if (m==0.999999d0) then
  k=0.490957d0   
  mu_bri=0.103646d1
  alf_darci=0.863146d0   
  bet_darci=0.187556d1 !p(h)
  !bet_darci=0.937446d0 !Vs
  !alf_darci=0.223824d1 !с учетом дельта (граница пористой среды по границе частиц) 
  !bet_darci=0.219312d1 
  
  !alf_darci=0.863146d2   !увеличили в 100 раз

  
  !alf_darci=0.852623d0
  !bet_darci=alf_darci
  !alf_darci=0.788083d0  !by q
  !bet_darci=0.758766d0
  !alf_darci=2.0d0  !альфа испортили
  !bet_darci=alf_darci
elseif (m==0.99d0) then
  k=0.125282d0  
  mu_bri=0.140460d1
  alf_darci=0.664998d0   
  bet_darci=0.177299d1
  !alf_darci=0.651209d0
  !bet_darci=alf_darci
  !alf_darci=0.398678d1 !с учетом дельта (граница пористой среды по границе частиц)
  !bet_darci=0.175716d1
elseif (m==0.9d0) then
  k=0.402716d-1  
  mu_bri=0.333101d1
  alf_darci=0.525692d0
  bet_darci=0.196474d1
  !alf_darci=0.512121d0
  !bet_darci=alf_darci
  !alf_darci=0.333365d1 !с учетом дельта (граница пористой среды по границе частиц)
  !bet_darci=0.576556d0
elseif (m==0.7d0) then
  k=0.972004d-2
  mu_bri=0.138695d2
  alf_darci=0.409262d0
  bet_darci=0.256730d1
elseif (m==0.5d0) then
  k=0.187780d-2  
  mu_bri=0.728562d2
  alf_darci=0.296368d0   
  bet_darci=0.377267d1
  !alf_darci=0.290740d0
  !bet_darci=alf_darci
  !alf_darci=0.296607d0  !by q
  !bet_darci=0.373006d1
  !alf_darci=0.959447d0   !с учетом дельта (граница пористой среды по границе частиц)
  !bet_darci=0.221426d1
elseif (m==0.2d0) then
  k=0.102469d-3   !S2 - квадратики
  mu_bri=0.136301d4
  alf_darci=0.188940d+00   
  bet_darci=0.517908d+01
  !alf_darci=0.128074d+02  !с учетом дельта (граница пористой среды по границе частиц)
  !bet_darci=d0
endif
k=k/kap**2
if (mode==0) then
  s_darci=d1/dsqrt(k)
else
  !s_darci=d1/dsqrt(k*mu_bri)
  !a_bri=d1/(bet_darci-alf_darci)
  !b_bri=a_bri*bet_darci
  k_bri=k
  s_darci2=d1/(k_bri*mu_bri)
  s_darci=dsqrt(s_darci2)
  sqrt_mu=dsqrt(mu_bri)
  a_bri=d1/(alf_darci-bet_darci+sqrt_mu)
  b_bri=a_bri*(sqrt_mu-bet_darci)*sqrt_mu
  !a_bri=alf_darci
  !b_bri=bet_darci-sqrt_mu
endif
end

subroutine get_alf_bet_brink(m,kap)
use pgmod2
real(8) m,kap,k
if (m==0.999999d0) then
  k=0.490957d0   
  mu_bri=0.103646d1
  a_bri=0.147035d+03   
  b_bri=0.128157d+03 
elseif (m==0.9d0) then
  k=0.402716d-1
  mu_bri=0.333101d1
  a_bri=0.287960d+01
  b_bri=0.939105d+00
  !a_bri=-0.894d0
  !b_bri=-2.727d0
  !a_bri=2.59029d0  !мое через дарси
  !b_bri=0.660129d0
elseif (m==0.7d0) then
  k=0.972004d-2        
  mu_bri=0.138695d2
  a_bri=0.707972d0
  b_bri=-0.264409d1
elseif (m==0.5d0) then
  k=0.187780d-02      
  mu_bri=0.728562d+02   
  a_bri=0.219367d+00  
  b_bri=-0.797989d+01  
endif
!a_bri=d0
!b_bri=-dsqrt(mu_bri)
!a_bri=a_bri*d2
!b_bri=b_bri*d5
k_bri=k/kap**2
s_darci2=d1/(k_bri*mu_bri)
s_darci=dsqrt(s_darci2)
end

subroutine get_alf_bet_brink2(m,kap)
!коэффициенты, пересчитанные Дамиром
!v_tau_e=v_tau_i
!a*om_e-b*om_i=v_tau_e/sqrt(k)
use pgmod2
real(8) m,kap,k
if (m==0.999999d0) then
  k=0.490957d0   
  mu_bri=0.103646d1
  a_bri=102.5d0
  b_bri=270d0
elseif (m==0.9d0) then
  k=0.402716d-1
  mu_bri=0.333101d1
  a_bri=66d0
  b_bri=141d0
elseif (m==0.7d0) then
  k=0.972004d-2        
  mu_bri=0.138695d2
  a_bri=10d0
  b_bri=201d0
elseif (m==0.5d0) then
  k=0.187780d-02      
  mu_bri=0.728562d+02   
  a_bri=0.005d0
  b_bri=270.0d0
endif
!a_bri=d0
!b_bri=-dsqrt(mu_bri)
!a_bri=a_bri*d2
!b_bri=b_bri*d5
k_bri=k/kap**2
s_darci2=d1/(k_bri*mu_bri)
s_darci=dsqrt(s_darci2)
end

subroutine get_alf_bet_brink3(m,kap)
!коэффициенты, пересчитанные Дамиром для новых ГУ 
!v_s=a*sqrt(k)*om_s+b*k*dom_s/dn
!v_b=c*v_s+d*k*dom_s/dn
use pgmod2
real(8) m,kap,k
if (m==0.999999d0) then
  k=0.490957d0   
  mu_bri=0.103646d1
  a_bri=0.115873d1
  b_bri=0.108645d1
  c_bri=0.920719d0
  d_bri=0.159728d-2
elseif (m==0.9d0) then
  k=0.402716d-1
  mu_bri=0.333101d1
  a_bri=0.190222d1
  b_bri=0.186807d1
  c_bri=0.537961d0
  d_bri=0.427099d-1
elseif (m==0.7d0) then
  k=0.972004d-2        
  mu_bri=0.138695d2
  a_bri=0.244343d1
  b_bri=0.313159d1
  c_bri=0.344007d0  
  d_bri=0.122289d0 
elseif (m==0.5d0) then
  k=0.187780d-02      
  mu_bri=0.728562d+02   
  a_bri=0.337424d1
  b_bri=0.633721d1
  c_bri=0.219908d0
  d_bri=0.455616d0
endif
!a_bri=d0
!b_bri=-dsqrt(mu_bri)
!a_bri=a_bri*d2
!b_bri=b_bri*d5
k_bri=k/kap**2
s_darci2=d1/(k_bri*mu_bri)
s_darci=dsqrt(s_darci2)
end

subroutine get_alf_bet_brink4(m,kap)
!коэффициенты, пересчитанные Дамиром для новых ГУ 
!a*v_s=sqrt(k)*(om_b-b*om_s) !a=alf, b=gam
!sqrt(k)*d*om_s=c*v_s-v_b    !c=bet, d=delt
use pgmod2
real(8) m,kap,k
if (m==0.999999d0) then
  k=0.490957d0   
  mu_bri=0.103646d1
  a_bri=0.159369d-2
  b_bri=0.104610d1
  c_bri=0.921975d0
  d_bri=0.142786d-2
elseif (m==0.9d0) then
  k=0.402716d-1
  mu_bri=0.333101d1
  a_bri=0.139784d-1
  b_bri=0.534101d0
  c_bri=0.560825d0
  d_bri=0.434908d-1
elseif (m==0.7d0) then
  k=0.972004d-2        
  mu_bri=0.138695d2
  a_bri=0.171127d-1
  b_bri=0.183889d0
  c_bri=0.383057d0
  d_bri=0.954159d-1
elseif (m==0.5d0) then
  k=0.187780d-02      
  mu_bri=0.728562d+02   
  a_bri=0.156996d-1
  b_bri=0.339586d-1
  c_bri=0.291803d0
  d_bri=0.242590d0
endif
!a_bri=d0
!b_bri=-dsqrt(mu_bri)
!a_bri=a_bri*d2
!b_bri=b_bri*d5
k_bri=k/kap**2
s_darci2=d1/(k_bri*mu_bri)
s_darci=dsqrt(s_darci2)
end

subroutine get_alf_bet_brink5(m,kap)
!коэффициенты, пересчитанные Дамиром для новых ГУ 
!om_s=mub*om_b !a=alf=0, b=gam=1/mub
!sqrt(k)*d*om_s=c*v_s-v_b    !c=bet, d=delt
use pgmod2
real(8) m,kap,k
if (m==0.999999d0) then
  k=0.490957d0   
  mu_bri=0.103646d1
  c_bri=0.853390d0
  d_bri=0.658806d-2
elseif (m==0.9d0) then
  k=0.402716d-1
  mu_bri=0.333101d1
  c_bri=0.364920d0
  d_bri=0.146224d0
elseif (m==0.7d0) then
  k=0.972004d-2        
  mu_bri=0.138695d2
  c_bri=0.189327d0
  d_bri=0.194067d0
elseif (m==0.5d0) then
  k=0.187780d-02      
  mu_bri=0.728562d+02   
  c_bri=0.850515d-1
  d_bri=0.169807d0
endif
a_bri=d0
b_bri=d1/mu_bri

c_bri=c_bri/(d1+d_bri*dsqrt(mu_bri))
d_bri=d0

k_bri=k/kap**2
s_darci2=d1/(k_bri*mu_bri)
s_darci=dsqrt(s_darci2)
end

subroutine main5
!обтекание пористого цилиндра (модель Дарси) течением стокса
use pgmod2
use kuv_mod
real(8) ee,mcell_m,kapp,mcell_dx,mcell_r
logical square_cell,GUBivers,gu_closing,update_h2_by_mcell
integer(4) ku1,ku2,periodic,mode
real(8) par(5)
integer(4) nng,nbndl1
nmain=14
ku1=3
ku2=1

square_cell=.false.
update_h2_by_mcell=.false.
periodic=1
GUBivers=.false. !ГУ Биверса или мое
gu_closing=.true.

s_darci=1.0d0
alf_darci=1.0d0 !d1
bet_darci=1.0d0 !d1

mcell_m=0.5d0
kapp=10.0d0 !10,25
call get_alf_bet_bivers(mcell_m,kapp,0) !(m,kap,mode)

!bet_darci=d0

if (square_cell) then
  mode=3
  ee=0.96d0 !плотность упаковки 0.96   0.6
  !h=d1/dsqrt(d1-ee)  
  !h=dsqrt(pi)*d5*h
else
  periodic=0
  mode=0
endif 

if ((mode==0.or.mode==3).and.update_h2_by_mcell) then
  mcell_dx=d1/kapp
  mcell_r=dsqrt((d1-mcell_m)/pi)*mcell_dx
  h2=d1-mcell_dx*d5+mcell_r
endif

if (.not.square_cell) then
  nfun=34
  if (GUBivers) then
    call init_kuv_darci(h,s_darci,alf_darci)
    !call init_kuv_darci_om3(h,s_darci,alf_darci,bet_darci,1,0)   !1,0
    gu_closing=.false.
  else
    !call init_kuv_darci_om2(h,s_darci)
    call init_kuv_darci_om3(h,s_darci,alf_darci,bet_darci,1,0)   !1,0
  endif
endif

call pg_allocate_problems(1)
call pg_bind_problem(1)
call pg_allocate_domains(2)
!area_e
call pg_bind_domain(1)
call pg_set_domain_equation(ku1)
call pg_allocate_bounds(1)
call pg_bind_bound(1)
call init_geom2mode(mode,square_cell,nng,par)
!call init_geom2(square_cell)
call pg_geom_postprocessor
call init_gu2(square_cell,periodic>0,.false.,.false.,.false.)
nbndl1=gsbnd%nline
!area_i
call pg_bind_domain(2)
call pg_set_domain_equation(ku2)
call pg_allocate_bounds(1)
call pg_bind_bound(1)
call init_geom4
call pg_geom_postprocessor
call init_gu4_2
if (gu_closing) call pg_init_boundline_gu_stoks_darcy(1,1,nbndl1,2,1,1,s_darci,alf_darci,bet_darci)
call pg_get_matrix
if (.not.gu_closing) then
  if (GUBivers) then
    call closing_matrix5
  else
    call closing_matrix5_3
  endif
endif
if (periodic>0) call closing_matrix2_per(periodic)
call pg_solve
call pg_get_psioml
!call drw_solv4circ
call drw_p
call drw_solv4mode(mode,par,nng,square_cell)
end

subroutine main6
!тестирование решения внешней задачи обтекания цилиндра
use pgmod2
integer(4) ku
real(8) psi_inf
external psi_inf
nmain=22
nfun=41
ku=1
if (nfun==41) ku=3
call pg_allocate_problems(1)
call pg_bind_problem(1)
call pg_allocate_domains(1)
call pg_bind_domain(1)
call pg_set_domain_equation(ku)
call pg_allocate_bounds(1)
call pg_bind_bound(1)
call ga_init_geom_circle(d0,d0,d1,.true.,2*nj_)
call pg_geom_postprocessor
call pg_allocate_bound_gu
call pg_bind_boundline(1)
call pg_init_boundline_gu_val_const(1,1,d0)
if (nfun==41) call pg_init_boundline_gu_val_const(2,2,d0)
call pg_get_matrix_finf(psi_inf)
call pg_solve
call pg_get_psioml
call drw_solv1
end

subroutine testDistanceToCircleArc
use gen_mod
real(8) r,gam,DistanceToCircleArc
r=DistanceToCircleArc(0.1d0,-0.1d0,d0,d0,d1,-pi2,-pi,gam)
end

subroutine calc_pmap
use pgmod2
type(areatype), pointer :: a
integer(4) i,j
real(8) dpdl,eps,dp,ga_gradp_oneline_ds,y
external dpdl,write_ppoint
complex(8) xy
logical ceq_eps,eq_eps
integer(4) mode !заданное давлние
                !1 - только угол
                !2 - боковые стороны
eps=1.0d-6
if (.false.) then
  call pgani_createarea(50000)
else
  call pg_allocate_area(1)
  call pg_bind_areapart(1)
  call init_mesh4_2(.false.,.true.)
endif
call pg_areageom_postprocessor
a=>gsarea%a
allocate(parr(a%n),pknow(a%n))
pknow=.false.
mode=2
if (mode==1) then
  xy=c1*h
  do i=1,a%n
    pknow(i)=ceq_eps(xy,a%zm(i),eps)
    if (pknow(i)) parr(i)=d0
  enddo
elseif (mode==2) then
  call pg_bind_boundline(2)
  do i=1,gsbndl%npanel+1
    if (gsbndl%y(i)>d1) then
      j=(gsbndl%npanel+1+i)/2
      exit
    endif
  enddo
  y=gsbndl%y(j)
  dp=ga_gradp_oneline_ds(-h,y,h,y,1.0d-2)
  do i=1,a%n
    pknow(i)=eq_eps(h,dreal(a%zm(i)),eps)
    if (pknow(i)) then
      parr(i)=d0
    else
      pknow(i)=eq_eps(-h,dreal(a%zm(i)),eps)
      if (pknow(i)) parr(i)=dp
    endif
  enddo
endif
call intf_dxdy(a%zm,a%trm,a%n,a%ntr,gsarea%a%npe,parr,pknow,dpdl,d0)
call ga_drw_trmesh_func('p.dat','"X","Y","P"',write_ppoint)
end

function dpdl(x,y,xdir,ydir)
use pgmod2
real(8) dpdl,x,y,xdir,ydir,ff(2),scal_p
complex(8) df,lder
call pg_get_fun_xy_gradient(x,y,6,0,ff)
df=dcmplx(ff(2),-ff(1)) !градиент давления
lder=dcmplx(xdir,ydir)  !направляющий вектор
lder=lder/cdabs(lder)
dpdl=scal_p(df,lder)
end

subroutine write_ppoint(nf,i)
use pgmod2
integer(4) i,nf
WRITE(nf,"(3(F13.5))") dreal(gsarea%a%zm(i)), dimag(gsarea%a%zm(i)), parr(i)
end

subroutine calc_dp
use pgmod2
real(8) dp,ds,ga_gradp_oneline_ds,y,F,a,F1
ds=1.0d-2
y=(h2+h)*d5
dp=ga_gradp_oneline_ds(-h,y,h,y,ds)
F=dp*d2*h
F1=d1/F
a=h2/h
end

subroutine test_fast_cash_solve
use pgmod2
real(8) angle,t,maxt,maxt_all
integer(4) i,j,k,i1,j1,ku
type(TArea), pointer :: a,a2
type(TCash), pointer :: c,c2
type(TCashIntegral), pointer :: ci,ci2
type(TCashValue), pointer :: cv,cv2
logical, pointer :: need_cash(:)

angle=-d1
ku=1
k_helm=d2
!gs_use_parallel_matrix_build=0

call pg_allocate_problems(1)
call pg_bind_problem(1)

!gs%const%ds0numint=1.0d-4

call pg_allocate_domains(3)
call pg_bind_domain(1)
call pg_set_domain_equation(ku)
call pg_set_areaconst(1,k_helm)
call pg_allocate_bounds(1)
call pg_bind_bound(1)
call init_geom3
call pg_geom_postprocessor
call pg_allocate_area(1)
call pg_bind_areapart(1)
call init_mesh3(d1,d1,nj_,nj_)
call pg_areageom_postprocessor
call pg_allocate_area_gu
call pg_copy_subdomain_dxdy_rotate(1,2,d0,d0,angle,d0,d0)
gsarea%type_rotate=-1
gsarea%a_ref=>null()
call pg_copy_subdomain_dxdy_rotate(1,3,d0,d0,angle,d0,d0)
gs_cash_is_solve=.false.
call init_gs_time
do i=1,gs%na
  a=>gs%a(i)
  a%need_cash_bb=.true.
  a%need_cash_ba=.true.
  forall (j=0:11) a%need_cash_ab(j)=int_type_area(j)==1
  a%need_cash_aa=a%need_cash_ab
  call initCash_area(i)
  !повторяем, так как в initCash_area все зануляется
  a%need_cash_bb=.true.
  a%need_cash_ba=.true.
  forall (j=0:11) a%need_cash_ab(j)=int_type_area(j)==1
  a%need_cash_aa=a%need_cash_ab
enddo
if (gs_cash_presolve.and.gs_cash_is_solve) call print_total_time
a=>gs%a(2)
a2=>gs%a(3)
maxt_all=d0
do i=1,2
  select case (i)
  case (1)
    c=>a%bnd(1)%self_cash
    c2=>a2%bnd(1)%self_cash
    print*, "bound cash"
  case (2)
    c=>a%a%self_cash
    c2=>a2%a%self_cash
    print*, "area_cash"
  endselect
  do j=1,2
    select case (j)
    case (1)
      ci=>c%bnd(1)
      ci2=>c2%bnd(1)
      print*, "bound cp"
      select case (i)
      case (1)
        need_cash=>a%need_cash_bb
      case (2)
        need_cash=>a%need_cash_ab
      endselect
    case (2)
      ci=>c%area
      ci2=>c2%area
      print*, "area cp"
      select case (i)
      case (1)
        need_cash=>a%need_cash_ba
      case (2)
        need_cash=>a%need_cash_aa
      endselect
    endselect
    do k=0,max_int
      if (.not.need_cash(k)) cycle
      maxt=d0
      cv=>ci%intvals(k)
      cv2=>ci2%intvals(k)
      do i1=1,ubound(cv%i,1)
        do j1=1,ubound(cv%i,2)
          t=dabs(cv%i(i1,j1)-cv2%i(i1,j1))
          if (t>maxt) maxt=t
        enddo
      enddo
      write(*,"('nf=',i0,' err=',F9.5)") k,maxt
      if (maxt>maxt_all) maxt_all=maxt
    enddo
  enddo 
enddo
print*,"errall=",maxt_all
end

subroutine test_fast_cash_solve_dual
use pgmod2
real(8) angle,t,maxt,maxt_all
integer(4) i,j,k,i1,j1,ku
type(TArea), pointer :: a,a2
type(TCash), pointer :: c,c2
type(TCashIntegral), pointer :: ci,ci2
type(TCashValue), pointer :: cv,cv2
logical, pointer :: need_cash(:)
logical is_dual(0:max_int)

angle=d0
ku=1
!gs_use_parallel_matrix_build=0

is_dual=int_type_dual/=-1

call pg_allocate_problems(1)
call pg_bind_problem(1)

!gs%const%ds0numint=1.0d-4

call pg_allocate_domains(2)
call pg_bind_domain(1)
call pg_set_domain_equation(ku)
call pg_allocate_bounds(1)
call pg_bind_bound(1)
call init_geom3
call pg_geom_postprocessor
call pg_allocate_area(1)
call pg_bind_areapart(1)
call init_mesh3(d1,d1,nj_,nj_)
call pg_areageom_postprocessor
call pg_allocate_area_gu
call pg_copy_subdomain_dxdy_rotate(1,2,d0,d0,angle,d0,d0)
gsarea%type_rotate=-1
gsarea%a_ref=>null()
gs_cash_is_solve=.false.
do i=1,gs%na
  gs_use_dual_integral_solve=i==1
  call init_gs_time
  a=>gs%a(i)
  a%need_cash_bb=is_dual
  a%need_cash_ba=is_dual
  forall (j=0:11) a%need_cash_ab(j)=int_type_area(j)==1.and.is_dual(j)
  a%need_cash_aa=a%need_cash_ab
  call initCash_area(i)
  !повторяем, так как в initCash_area все зануляется
  a%need_cash_bb=is_dual
  a%need_cash_ba=is_dual
  forall (j=0:11) a%need_cash_ab(j)=int_type_area(j)==1.and.is_dual(j)
  a%need_cash_aa=a%need_cash_ab
  if (gs_cash_presolve.and.gs_cash_is_solve) call print_total_time
enddo
a=>gs%a(1)
a2=>gs%a(2)
maxt_all=d0
do i=1,2
  select case (i)
  case (1)
    c=>a%bnd(1)%self_cash
    c2=>a2%bnd(1)%self_cash
    print*, "bound cash"
  case (2)
    c=>a%a%self_cash
    c2=>a2%a%self_cash
    print*, "area_cash"
  endselect
  do j=1,2
    select case (j)
    case (1)
      ci=>c%bnd(1)
      ci2=>c2%bnd(1)
      print*, "bound cp"
      select case (i)
      case (1)
        need_cash=>a%need_cash_bb
      case (2)
        need_cash=>a%need_cash_ab
      endselect
    case (2)
      ci=>c%area
      ci2=>c2%area
      print*, "area cp"
      select case (i)
      case (1)
        need_cash=>a%need_cash_ba
      case (2)
        need_cash=>a%need_cash_aa
      endselect
    endselect
    do k=0,max_int
      if (.not.need_cash(k)) cycle
      maxt=d0
      cv=>ci%intvals(k)
      cv2=>ci2%intvals(k)
      do i1=1,ubound(cv%i,1)
        do j1=1,ubound(cv%i,2)
          t=dabs(cv%i(i1,j1)-cv2%i(i1,j1))
          if (t>maxt) maxt=t
        enddo
      enddo
      write(*,"('nf=',i0,' err=',F9.5)") k,maxt
      if (maxt>maxt_all) maxt_all=maxt
    enddo
  enddo 
enddo
print*,"errall=",maxt_all
end

subroutine testselectif
use pgmod
integer(8) i
integer(4) j
real(8) s,f,testselectif1 !,testselectif2
s=d0
call init_gs_time
do i=1,1000000000
  j=mod(i,10)
  f=testselectif1(j)
  s=s+f
enddo
call print_total_time
end

function testselectif1(i) result (f)
real(8) f
integer(4) i
if (i==0) then
  f=dsin(0.0d0)
elseif (i==1) then
  f=dsin(1.0d0)
elseif (i==2) then
  f=dsin(2.0d0)
elseif (i==3) then
  f=dsin(3.0d0)
elseif (i==4) then
  f=dsin(4.0d0)
elseif (i==5) then
  f=dsin(5.0d0)
elseif (i==6) then
  f=dsin(6.0d0)
elseif (i==7) then
  f=dsin(7.0d0)
elseif (i==8) then
  f=dsin(8.0d0)
elseif (i==9) then
  f=dsin(9.0d0)
endif
end

function testselectif2(i) result(f)
real(8) f
integer(4) i
select case (i)
case (0) 
  f=dsin(0.0d0)
case (1)
  f=dsin(1.0d0)
case (2)
  f=dsin(2.0d0)
case (3)
  f=dsin(3.0d0)
case (4)
  f=dsin(4.0d0)
case (5)
  f=dsin(5.0d0)
case (6)
  f=dsin(6.0d0)
case (7)
  f=dsin(7.0d0)
case (8)
  f=dsin(8.0d0)
case (9)
  f=dsin(9.0d0)
end select
end

subroutine test_cash_2_problem
use pgmod2

!gs_use_parallel_matrix_build=0

call pg_allocate_problems(3)
!уравнение Лапласа
nmain=1
nfun=1 !1 - xx
call pg_bind_problem(1)
call pg_allocate_domains(1)
call pg_bind_domain(1)
call pg_set_domain_equation(1)
call pg_allocate_bounds(1)
call pg_bind_bound(1)
call init_geom3
call pg_geom_postprocessor
call init_gu3
call pg_get_matrix
call pg_solve
call pg_get_psioml
call drw_solv3
!уравнение Лапласа
nfun=2 !2 - d1-xx-yy+xx*yy
call pg_bind_problem(2)
call pg_allocate_domains(1)
call pg_bind_domain(1)
call pg_set_domain_equation(1)
call pg_copy_geom_domain_dxdy_rotate(1,1,d0,d0,d0,d0,d0)
call pg_set_cash_ref(1,1)
call init_gu3
call pg_get_matrix
call pg_solve
call pg_get_psioml
call drw_solv3
!багармоническое уравнение
nmain=4
nfun=19 !19 - get_fund=xx**3
call pg_bind_problem(3)
call pg_allocate_domains(1)
call pg_bind_domain(1)
call pg_set_domain_equation(3)
call pg_copy_geom_domain_dxdy_rotate(1,1,d0,d0,d0,d0,d0)
call pg_set_cash_ref(1,1)
call init_gu3_2
call pg_get_matrix
call pg_solve
call pg_get_psioml
call drw_solv3
end
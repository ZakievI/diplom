subroutine main_kuvbri
!обтекание пористого цилиндра (уравнение бринкмана) течением стокса
use pgmod2
integer(4) ku1,ku2,ku
nmain=18
nfun=35
ku=3
ku1=8
ku2=14 !13,14
!ku1=11
!ku2=13
s_darci=2.0d0
s_darci2=s_darci**2
!call init_kuv_bri_stech(h,s_darci)
!call init_kuv_bri3(h,s_darci) 

call pg_allocate_problems(1)
call pg_bind_problem(1)
call pg_allocate_domains(2)
!area_e
call pg_bind_domain(1)
call pg_set_domain_equation(ku)
call pg_allocate_bounds(1)
call pg_bind_bound(1)
call init_geom2(.false.)
call pg_geom_postprocessor
call init_gu2(.false.,.false.,.false.,.false.,.false.)
!area_i
call pg_bind_domain(2)
call pg_set_domain_equation_syst(ku1,ku2)
call pg_set_areaconst(1,s_darci) 
call pg_allocate_bounds(1)
call pg_bind_bound(1)
call init_geom4
call pg_geom_postprocessor
call init_gu4_3
call pg_allocate_area(1)
call pg_bind_areapart(1)
call init_mesh4
call pg_areageom_postprocessor
call pg_allocate_area_gu
call init_meshval_bri
call pg_get_matrix
call closing_matrix_bri
call pg_solve
call pg_get_psioml
call drw_solv4circ
end

subroutine main_kuvbri5
!обтекание пористого цилиндра (уравнение бринкмана) течением стокса через единую функцию √рина
use pgmod2
integer(4) ku1,ku2,nng,mode,periodic
real(8) ee,mcell_m,kapp !,psi1,pg_get_fun_xy
logical square_cell
real(8) par(5)
nmain=18

!hv=d2

solve_only_psi=.true.

square_cell=.false.
periodic=1
if (square_cell) then 
  ee=0.96d0 !плотность упаковки 0.96   0.6
  !h=d1/dsqrt(d1-ee)  
  !h=dsqrt(pi)*d5*h
else
  periodic=0
endif

mode=0 !0 - circle
       !1 - ellipse
	   !2 - virus
	   !3 - square
	   !4 - triangle (poly)

ku1=3  !стокс
ku2=15 !бринкман
!ku2=1 !дарси-мое √”

s_darci=2.0d0
!s_darci=dsqrt(10.0d0) 
s_darci2=s_darci**2

mcell_m=0.9d0
kapp=10.0d0 !10,25

!call get_alf_bet_bivers(mcell_m,kapp,1) !(m,kap,mode)
!call get_alf_bet_brink(mcell_m,kapp)
!call get_alf_bet_brink2(mcell_m,kapp)
call get_alf_bet_brink5(mcell_m,kapp)

if (mode==0.and.(.not.square_cell)) then
  if (ku2==15) then
    !call init_kuv_bri_stech(h,s_darci)
    call init_kuv_bri3(h,s_darci,k_bri,a_bri,b_bri,c_bri,d_bri) 
  else
    call init_kuv_darci_om2(h,s_darci)
  endif
  nfun=35
endif

call pg_allocate_problems(1)
call pg_bind_problem(1)

!gs%const%ds0numint=1.0d-4

call pg_allocate_domains(2)
!area_e
call pg_bind_domain(1)
call pg_set_domain_equation(ku1)
call pg_allocate_bounds(1)
call pg_bind_bound(1)
call init_geom2mode(mode,square_cell,nng,par)
call pg_geom_postprocessor
call init_gu2(square_cell,periodic>0,.false.,.false.,.false.)
!area_i
call pg_bind_domain(2)
call pg_set_domain_equation(ku2)
if (ku2==15) then
  call pg_set_areaconst(1,s_darci)
  call pg_set_areaconst(4,k_bri)
  call pg_set_areaconst(5,mu_bri)
endif
call pg_allocate_bounds(1)
call pg_bind_bound(1)
call init_geom4
call pg_geom_postprocessor
if (ku2==15) then
  call init_gu4_3
else
  call init_gu4_2
endif

!only_analit=.true.
!call drw_solv4mode(mode,par,nng,square_cell)

call main_kuvbri5_2(square_cell,periodic,mode,par,nng,ku2,mcell_m)
!call main_kuvbri5_4(periodic,ku2)
!psi1=pg_get_fun_xy(d0,d1,1,d0,d0,2)
end

subroutine main_kuvbri5_1(periodic,ku2)
!пересоставление системы
integer(4) ku2,periodic
call pg_get_matrix
if (ku2==15) then
  call closing_matrix_bri
else
  call closing_matrix5_3
endif
if (periodic>0) call closing_matrix2_per(periodic)
call pg_solve
end

subroutine main_kuvbri5_2(square_cell,periodic,mode,par,nng,ku2,mcell_m)
!расчет течени€ в €чейке с выводом в файл
use pgmod2
logical square_cell
integer(4) periodic
integer(4) ku2
real(8) par(5),r_ellipse,r_circle,r_virus,r_poly,mcell_m
external r_ellipse,r_circle,r_virus,r_poly
integer(4) nng,mode
call main_kuvbri5_1(periodic,ku2)
call pg_get_psioml
!call save_om_domdn(mcell_m,mode)
call drw_solv4mode(mode,par,nng,square_cell)
call get_f_kuvbri5(square_cell,mode)
end

subroutine main_kuvbri5_4(periodic,ku2)
!построение графика Q(s)
use pgmod2
integer(4) i0,nds0,ku2
integer(4) periodic
type(TBound), pointer :: b
b=>gs%a(1)%bnd(1)
nds0=50
OPEN (1,FILE='test5s.dat')
write(1,*) 'TITLE = "Q(s)"'
write(1,*) 'VARIABLES = "s", "Q"'
write(1,"('ZONE T=""area"", I=', i0, ', F=POINT')") nds0+1
write(1,"(F4.1, ' ', F9.6)") d0, d1
do i0=1,nds0 !шаг по s
  s_darci=i0*0.2d0
  s_darci2=s_darci**2
  gs%a(2)%const%k_helm=s_darci
  call main_kuvbri5_1(periodic,ku2)
  write(1,"(F4.1, ' ', F9.6)") s_darci,b%psiom(b%npanel-nj_/2,1)
enddo
close(1)
end

subroutine main_kuvbri6
!обтекание пористого цилиндра (уравнение бринкмана) в пористой среде (уравнение Ѕринкмана) через единую функцию √рина
use pgmod2
integer(4) ku1,ku2
nmain=19
ku1=15
!ku1=3
!ku2=15
ku2=3
s_darci=5.0d0
s_darci_2=2.00d0
if (ku1==15.and.ku2==15) then
  call init_kuv_bribri(h,s_darci,s_darci_2)
  nfun=36
endif
if (ku1==3) s_darci=d0
if (ku2==3) s_darci_2=d0
s_darci2=s_darci**2
s_darci_22=s_darci_2**2

call pg_allocate_problems(1)
call pg_bind_problem(1)
call pg_allocate_domains(2)
!area_e
call pg_bind_domain(1)
call pg_set_domain_equation(ku1)
if (ku1==15) call pg_set_areaconst(1,s_darci)
call pg_allocate_bounds(1)
call pg_bind_bound(1)
call init_geom2(.false.)
call pg_geom_postprocessor
call init_gu2(.false.,.false.,.false.,.false.,.false.)
!area_i
call pg_bind_domain(2)
call pg_set_domain_equation(ku2)
if (ku2==15) call pg_set_areaconst(1,s_darci_2)
call pg_allocate_bounds(1)
call pg_bind_bound(1)
call init_geom4
call pg_geom_postprocessor
call init_gu4_3
call pg_get_matrix
call closing_matrix_bribri
call pg_solve
call pg_get_psioml
call drw_solv4circ
end

subroutine main_kuvbri7
!тестирование перехода течени€ бринкмана в стокса в €чейке кувабара (обтекание непроницаемого цилиндра)
use pgmod2
real(8), parameter :: nn=101
real(8) qpsi(2,nn),y,u,v,om,psi1,qu(2,nn),qom(2,nn)
integer(4) i
s_darci=0.01d0
call init_kuv(h)
!s_darci=20.0d0
!call init_kuv_d(h)
do i=1,nn
  y=(i-d1)/(nn-d1)*(h-d1)+d1
  call Get_kuv_uvom_psi(d0,y,u,v,om,psi1)
  qpsi(1,i)=psi1
  qu(1,i)=u
  qom(1,i)=om
enddo
call init_kuv_brie2(h,s_darci)
do i=1,nn
  y=(i-d1)/(nn-d1)*(h-d1)+d1
  call Get_kuv_uvom_psi_bri(d0,y,u,v,om,psi1)
  qpsi(2,i)=psi1
  qu(2,i)=u
  qom(2,i)=om
enddo
end

subroutine save_om_domdn(mcellm,mode)
use pgmod2
integer(4) i,mode,j,j1
real(8) k,cc(0:1),mcellm
real(8), allocatable :: om(:), dom(:), domk(:), om_(:), dom_(:)
logical, allocatable :: maska(:)
k=dsqrt(k_bri*mu_bri)
call pg_bind_domain(2)
call pg_bind_bound(1)
call pg_bind_boundline(1)
allocate(om(gsbndl%npanel),dom(gsbndl%npanel),domk(gsbndl%npanel),om_(gsbndl%npanel),dom_(gsbndl%npanel))
om=gsbnd%psiom(1:gsbndl%npanel,3)
dom=gsbnd%psiom(1:gsbndl%npanel,4)*k
if (mode==3) then
  allocate(maska(gsbndl%npanel))
  j=gsbndl%npanel/4
  maska=.true.
  j1=6
  maska(j-j1+1:j+j1)=.false.
  j=j*3
  maska(j-j1+1:j+j1)=.false.
  j=0
  do i=1,gsbndl%npanel
    if (maska(i)) then
      j=j+1
      om_(j)=om(i)
      dom_(j)=dom(i)
    endif
  enddo
  deallocate(maska)
else
  j=gsbndl%npanel
  om_=om
  dom_=dom
endif
call LS_approx(1,j,om_,dom_,cc)
domk=dom/cc(1)
OPEN (1,FILE='om.dat')
write(1,"('TITLE = ""om(s) m=',F10.6,'""')") mcellm
write(1,*) 'VARIABLES = "s", "om", "domdn", "domdn_k"'
write(1,"('ZONE T=""om(s)"", I=', i4, ', F=POINT')") gsbndl%npanel
do i=1,gsbndl%npanel
  write(1,"(4(E15.6))") gsbnd%sc(i),om(i),dom(i),domk(i)
enddo
close(1)
deallocate(om,dom,domk,om_,dom_)
OPEN (1,FILE='kdom.dat')
write(1,"(2(E15.6))") (cc(i),i=0,1)
close(1)
end

subroutine get_f_kuvbri5(square_cell,mode)
use pgmod2
real(8) f,fy,yy,ga_gradp_oneline_ds,f1
logical square_cell
integer(4) mode
call pg_bind_domain(1)
if (square_cell) then
  yy=h*d5
  f=ga_gradp_oneline_ds(-h,yy,h,yy,1d-3)*h
else
  call ga_get_force_circle_bl(1,2,f,fy)
endif
f=f*2
if (mode==0) then
  call pg_bind_domain(2)
  call ga_get_force_circle_bl(1,1,f1,fy)
  f1=f1*2
endif
end
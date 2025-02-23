subroutine main_p1
!уравнение пуассона
use pgmod2
integer(4) ku
ku=4
nmain=2
nfun=4   !3 - get_funb=xx*xx
         !4 - get_funb=dsin(xx)

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
!call ga_drw_trmesh(1)
call pg_allocate_area_gu
call init_meshval3
call pg_get_matrix
call pg_solve
call drw_solv3
end

subroutine main_p1_2_first
use pgmod2
integer(4) ku
ku=4
nmain=2
nfun=4   !3 - get_funb=xx*xx
         !4 - get_funb=dsin(xx)
call pg_set_domain_equation(ku)
end

subroutine main_p1_2
!уравнение пуассона
!тестирование сабдоменов
use pgmod2
external main_p1_2_first
call main_subdomain(main_p1_2_first,1)
end

subroutine main_p2
!уравнение Гельмгольца
use pgmod2
integer(4) ku
ku=6 !5,6
nmain=3
nfun=13
       !ku=5
       !4 - get_funb=dsin(xx)
       !5 - get_funb=dexp(xx)
	   !6 - get_funb=dsin(k_helm*xx)
       !7 - get_funb=yy*dsin(xx)
       !8 - get_funb=yy*dsin(k_helm*xx)
       !9 - get_funb=yy*dexp(k_helm*xx) 
       !10- get_funb=dsin(xx*yy) 
       !17- get_funb=(xx+d1)**3
       
       !ku=6
       !11 - get_funb=(xx*yy)**2
       !12 - get_funb=dsin(xx*yy) 
       !13 - get_funb=xx**2*dsin(yy)
       !14 - get_funb=xx*yy*(d1-dexp((xx-d1)/e_helm))*(d1-dexp((yy-d1)/e_helm))
       !14      f_helmgolz=d1/e_helm
       !15      f_helmgolz=(xx+yy)/e_helm
       !16      f_helmgolz=dcos(xx*yy)/e_helm
       !46 - get_funb=sin(x)
             !f_helmgolz={0 : x<0.5; 1 x>0.5}
             !f_puasson={-sin(x) : x<0.5; 0 x>0.5}

e_helm=1.0d0
k_helm=d2

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
call drw_solv3
end


subroutine main_p4
!бигармоническое уравнение, старый подход
use pgmod2
integer(4) ku
nmain=4
nfun=18 !18 - get_fund=xx
        !19 - get_fund=xx**3
ku=3
call pg_allocate_problems(1)
call pg_bind_problem(1)
call pg_allocate_domains(1)
call pg_bind_domain(1)
call pg_set_domain_equation(ku)
call pg_allocate_bounds(1)
call pg_bind_bound(1)
call init_geom3
call pg_geom_postprocessor
call init_gu3_2

gs%m%matrix_store_type=2

call pg_get_matrix
call pg_solve
call pg_get_psioml
call drw_solv3
end

subroutine main_p4_2_first
use pgmod2
integer(4) ku
nmain=4
nfun=19 !18 - get_fund=xx
        !19 - get_fund=xx**3
ku=3
call pg_set_domain_equation(ku)
end

subroutine main_p4_2
!уравнение пуассона
!тестирование сабдоменов
use pgmod2
external main_p4_2_first
call main_subdomain(main_p4_2_first,1)
end

subroutine main_p5
!уравнение Бринкмана через единую функцию Грина
use pgmod2
integer(4) ku
nmain=7
nfun=25 !25 - yy*dexp(k_helm*xx) 

ku=15
k_helm=d2

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
call init_gu3_2
call pg_get_matrix
call pg_solve
call pg_get_psioml
call drw_solv3
end

subroutine main_p6
!неоднородное бигармоническое уравнение
use pgmod2
integer(4) ku
ku=21
nmain=23
nfun=22   !22 - get_fund=xx**4   get_fund2=12*xx**2

call pg_allocate_problems(1)
call pg_bind_problem(1)
call pg_allocate_domains(1)
call pg_bind_domain(1)
call pg_set_domain_equation(ku)
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

subroutine main_p7
!неоднородное бигармоническое уравнение с неизвестными в правой части
use pgmod2
integer(4) ku
logical eq_var(7)
ku=22
nmain=24
nfun=48   !42 - get_fund=xx**4*y   get_fund2=12*xx**2*y
            !(f(1)=24*y-x**4*y; f(2)=1; f(3)=-x**2; f(4)=x; f(5)=-y; f(6)=x**3; f(7)=-x**2*y)
          !43 - get_fund=x**2+y**2   get_fund2=4
            !(f(1)=-4x*y*; f(4)=y; f(5)=x);
          !44 - get_fund=x**2+y**2   get_fund2=4
            !(f(4)=y; f(5)=-x)		   
          !45 - get_fund=y*sin(x) get_fund2=-y*sin(x)
            !(f(1)=x**2*y**2*cos(x); f(2)=x; f(3)=cos(x); f(4)=sin(x); f(5)=y; f(6)=x**2*y; f(7)=x*y)
          !48 - get_fund=y*sin(x) get_fund2=-y*sin(x) (часть коэффициентов при неизвестных в области = 0)
             !(f(1)={0 : x<0.5; y*(sin(x)-x*cos(x)) : x>0.5} f(4)={0 : x<0.5; x : x>0.5} f(5)={y : x<0.5; 0 : x>0.5})
eq_var=.false.
select case (nfun)
case (42,45)
  eq_var=.true.
case (43,48)
  eq_var(1)=.true.
  eq_var(4)=.true.
  eq_var(5)=.true.
case (44)
  eq_var(4)=.true.
  eq_var(5)=.true.
endselect

call pg_allocate_problems(1)
call pg_bind_problem(1)
call pg_allocate_domains(1)
call pg_bind_domain(1)
call pg_set_domain_equation(ku)
call pg_set_eq_var(eq_var,7)
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
call find_error3
end

subroutine main_p7_2_first
use pgmod2
integer(4) ku
logical eq_var(7)
ku=22
nmain=24
nfun=48   !42 - get_fund=xx**4*y   get_fund2=12*xx**2*y
            !(f(1)=24*y-x**4*y; f(2)=1; f(3)=-x**2; f(4)=x; f(5)=-y; f(6)=x**3; f(7)=-x**2*y)
          !43 - get_fund=x**2+y**2   get_fund2=4
            !(f(1)=-4x*y*; f(4)=y; f(5)=x);
          !44 - get_fund=x**2+y**2   get_fund2=4
            !(f(4)=y; f(5)=-x)		   
          !45 - get_fund=y*sin(x) get_fund2=-y*sin(x)
            !(f(1)=x**2*y**2*cos(x); f(2)=x; f(3)=cos(x); f(4)=sin(x); f(5)=y; f(6)=x**2*y; f(7)=x*y)
          !48 - get_fund=y*sin(x) get_fund2=-y*sin(x) (часть коэффициентов при неизвестных в области = 0)
             !(f(1)={0 : x<0.5; y*(sin(x)-x*cos(x)) : x>0.5} f(4)={0 : x<0.5; x : x>0.5} f(5)={y : x<0.5; 0 : x>0.5})
eq_var=.false.
select case (nfun)
case (42,45)
  eq_var=.true.
case (43,48)
  eq_var(1)=.true.
  eq_var(4)=.true.
  eq_var(5)=.true.
case (44)
  eq_var(4)=.true.
  eq_var(5)=.true.
endselect

call pg_set_domain_equation(ku)
call pg_set_eq_var(eq_var,7)
end

subroutine main_p7_2
!неоднородное бигармоническое уравнение с неизвестными в правой части
!тестирование сабдоменов
use pgmod2
external main_p7_2_first
call main_subdomain(main_p7_2_first,1)
end

subroutine find_error3
use pgmod2
integer(4) i,j,n,m
real(8) x,y,f,f1,om,om1,eps,epsom,epsa,epsoma,epse,epsome
real(8) get_fund,get_fund2,pg_get_fun,pg_get_fun2
n=10
epsa=d0
epse=d0
epsoma=d0
epsome=d0
m=0
do i=1,n-1
  x=(i-d0)/n
  do j=1,n-1
    y=(j-d0)/n
    f=get_fund(x,y)
    f1=pg_get_fun(x,y)
    om=get_fund2(x,y)
    om1=pg_get_fun2(x,y)
    eps=dabs(f-f1)
    epsom=dabs(om-om1)
    if (eps>epsa) epsa=eps
    if (epsom>epsoma) epsoma=epsom
    epse=epse+eps
    epsome=epsome+epsom
    m=m+1
  enddo
enddo
epse=epse/m
epsome=epsome/m
write(*,"('epsa=',E13.5)") epsa
write(*,"('epsoma=',E13.5)") epsoma
write(*,"('epse=',E13.5)") epse
write(*,"('epsome=',E13.5)") epsome
OPEN (1,FILE='err.dat')
write(1,"(I0,4(E13.5))") nj_,epsa,epsoma,epse,epsome
close(1)
end

subroutine main_p8
!уравнение Пуассона с неизвестными в правой части
use pgmod2
integer(4) ku
logical eq_var(4)
ku=26
nmain=25
nfun=28   !49 - get_funb=y*sin(x)
            !(f(1)=y**2*cos(x)**2-y*sin(x); f(2)=x; f(3)=-y*cos(x); f(4)=-x*y)
          !8 - get_funb=yy*dsin(k_helm*xx) - однородное уравнение Гельмгольца 
            !(f(1)=0; f(2)=-k_helm**2; f(3)=0; f(4)=0)
          !13 - get_funb=xx**2*dsin(yy) - неоднородное уравнение Гельмгольца 
            !(f(1)=(d2+xx**2*(yy-d1))*dsin(yy); f(2)=-y; f(3)=0; f(4)=0)
          !26 - y*y+2*x      (vx=d1;vy=d0) - однородное уравнение переноса  
            !(f(1)=0; f(2)=0; f(3)=d1; f(4)=d0)
          !28 - xx**2+yy**2  (vx=xx;vy=yy**2) - неоднородное уравнение переноса  
            !(f(1)=d2*(d2-xx**2-yy**3); f(2)=0; f(3)=x; f(4)=y**2)
eq_var=.false.
select case (nfun)
case (49)
  eq_var=.true.
case (8)
  k_helm=d2
  eq_var(2)=.true.
case (13)
  eq_var(1:2)=.true.
case (26)
  eq_var(3)=.true.
case (28)
  eq_var(1)=.true.
  eq_var(3:4)=.true.
endselect

call pg_allocate_problems(1)
call pg_bind_problem(1)
call pg_allocate_domains(1)
call pg_bind_domain(1)
call pg_set_domain_equation(ku)
call pg_set_eq_var(eq_var,4)
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
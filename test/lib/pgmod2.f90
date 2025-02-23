module pgmod2
use pgmod
use dnsmod_base
integer(4),PARAMETER :: nmax3=10
integer(4),PARAMETER :: nmax=1521
integer(4),PARAMETER :: nmax2=nmax*2
real(8) h,h2,hv
integer(4) nj_,njr_,npe_,nj_2
real(8) e_helm,kkr,s_darci,alf_darci,s_darci2,r2,bet_darci,a_bri,b_bri,mu_bri,c_bri,d_bri
real(8) s_darci_2,s_darci_22,k_bri
real(8) Pe         !число Пекле
real(8) k_helm
real(8) k_oss
real(8) square_shift_dy
real(8) maxerre
real(8) rpoly_x(nmax3),rpoly_y(nmax3),rpoly_g(nmax3),rpoly_ksi(nmax3),rpoly_eta(nmax3),rpoly_tt(nmax3)
real(8), allocatable :: parr(:)
logical, allocatable :: pknow(:)
integer(4) rpoly_n
logical drw_e_tr
logical need_thickening
integer(4) n_thickening
real(8) slip_k
integer nsx,nsy
logical use_subdomain,only_analit
type(mcelltype) mcell
logical solve_only_psi
integer(4) nfun
           !1 - xx
		   !2 - d1-xx-yy+xx*yy
		   !3 - xx*xx
		   !4 - dsin(xx)
		   !5 - dexp(xx)
		   !6 - dsin(k_helm*xx)
           !7 - yy*dsin(xx)
           !8 - yy*dsin(k_helm*xx)
           !9 - yy*dexp(k_helm*xx) 
           !10 - dsin(xx*yy)  в однородном уравнении гельмгольца
           !11 - (xx*yy)**2
           !12 - dsin(xx*yy)  в неоднородном уравнении гельмгольца
           !13 - xx**2*dsin(yy)
           !14 - xx*yy*(d1-dexp((xx-d1)/e_helm))*(d1-dexp((yy-d1)/e_helm)) 
           !14      f_helmgolz=d1/e_helm
           !15      f_helmgolz=(xx+yy)/e_helm
           !16      f_helmgolz=dcos(xx*yy)/e_helm
           !17 - (xx+d1)**3
		   !18 - xx
		   !19 - xx**3
		   !20 - yy*dsin(k_helm*xx)
           !21 - yy*dexp(k_helm*xx) 
		   !22 - xx**4   12*xx**2
		   !23 - xx**4   6*xx
		   !24 - dsin(xx*yy)   dsin(xx*yy)
		   !25 - yy*dexp(k_helm*xx) 
		   !26 - y*y+2*x      (vx=d1;vy=d0)
		   !27 - y*dexp(x)    (vx=d1;vy=d0)
		   !28 - xx**2+yy**2  (vx=xx;vy=yy**2)
		   !29 - xx+yy        (vx=xx;vy=yy)
		   !30 - осаждение взвеси в канале    (vx=d1;vy=d0)
		   !31 - осаждение взвеси в канале    (vx=течение Пуазеля;vy=d0)
		   !32 - Дарси с переменной проницаемостью (S^2(r)=r**(-s_darci))
		   !33 - аналитическое решение течение стокса в ячейке кувабара
		   !34 - аналитическое решение течение стокса в ячейке кувабара с пористым цилиндром модель Дарси
		   !35 - аналитическое решение течение стокса в ячейке кувабара с пористым цилиндром модель Бринкмана
		   !36 - аналитическое решение течение Бринкмана в ячейке кувабара с пористым цилиндром модель Бринкмана
		   !37 - y*y+x*x      (vx=0,vy=x*y)
		   !38 - x            (vx=0,vy=-1/y)(k_oss=d1)
		   !39 - y**2/2       (k_oss=-d1)
		   !40 - обтекание цилиндра поленциальным безграничным потоком
		   !41 - обтекание цилиндра безграничным потоком Стокса
       !42 - x**4*y (f(1)=24*y-x**4*y; f(2)=1; f(3)=-x**2; f(4)=x; f(5)=-y; f(6)=x**3; f(7)=-x**2*y)
       !43 - x**2+y**2 (f(1)=-4x*y*; f(4)=y; f(5)=x)		   
       !44 - x**2+y**2 (f(4)=y; f(5)=-x)	
       !45 - y*sin(x) (f(1)=x**2*y**2*cos(x); f(2)=x; f(3)=cos(x); f(4)=sin(x); f(5)=y; f(6)=x**2*y; f(7)=x*y)
       !46 - sin(x) - в неоднородном уравнении гельмгольца (часть коэффициентов при неизвестных в области = 0)
             !f_helmgolz={0 : x<0.5; 1 : x>0.5}
             !f_puasson={-sin(x) : x<0.5; 0 : x>0.5}
       !47 - x**2+y**2 - в неоднородном уравнении переноса (часть коэффициентов при неизвестных в области = 0)
             !vx={0 : x<0.5; x : x>0.5} vy={y : x<0.5; 0 : x>0.5}
       !48 - y*sin(x) - бигарм.ур с неизвестными в правой части (часть коэффициентов при неизвестных в области = 0)
             !(f(1)={0 : x<0.5; y*(sin(x)-x*cos(x)) : x>0.5} f(4)={0 : x<0.5; x : x>0.5} f(5)={y : x<0.5; 0 : x>0.5})
       !49 - y*sin(x) -(f(1)=y**2*cos(x)**2-y*sin(x); f(2)=x; f(3)=-y*cos(x); f(4)=-x*y)
integer(4) nmain
           !1 - main3 main3_2
		   !2 - main_p1 main_p1_2
		   !3 - main_p2
		   !4 - main_p4 main_p4_2
		   !5 - main_s1 
		   !6 - main_s2
		   !7 - main_p5
		   !8 - main_per1
		   !9 - main_per2
		   !10 - main_per5
		   !11 - main1
		   !12 - main2
		   !13 - main4
		   !14 - main5
		   !15 - main_h
		   !16 - main_konform
		   !17 - main_mc1
		   !18 - main_kuvbri,main_kuvbri5
		   !19 - main_kuvbri6
		   !20 - main_oss1
		   !21 - main_oss2
		   !22 - main6
       !23 - main_p6
       !24 - main_p7 main_p7_2
       !25 - main_p8


end

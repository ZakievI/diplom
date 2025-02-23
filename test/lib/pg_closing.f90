subroutine pg_closing_matrix_periodic(ibndl1,ibndl2,is_direct,dpsi1,dom)
!dec$ attributes dllexport :: pg_closing_matrix_periodic
!замыкаем получившуюся СЛАУ.
!периодические условия для квадратной ячейки
use pgmod
integer(4) i,j,k
integer(4) ibndl1 !номер первого участка текущей границы
integer(4) ibndl2 !номер второго участка текущей границы
logical is_direct  !прямая или обратная (шахматный порядок) периодичность 
integer(4) i1   !первая панель первого участка
integer(4) dir1 !=+1,-1 направление возрастания интекса первого участка
integer(4) i2   !первая панель второго участка
integer(4) dir2 !=+1,-1 направление возрастания интекса второго участка
real(8) dpsi1   !=0 для прямой периодичности (одно тело над другим)
                !=h для обратной периодичности (тела в шахматном порядке)
real(8) dom
real(8) b
real(8), allocatable :: m(:)
type(TBoundline), pointer :: bl1,bl2
type(col_info) ci
allocate(m(gs%m%nx))
bl1=>gsbnd%line(ibndl1)
bl2=>gsbnd%line(ibndl2)
if (is_direct) then
  i1=bl1%i_begin
  dir1=1
  i2=bl2%i_end
  dir2=-1
else
  i1=bl1%i_begin
  dir1=1
  i2=bl2%i_begin
  dir2=1
endif
call init_ind_psiom_all(gsarea,gsbnd,bl1,ci,gsarea,gsbnd,bl2,ci) 
do i=1,bl1%npanel
  do j=1,gsarea%nu
    k=j*2-1
    m=d0
	  if (is_direct) then
	    !psi_left=psi_right      !om_left=om_right 
      m(gsbnd%psiind(i1,k))=d1
      m(gsbnd%psiind(i2,k))=-d1
    else
	    !psi_left=h-psi_right    !om_left=dom-om_right 
      m(gsbnd%psiind(i1,k))=d1
      m(gsbnd%psiind(i2,k))=d1
    endif
    if (j==1) then
	    b=dpsi1
	  else
	    b=dom
    endif
    call add_closing_eq_to_area(gsarea%i,m,b)
    m=d0
	  if (is_direct) then
	    !psi'_left=-psi'_right   !om'_left=-om'_right
      m(gsbnd%psiind(i1,k+1))=d1
      m(gsbnd%psiind(i2,k+1))=d1
	  else
	    !psi'_left=psi'_right    !om'_left=om'_right
	    m(gsbnd%psiind(i1,k+1))=d1
      m(gsbnd%psiind(i2,k+1))=-d1
    endif
    call add_closing_eq_to_area(gsarea%i,m,d0)
  enddo
  i1=i1+dir1
  i2=i2+dir2
enddo
deallocate(m)
call deallocate_col_info(ci)
end

subroutine init_col_info_closing(a,ci,b1,bl1,k1,b2,bl2,k2)
use pgmod
Type(TArea) a
type(col_info), target :: ci
Type(TBound) b1,b2
type(TBoundline) bl1,bl2
integer(4), allocatable :: buff_c(:)
logical, allocatable :: ind(:)
integer(4) k1(:),k2(:)
integer(4) i
allocate(ind(gs%m%nx),buff_c(gs%m%nx))
ind=.false.
do i=1,ubound(k1,1)
  call init_ind_psiom(b1,bl1,ind,k1(i))
enddo
do i=1,ubound(k2,1)
  call init_ind_psiom(b2,bl2,ind,k2(i))
enddo
call init_col_info(a,ci,ind,buff_c,gs%m%nx,.false.)
a%m%ci=>ci
deallocate(ind,buff_c)
end

subroutine init_ind_psiom(b,bl,ind,k)
use pgmod
Type(TBound) b
type(TBoundline) bl
logical ind(gs%m%nx)
integer(4) k,i
do i=bl%i_begin,bl%i_end
  ind(b%psiind(i,k))=.true.
enddo
end

subroutine init_ind_psiom_all(a1,b1,bl1,ci1,a2,b2,bl2,ci2)
use pgmod
Type(TArea) a1,a2
type(col_info) ci1,ci2
Type(TBound) b1,b2
type(TBoundline) bl1,bl2
integer(4), allocatable :: k1(:),k2(:)
integer(4) i
allocate(k1(a1%umax),k2(a2%umax))
forall (i=1:a1%umax) k1(i)=i
forall (i=1:a2%umax) k2(i)=i
call init_col_info_closing(a1,ci1,b1,bl1,k1,b2,bl2,k2)
if (a1%i/=a2%i) call init_col_info_closing(a2,ci2,b2,bl2,k2,b1,bl1,k1)
deallocate(k1,k2)
end

subroutine pg_closing_matrix_darcidarci(ia1,ibnd1,ibndl1,ia2,ibnd2,ibndl2,kkr)
!dec$ attributes dllexport:: pg_closing_matrix_darcidarci
!замыкаем получившуюся СЛАУ
!цилиндр в пористой среде
use pgmod
integer(4) i,j,i0,k,i1
!первая область
integer(4) ia1    !номер области
integer(4) ibnd1  !номер границы
integer(4) ibndl1 !номер участка границы
!вторая область
integer(4) ia2    !номер области
integer(4) ibnd2  !номер границы
integer(4) ibndl2 !номер участка границы
real(8), allocatable :: m(:)
type(TArea), pointer :: a1,a2
type(TBound), pointer :: b1,b2
type(TBoundline), pointer :: bl1,bl2
real(8) kkr  !отношение проницаемостей k2/k1
type(col_info) ci1,ci2
allocate(m(gs%m%nx))
a1=>gs%a(ia1)
a2=>gs%a(ia2)
b1=>a1%bnd(ibnd1)
b2=>a2%bnd(ibnd2)
bl1=>b1%line(ibndl1)
bl2=>b2%line(ibndl2)
call init_col_info_closing(a1,ci1,b1,bl1,[1],b2,bl2,[1])
call init_col_info_closing(a2,ci2,b1,bl1,[2],b2,bl2,[2])
do i0=1,2
  j=bl1%i_begin
  i=bl2%i_end
  i1=ia1
  if (i0>1) i1=ia2
  do k=1,bl1%npanel
    m=d0
    if (i0==1) then
	  !ГУ psi_e=psi_i
      m(b1%psiind(j,1))=d1 !psi_e
      m(b2%psiind(i,1))=-d1 !psi_i
    elseif (i0==2) then
      !ГУ dpsi_e/dn+dpsi_i/dn=0
	    m(b1%psiind(j,2))=kkr !dpsi_e/dn
      m(b2%psiind(i,2))=d1 !dpsi_i/dn
    endif
    call add_closing_eq_to_area(i1,m,d0)
    j=j+1
    i=i-1
  enddo
enddo
deallocate(m)
call deallocate_col_info(ci1)
call deallocate_col_info(ci2)
end

subroutine pg_closing_matrix_darci(ia1,ibnd1,ibndl1,ia2,ibnd2,ibndl2,s_darci,alf,bet)
!dec$ attributes dllexport:: pg_closing_matrix_darci
!замыкаем получившуюся СЛАУ. ГУ Стокс-Дарси
!параметры с индексом 1 - стокс, 2 - Дарси
!omega_e=S(alf*u_e_tau-bet*u_i_tau)
use pgmod
!первая область
integer(4) ia1    !номер области
integer(4) ibnd1  !номер границы
integer(4) ibndl1 !номер участка границы
!вторая область
integer(4) ia2    !номер области
integer(4) ibnd2  !номер границы
integer(4) ibndl2 !номер участка границы
real(8) s_darci !параметр S для пористой (второй) области =1/sqrt(k)
real(8) alf,bet 
real(8) alf_s,bet_s
integer(4) i,j,k,i0,i1
real(8) s_darci2
real(8), allocatable :: m(:)
type(TArea), pointer :: a1,a2
type(TBound), pointer :: b1,b2
type(TBoundline), pointer :: bl1,bl2
type(col_info) ci1,ci2
allocate(m(gs%m%nx))
a1=>gs%a(ia1)
a2=>gs%a(ia2)
b1=>a1%bnd(ibnd1)
b2=>a2%bnd(ibnd2)
bl1=>b1%line(ibndl1)
bl2=>b2%line(ibndl2)
s_darci2=s_darci**2
alf_s=alf*s_darci
bet_s=bet*s_darci
call init_col_info_closing(a1,ci1,b1,bl1,[1,2,3],b2,bl2,[1,2])
call init_col_info_closing(a2,ci2,b1,bl1,[4],b2,bl2,[2])
do i0=1,3
  j=bl1%i_begin
  i=bl2%i_end
  i1=ia1
  if (i0>2) i1=ia2
  do k=1,bl1%npanel
    m=d0
    select case (i0)
    case (1)
      !ГУ psi_e=psi_i
      m(b1%psiind(j,1))=d1 !psi_e
      m(b2%psiind(i,1))=-d1 !psi_i
    case (2)
      !ГУ om_e=S(alf*v_e_tau-bet*v_i_tau)  аналог Биверса
	    m(b1%psiind(j,2))=alf_s !dpsi_e/dn
	    m(b1%psiind(j,3))=d1      !\Delta\psi_e
      m(b2%psiind(i,2))=bet_s !dpsi_i/dn
    case (3)
      !ГУ равенства давлений
	    m(b1%psiind(j,4))=d1 !d\Delta\psi_e/dn
      m(b2%psiind(i,2))=-s_darci2 !dpsi_i/dn
    endselect
    call add_closing_eq_to_area(i1,m,d0)
    j=j+1
    i=i-1
  enddo
enddo
deallocate(m)
call deallocate_col_info(ci1)
call deallocate_col_info(ci2)
end

subroutine pg_closing_matrix_bri(ia1,ibnd1,ibndl1,ia2,ibnd2,ibndl2,s_darci)
!dec$ attributes dllexport:: pg_closing_matrix_bri
!замыкаем получившуюся СЛАУ. ГУ Стокс-Бринкман
!параметры с индексом 1 - стокс, 2 - бринкман
use pgmod
!первая область
integer(4) ia1    !номер области
integer(4) ibnd1  !номер границы
integer(4) ibndl1 !номер участка границы
!вторая область
integer(4) ia2    !номер области
integer(4) ibnd2  !номер границы
integer(4) ibndl2 !номер участка границы
real(8) s_darci !параметр S для пористой (второй) области
integer(4) i,j,k,i0,i1
real(8) s_darci2
real(8), allocatable :: m(:)
type(TArea), pointer :: a1,a2
type(TBound), pointer :: b1,b2
type(TBoundline), pointer :: bl1,bl2
type(col_info) ci1,ci2
allocate(m(gs%m%nx))
a1=>gs%a(ia1)
a2=>gs%a(ia2)
b1=>a1%bnd(ibnd1)
b2=>a2%bnd(ibnd2)
bl1=>b1%line(ibndl1)
bl2=>b2%line(ibndl2)
s_darci2=s_darci**2
call init_col_info_closing(a1,ci1,b1,bl1,[1,2],b2,bl2,[1,2])
call init_col_info_closing(a2,ci2,b1,bl1,[3,4],b2,bl2,[2,4])
do i0=1,4
  j=bl1%i_begin
  i=bl2%i_end
  i1=ia1
  if (i0>2) i1=ia2
  do k=1,bl1%npanel
    m=d0
    select case (i0)
    case (1)
      !ГУ psi_e=psi_i
      m(b1%psiind(j,1))=d1 !psi_e
      m(b2%psiind(i,1))=-d1 !psi_i
    case (2)
      !ГУ dpsi_e/dn+dpsi_i/dn=0
	    m(b1%psiind(j,2))=d1 !dpsi_e/dn
      m(b2%psiind(i,2))=d1 !dpsi_i/dn
    case (3)
      !ГУ равенства давлений
	    m(b1%psiind(j,4))=d1 !d\Delta\psi_e/dn
      m(b2%psiind(i,2))=-s_darci2 !dpsi_i/dn
      m(b2%psiind(i,4))=d1 !d\Delta\psi_i/dn  
    case (4)
    !omega_e=omega_i
	    m(b1%psiind(j,3))=-d1 !\Delta\psi_e
      m(b2%psiind(i,3))=d1 !\Delta\psi_i
    end select
    call add_closing_eq_to_area(i1,m,d0)
    j=j+1
    i=i-1
  enddo
enddo
deallocate(m)
call deallocate_col_info(ci1)
call deallocate_col_info(ci2)
end

subroutine pg_closing_matrix_bri2(ia1,ibnd1,ibndl1,ia2,ibnd2,ibndl2,aa,bb,cc,dd,k_bri,mu_bri)
!dec$ attributes dllexport:: pg_closing_matrix_bri2
!замыкаем получившуюся СЛАУ. ГУ Стокс-Бринкман
!параметры с индексом 1 - стокс, 2 - бринкман
use pgmod
!первая область
integer(4) ia1    !номер области
integer(4) ibnd1  !номер границы
integer(4) ibndl1 !номер участка границы
!вторая область
integer(4) ia2    !номер области
integer(4) ibnd2  !номер границы
integer(4) ibndl2 !номер участка границы
real(8) k_bri !проницаемость для пористой (второй) области 
real(8) mu_bri !эф. вязкость для пористой (второй) области 
real(8) sqk !1/sqrt(k)
real(8) sqk2 !sqrt(k)
real(8) k_bri1 !1/k
real(8) aa,bb !коэффициенты в ГУ a*om_e-b*om_i=v_tau_e/sqrt(k)
real(8) cc,dd
integer(4) i,j,k,i0,i1
real(8), allocatable :: m(:)
type(TArea), pointer :: a1,a2
type(TBound), pointer :: b1,b2
type(TBoundline), pointer :: bl1,bl2
type(col_info) ci1,ci2
allocate(m(gs%m%nx))
a1=>gs%a(ia1)
a2=>gs%a(ia2)
b1=>a1%bnd(ibnd1)
b2=>a2%bnd(ibnd2)
bl1=>b1%line(ibndl1)
bl2=>b2%line(ibndl2)
k_bri1=d1/k_bri
sqk=dsqrt(k_bri1)
sqk2=dsqrt(k_bri)
!call init_col_info_closing(a1,ci1,b1,bl1,[1,2],b2,bl2,[1,2])
!call init_col_info_closing(a2,ci2,b1,bl1,[2,3,4],b2,bl2,[2,3,4])

!call init_col_info_closing(a1,ci1,b1,bl1,[1,2,3,4],b2,bl2,[1])
!call init_col_info_closing(a2,ci2,b1,bl1,[3,4],b2,bl2,[2,3,4])

call init_col_info_closing(a1,ci1,b1,bl1,[1,2,3,4],b2,bl2,[1,2,3,4])
call init_col_info_closing(a2,ci2,b1,bl1,[1,2,3,4],b2,bl2,[1,2,3,4])
do i0=1,4
  j=bl1%i_begin
  i=bl2%i_end
  i1=ia1
  if (i0>2) i1=ia2
  do k=1,bl1%npanel
    m=d0
    select case (i0)
    case (1)
      !ГУ psi_e=psi_i
      m(b1%psiind(j,1))=d1 !psi_e
      m(b2%psiind(i,1))=-d1 !psi_i
    case (2)
      !ГУ dpsi_e/dn+dpsi_i/dn=0
	    !m(b1%psiind(j,2))=d1 !dpsi_e/dn
      !m(b2%psiind(i,2))=d1 !dpsi_i/dn
      
      !omega_e=aa/sqrt(k)*v_tau_e (переносим вправо)
      !m(b1%psiind(j,3))=d1 !Delta\psi_e
      !m(b1%psiind(j,2))=aa*sqk !dpsi_e/dn
      
      !omega_e=aa/sqrt(k)*v_tau_e-bb*sqrt(k)*dOmega_e/dn (переносим вправо)
      !m(b1%psiind(j,3))=d1 !Delta\psi_e
      !m(b1%psiind(j,2))=aa*sqk !dpsi_e/dn
      !m(b1%psiind(j,4))=bb*sqk2 !d\Delta\psi_e/dn
      
      !omega_e-b*sqrt(mub)*omega_i=(aa*v_tau_e-bb*v_tau_i)/sqrt(k) (переносим вправо)
      !m(b1%psiind(j,3))=d1 !Delta\psi_e
      !m(b1%psiind(j,2))=aa*sqk !dpsi_e/dn
      !m(b2%psiind(i,3))=-bb*dsqrt(mu_bri) !Delta\psi_i
      !m(b2%psiind(i,2))=bb*sqk !dpsi_i/dn
      
      !omega_i=cc*omega_e+dd/sqrt(k)*v_tau_i-dd*sqrt(mub)*Omega_i (переносим вправо)
      !m(b1%psiind(j,3))=-aa !Delta\psi_e
      !m(b2%psiind(i,3))=d1+bb*dsqrt(mu_bri) !Delta\psi_i
      !m(b2%psiind(i,2))=-bb*sqk !dpsi_i/dn
      
      !vtau_e=a*sqrt(k)*om_e+b*k*dom_e/dn (переносим влево)
      !m(b1%psiind(j,2))=d1 !dpsi_e/dn
      !m(b1%psiind(j,3))=aa*sqk2 !Delta\psi_e
      !m(b1%psiind(j,4))=bb*k_bri !dDelta\psi_e/dn
      
      !vtau_e=a*sqrt(k)*om_e+b*(vtau_i-sqrt(k*mu)*om_i) (переносим влево)
      !m(b1%psiind(j,2))=d1 !dpsi_e/dn
      !m(b1%psiind(j,3))=aa*sqk2 !Delta\psi_e
      !m(b2%psiind(i,2))=bb !dpsi_i/dn
      !m(b2%psiind(i,3))=-bb*dsqrt(k_bri*mu_bri) !Delta\psi_i
      
      !!a*vtau_e=sqrt(k)*(b*om_e-om_i) (переносим влево)
      !m(b1%psiind(j,2))=aa !dpsi_e/dn
      !m(b1%psiind(j,3))=bb*sqk2 !Delta\psi_e
      !m(b2%psiind(i,3))=-sqk2 !Delta\psi_i
      
      !a*vtau_e=sqrt(k)*(om_i-b*om_e) (переносим влево)
      m(b1%psiind(j,2))=aa !dpsi_e/dn
      m(b1%psiind(j,3))=-bb*sqk2 !Delta\psi_e
      m(b2%psiind(i,3))=sqk2 !Delta\psi_i
      
    case (3)
      !ГУ равенства давлений
	    m(b1%psiind(j,4))=d1 !d\Delta\psi_e/dn 
      m(b2%psiind(i,2))=-k_bri1 !dpsi_i/dn
      m(b2%psiind(i,4))=mu_bri !d\Delta\psi_i/dn  
    case (4)
      !a*om_e-b*om_i=v_tau_e/sqrt(k)
	    !m(b1%psiind(j,3))=-aa !\Delta\psi_e
      !m(b2%psiind(i,3))=bb !\Delta\psi_i
      !m(b1%psiind(j,2))=-sqk !dpsi_e/dn
      
      !a*v_tau_e-b*k*dom_e/dn=sqrt(k)*om_e
	    !m(b1%psiind(j,2))=aa !dpsi_e/dn
      !m(b1%psiind(j,4))=bb*k_bri !d\Delta\psi_e/dn
      !m(b1%psiind(j,3))=sqk !\Delta\psi_e
      
      !omega_i=bb*omega_e
      !m(b1%psiind(j,3))=bb !Delta\psi_e
      !m(b2%psiind(i,3))=-d1 !Delta\psi_i
      
      !omega_i=cc*omega_e+dd*sqrt(k)*dOmega_e/dn
      !m(b1%psiind(j,3))=cc !Delta\psi_e
      !m(b2%psiind(i,3))=-d1 !Delta\psi_i
      !m(b1%psiind(j,4))=dd*sqk2 !d\Delta\psi_e/dn
      
      !omega_i=cc*omega_e+dd/sqrt(k)*v_tau_i-dd*sqrt(mub)*Omega_i (переносим вправо)
      !m(b1%psiind(j,3))=-cc !Delta\psi_e
      !m(b2%psiind(i,3))=d1+dd*dsqrt(mu_bri) !Delta\psi_i
      !m(b2%psiind(i,2))=-dd*sqk !dpsi_i/dn
      
      !vtau_i=cc*vtau_e (переносим вправо)
      !m(b2%psiind(i,2))=d1 !dpsi_i/dn
      !m(b1%psiind(j,2))=cc !dpsi_e/dn
      
      !vtau_i=cc*vtau_e+dd*(vtau_i-sqrt(k*mub)*omega_i) (переносим вправо)
      !m(b2%psiind(i,2))=d1-dd !dpsi_i/dn
      !m(b1%psiind(j,2))=cc !dpsi_e/dn
      !m(b2%psiind(i,3))=dd*dsqrt(k_bri*mu_bri) !Delta\psi_i
      
      !vtau_i=cc*vtau_e+dd*k*dom_e/dn (переносим вправо)
      !m(b2%psiind(i,2))=d1 !dpsi_i/dn
      !m(b1%psiind(j,2))=cc !dpsi_e/dn
      !m(b1%psiind(j,4))=-dd*k_bri !dDelta\psi_e/dn

      !om_i=sqrt(k*mub)*domega_i/dn (переносим вправо)
      !m(b2%psiind(i,3))=d1 !Delta\psi_i
      !m(b2%psiind(i,4))=-dsqrt(k_bri*mu_bri) !dDelta\psi_i/dn
      
      !cc*vtau_e-vtau_i=dd*sqrt(k)*omega_e (переносим влево)
      m(b2%psiind(i,2))=d1 !dpsi_i/dn
      m(b1%psiind(j,2))=cc !dpsi_e/dn
      m(b1%psiind(j,3))=dd*sqk2 !Delta\psi_e
      
    end select
    call add_closing_eq_to_area(i1,m,d0)
    j=j+1
    i=i-1
  enddo
enddo
deallocate(m)
call deallocate_col_info(ci1)
call deallocate_col_info(ci2)
end

subroutine pg_closing_matrix_bribri(ia1,ibnd1,ibndl1,ia2,ibnd2,ibndl2,s_darci1,s_darci2)
!dec$ attributes dllexport:: pg_closing_matrix_bribri
!замыкаем получившуюся СЛАУ
!ГУ бринкман-бринкман
use pgmod
!первая область
integer(4) ia1    !номер области
integer(4) ibnd1  !номер границы
integer(4) ibndl1 !номер участка границы
!вторая область
integer(4) ia2    !номер области
integer(4) ibnd2  !номер границы
integer(4) ibndl2 !номер участка границы
real(8) s_darci1 !параметр S для первой области
real(8) s_darci2 !параметр S для второй области
integer(4) i,j,k,i0,i1
real(8) s_darci12,s_darci22
real(8), allocatable :: m(:)
type(TArea), pointer :: a1,a2
type(TBound), pointer :: b1,b2
type(TBoundline), pointer :: bl1,bl2
type(col_info) ci1,ci2
allocate(m(gs%m%nx))
a1=>gs%a(ia1)
a2=>gs%a(ia2)
b1=>a1%bnd(ibnd1)
b2=>a2%bnd(ibnd2)
bl1=>b1%line(ibndl1)
bl2=>b2%line(ibndl2)
s_darci12=s_darci1**2
s_darci22=s_darci2**2
call init_col_info_closing(a1,ci1,b1,bl1,[1,2],b2,bl2,[1,2])
call init_col_info_closing(a2,ci2,b1,bl1,[2,3,4],b2,bl2,[2,3,4])
do i0=1,4
  j=bl1%i_begin
  i=bl2%i_end
  i1=ia1
  if (i0>2) i1=ia2
  do k=1,bl1%npanel
    m=d0
    select case (i0)
    case (1)
      !ГУ psi_e=psi_i
      m(b1%psiind(j,1))=d1 !psi_e
      m(b2%psiind(i,1))=-d1 !psi_i
    case (2)
      !ГУ dpsi_e/dn+dpsi_i/dn=0
	    m(b1%psiind(j,2))=d1 !dpsi_e/dn
      m(b2%psiind(i,2))=d1 !dpsi_i/dn
    case (3)
      !ГУ равенства давлений
	    m(b1%psiind(j,2))=-s_darci12 !dpsi_e/dn
      m(b1%psiind(j,4))=d1 !d\Delta\psi_e/dn
      m(b2%psiind(i,2))=-s_darci22 !dpsi_i/dn
      m(b2%psiind(i,4))=d1 !d\Delta\psi_i/dn
    case (4)
      !omega_e=omega_i
	    m(b1%psiind(j,3))=-d1 !\Delta\psi_e
      m(b2%psiind(i,3))=d1 !\Delta\psi_i
    end select
    call add_closing_eq_to_area(i1,m,d0)
    j=j+1
    i=i-1
  enddo
enddo
deallocate(m)
call deallocate_col_info(ci1)
call deallocate_col_info(ci2)
end

subroutine pg_add_closing_eq_to_area_with_fict(ia,eq,b)
!dec$ attributes dllexport:: pg_add_closing_eq_to_area_with_fict
!добавлнение замыкающего уравнения к блоку области 
!уравнение может содержать фиктивные переменные
use pgmod
integer(4) ia !индекс области
real(8) eq(gs%m%nx_all),b,b1 !уравнение и свободный член
b1=b !входное b может быть константой, но она будет меняться в update_fict_vars
call update_fict_vars(eq,b1,ia)
call pg_add_closing_eq_to_area(ia,eq,b1)
end

subroutine pg_add_closing_eq_to_area(ia,eq,b)
!dec$ attributes dllexport:: pg_add_closing_eq_to_area
!добавлнение замыкающего уравнения к блоку области
use pgmod
integer(4) ia !индекс области
real(8) eq(gs%m%nx),b !уравнение и свободный член
integer(4), allocatable :: buff_c(:)
logical, allocatable :: ind(:)
type(col_info), target :: ci
type(TArea), pointer :: a
if (ia>0) then
  allocate(ind(gs%m%nx),buff_c(gs%m%nx))
  a=>gs%a(ia)
  ind=.false.
  if (a%m%ix>1) ind(1:a%m%ix-1)=.true.
  if (a%m%ix2<gs%m%nx) ind(a%m%ix2+1:gs%m%nx)=.true.
  call init_col_info(a,ci,ind,buff_c,gs%m%nx,.false.)
  a%m%ci=>ci
  call add_closing_eq_to_area(ia,eq,b)
  deallocate(ind,buff_c)
  call deallocate_col_info(ci)
else
  call gs_print_stop("Error pg_add_closing_eq_to_area!!!")
endif
end

subroutine add_closing_eq_to_area(ia,eq,b)
!добавлнение замыкающего уравнения к блоку области
use pgmod
integer(4) ia !индекс области
real(8) eq(gs%m%nx),b !уравнение и свободный член
type(matrix_mb), pointer :: m
integer(4) k
m=>gs%a(ia)%m
m%nnt=m%nnt+1
if (m%nnt>m%nu_all) call gs_print_stop("Error add_closing_eq_to_area!!!")
k=m%nnt+m%iu-1
call add_equation_to_matrix(k,eq,b,ia)
end

!!!нужно указывать область, к которой добавляется замыкающее уравнение для корректной работы add_equation_to_matrix

!subroutine pg_add_closing_eq(eq,b)
!!dec$ attributes dllexport:: pg_add_closing_eq
!!добавлнение замыкающего уравнения
!use pgmod
!real(8) eq(gs%m%nx),b !уравнение и свободный член
!call pg_add_closing_eq_to_area(0,eq,b)
!end

!subroutine pg_replace_eq(eq,b,neq)
!!dec$ attributes dllexport:: pg_replace_eq
!!замена уравнения
!use pgmod
!real(8) eq(gs%m%nx),b !уравнение и свободный член
!integer(4) neq !номер заменяемого уравнения
!call add_equation_to_matrix(neq,eq,b)
!end

subroutine pg_closing_psiconst(ibnd)
!dec$ attributes dllexport:: pg_closing_psiconst
!замыкаем СЛАУ задачи для жидкой муфты
use pgmod
integer(4) ibnd !номер границы,для которой записывается условие
integer(4) abnd(1),abndl(1)
abnd=ibnd
abndl=0
call pg_closing_psiconst_arr(abnd,abndl,1)
end

subroutine pg_closing_psiconst_line(ibnd,ibndl)
!dec$ attributes dllexport:: pg_closing_psiconst_line
!замыкаем СЛАУ задачи для жидкой муфты (один участок)
use pgmod
integer(4) ibnd !номер границы,для которой записывается условие
integer(4) ibndl !участок границы,для которой записывается условие
integer(4) abnd(1),abndl(1)
abnd=ibnd
abndl=ibndl
call pg_closing_psiconst_arr(abnd,abndl,1)
end

subroutine pg_closing_psiconst2(ibnd,i1,i2)
!dec$ attributes dllexport:: pg_closing_psiconst2
!замыкаем СЛАУ задачи для жидкой муфты
!два участка одной границы (обтекаемая поверхность разрезана границей на две части)
use pgmod
integer(4) ibnd !номер границы
integer(4) i1,i2 !номера участков границы, для которой записывается условие
integer(4) abnd(2),abndl(2)
abnd=ibnd
abndl(1)=i1
abndl(2)=i2
call pg_closing_psiconst_arr(abnd,abndl,2)
end

subroutine pg_closing_psiconst_arr(ibnd,ibndl,n)
!dec$ attributes dllexport:: pg_closing_psiconst_arr
!замыкаем СЛАУ задачи для жидкой муфты
!для разных участков разных границ
!обтекаемая поверхность разрезана границами на несколько частей
!если ibndl(j)=0, то берется вся граница
use pgmod
integer(4) n !число участков/границ
integer(4) ibnd(n) !номера границы
integer(4) ibndl(n) !номера участков границы, для которой записывается условие
call pg_closing_psiconst_arr_dp(ibnd,ibndl,n,d0)
end

subroutine pg_closing_psiconst_arr_dp(ibnd,ibndl,n,dp)
!dec$ attributes dllexport:: pg_closing_psiconst_arr
!замыкаем СЛАУ задачи для жидкой муфты
!для разных участков разных границ
!обтекаемая поверхность разрезана границами на несколько частей
!если ibndl(j)=0, то берется вся граница
!случай заданного перепада давления
use pgmod
integer(4) n !число участков/границ
integer(4) ibnd(n) !номера границы
integer(4) ibndl(n) !номера участков границы, для которой записывается условие
integer(4) ia(n) !номера областей
real(8) dp   !перепад давления
ia=gsarea%i
call pg_closing_psiconst_arr2_dp(ia,ibnd,ibndl,n,gsarea%i,dp)
end

subroutine pg_closing_psiconst_arr2(ia,ibnd,ibndl,n,ia_add)
!dec$ attributes dllexport:: pg_closing_psiconst_arr2
!замыкаем СЛАУ задачи для жидкой муфты
!для разных участков разных границ и разных областей
!обтекаемая поверхность разрезана границами на несколько частей
!если ibndl(j)=0, то берется вся граница
use pgmod
integer(4) n !число участков/границ
integer(4) ia(n) !номера областей
integer(4) ibnd(n) !номера границы
integer(4) ibndl(n) !номера участков границы, для которой записывается условие
integer(4) ia_add !номер области, к которой добавляется уравнение
call pg_closing_psiconst_arr2_dp(ia,ibnd,ibndl,n,ia_add,d0)
end

subroutine pg_closing_psiconst_arr2_dp(ia,ibnd,ibndl,n,ia_add,dp)
!dec$ attributes dllexport:: pg_closing_psiconst_arr2
!замыкаем СЛАУ задачи для жидкой муфты
!для разных участков разных границ и разных областей
!обтекаемая поверхность разрезана границами на несколько частей
!если ibndl(j)=0, то берется вся граница
!случай заданного перепада давления
use pgmod
integer(4) n !число участков/границ
integer(4) ia(n) !номера областей
integer(4) ibnd(n) !номера границы
integer(4) ibndl(n) !номера участков границы, для которой записывается условие
integer(4) ia_add !номер области, к которой добавляется уравнение
real(8) dp   !перепад давления
integer(4) i,j,i1,i2,ib,ie
real(8), allocatable :: eq(:)
type(TArea), pointer :: a
type(TBound), pointer :: b
type(TBoundline), pointer :: bl
type(col_info), target :: ci
integer(4), allocatable :: buff_c(:)
logical, allocatable :: ind(:)
integer(4) type_eq
real(8) s2
!проверяем тип уравнения в области
type_eq=-1
s2=-d1
do j=1,n
  a=>gs%a(ia(j))
  if (type_eq<0) then
    type_eq=a%type_eq(1)
  else
    if (type_eq/=a%type_eq(1)) call gs_print_stop("Error 1 pg_closing_psiconst_arr2_dp!")
  endif
  if (type_eq==15) then
    if (s2<0) then
      s2=a%const%k_helm
    else
      if (s2/=a%const%k_helm) call gs_print_stop("Error S2 pg_closing_psiconst_arr2_dp!")
    endif
  endif
enddo
select case (type_eq)
case (1,3,21,24,25)
case (15)
  s2=s2**2
case default
  call gs_print_stop("Error 2 pg_closing_psiconst_arr2_dp!")
endselect
!составляем уравнение 
allocate(eq(gs%m%nx))
allocate(ind(gs%m%nx),buff_c(gs%m%nx))
ind=.false. 
eq=d0
do j=1,n
  a=>gs%a(ia(j))
  b=>a%bnd(ibnd(j))
  ib=ibndl(j)
  if (ib==0) then
    i1=1
	  i2=b%npanel
    ib=1
    ie=b%nline
  else
    bl=>b%line(ib)
	  i1=bl%i_begin
	  i2=bl%i_end
    ie=ib
  endif
  select case (type_eq)
  case (1)
    do i=i1,i2
      eq(b%psiind(i,2))=b%l(i)
    enddo
  case (3,21,24,25) !!!для Навье-Стокса предполагаем прилипание (v_n=v_tau=0)
    do i=i1,i2
      eq(b%psiind(i,4))=b%l(i)
    enddo
  case (15)
    do i=i1,i2
      !eq(b%psiind(i,2))=b%l(i)*s2  !!!предполагается условие прилипания (v_tau=d\psi/dn=0)
      eq(b%psiind(i,4))=b%l(i)
    enddo
  endselect
  do i=ib,ie
    bl=>b%line(i)
    !if (type_eq==1.or.type_eq==15) call init_ind_psiom(b,bl,ind,2) !!!!см выше
    select case (type_eq)
    case (1)
      call init_ind_psiom(b,bl,ind,2)
    case (3,15,21,24,25)
      call init_ind_psiom(b,bl,ind,4)
    endselect
  enddo
enddo
a=>gs%a(ia_add)
call init_col_info(a,ci,ind,buff_c,gs%m%nx,.false.)
a%m%ci=>ci
call add_closing_eq_to_area(a%i,eq,dp)
deallocate(eq)
deallocate(ind,buff_c)
call deallocate_col_info(ci)
end

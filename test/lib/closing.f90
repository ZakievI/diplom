subroutine closing_matrix2_per(periodic)
!периодические условия для квадратной ячейки
use pgmod2
integer(4) periodic
call pg_bind_domain(1)
call pg_bind_bound(1)
if (periodic==1) then
  !прямая периодичность
  call pg_closing_matrix_periodic(2,4,.true.,d0,d0)
else
  !обратная периодичность
  call pg_closing_matrix_periodic(2,4,.false.,h,d0)
endif
end

subroutine closing_matrix5
!замыкаем получившуюся СЛАУ
!пористый цилиндр в течении стокса с условием Биверса
use pgmod2
integer(4) i,j,j1,j2,np1,np2,pg_get_int,nx
real(8) dg2,kkk
real(8), allocatable :: eq(:)
integer(4), allocatable :: b1psiind(:,:), b2psiind(:,:)
type(col_info) ci1,ci2
call init_ind_psiom_all(gs%a(1),gs%a(1)%bnd(1),gs%a(1)%bnd(1)%line(4),ci1,gs%a(2),gs%a(2)%bnd(1),gs%a(2)%bnd(1)%line(1),ci2) 
call pg_bind_domain(1)
call pg_bind_bound(1)
np1=pg_get_int(1,1)
allocate(b1psiind(np1,4))
do i=1,4
  call pg_get_array_int2(1, 1, i, b1psiind(:,i), np1)
enddo
call pg_bind_domain(2)
call pg_bind_bound(1)
np2=pg_get_int(1,1)
allocate(b2psiind(np2,2))
do i=1,2
  call pg_get_array_int2(1, 1, i, b2psiind(:,i), np2)
enddo
nx=pg_get_int(3,1)
allocate(eq(nx))
!ГУ psi1=psi2
j=np1+1
do i=1,nj_
  j=j-1
  eq=d0
  eq(b1psiind(j,1))=d1 !psi1
  eq(b2psiind(i,1))=-d1 !psi2
  call pg_add_closing_eq_to_area(1,eq,d0)
enddo
!ГУ равенства давлений
s_darci2=s_darci**2
j=np1+1
do i=1,nj_
  j=j-1
  eq=d0
  eq(b1psiind(j,4))=d1 !domega_e/dn
  eq(b2psiind(i,2))=-s_darci2 !dpsi_i/dn
  call pg_add_closing_eq_to_area(1,eq,d0)
enddo
!ГУ условие Биверса
dg2=(pi/nj_)**(-2)
j=np1+1
kkk=d1 
do i=1,nj_
  j=j-1
  eq=d0
  eq(b1psiind(j,3))=d1 !\Delta\psi_e
  eq(b1psiind(j,2))=d1+alf_darci*s_darci !dpsi_e/dn
  eq(b2psiind(i,2))=alf_darci*s_darci !dpsi_i/dn
  !-d2psi_e/dtt2
  j1=j-1
  j2=j+1
  if (i==1) then
    eq(b1psiind(j1,1))=-dg2*kkk
	  eq(b1psiind(j,1))=3.0d0*dg2*kkk
  elseif (i==nj_) then
	  eq(b1psiind(j,1))=3.0d0*dg2*kkk
	  eq(b1psiind(j2,1))=-dg2*kkk
  else
    eq(b1psiind(j1,1))=-dg2*kkk
	  eq(b1psiind(j,1))=d2*dg2*kkk
	  eq(b1psiind(j2,1))=-dg2*kkk
  endif
  call pg_add_closing_eq_to_area(2,eq,d0)
enddo
deallocate(eq,b1psiind,b2psiind)
call deallocate_col_info(ci1)
call deallocate_col_info(ci2)
end

subroutine closing_matrix5_3
!замыкание СЛАУ Стокс-Дарси с условием Биверса через omega (мое ГУ)
use pgmod2
call pg_closing_matrix_darci(1,1,gs%a(1)%bnd(1)%nline,2,1,1,s_darci,alf_darci,bet_darci)
end

subroutine closing_matrix_bri
!замыкание СЛАУ Стокс-Бринкман
use pgmod2
!call pg_closing_matrix_bri(1,1,gs%a(1)%bnd(1)%nline,2,1,1,s_darci)
call pg_closing_matrix_bri2(1,1,gs%a(1)%bnd(1)%nline,2,1,1,a_bri,b_bri,c_bri,d_bri,k_bri,mu_bri)
end

subroutine closing_matrix_bribri
!замыкание СЛАУ Стокс-Бринкман
use pgmod2
call pg_closing_matrix_bribri(1,1,gs%a(1)%bnd(1)%nline,2,1,1,s_darci,s_darci_2)
end



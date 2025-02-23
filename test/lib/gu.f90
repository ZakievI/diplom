subroutine init_gu1
!������������� ��������� ������� ��� ������ �������������� ���������
!������� � ��������������� ������
use pgmod2
real(8) k,get_funa,gg,rr,get_funa_oss,vx,vy
real(8) bndval(nmax),xc(nmax),yc(nmax)
integer(4) i,j,n,pg_get_int
call pg_allocate_bound_gu
k=(h**2+d1)/(h**2-d1)
do i=1,4
  call pg_bind_boundline(i)
  if (i==1.or.i==3.or.i==4) then
    call pg_init_boundline_gu_val_const(1,1,d0)
  else
    call pg_get_array_real(2,3,xc,nmax)
    call pg_get_array_real(2,4,yc,nmax)
	  n=pg_get_int(2,1)
    do j=1,n
  	  gg=datan2(yc(j),xc(j))
      rr=dsqrt(yc(j)*yc(j)+xc(j)*xc(j))
	  if (nmain==11) then
	    bndval(j)=get_funa(rr,gg) 
      elseif (nmain==21) then
	    bndval(j)=get_funa_oss(rr,gg,vx,vy) 
	  endif
	enddo
	call pg_init_boundline_gu_val(1,1,bndval)
  endif
enddo
end

subroutine init_gu2(square,periodic,need_gu_circle,for_konform,slip_gu_circle)
!������������� ��������� ������� ��� ������ ������� ������ � ������ ��������
use pgmod2
logical square,periodic,need_gu_circle,for_konform,slip_gu_circle
real(8) yc(nmax)
integer(4) n,k,pg_get_int
call pg_allocate_bound_gu
!������ ����� ���������
call pg_bind_boundline(1)
call pg_init_boundline_gu_val_const(1,1,d0)
call pg_init_boundline_gu_val_const(2,3,d0)
if (square) then
  !���� ������� ������ �������
  call pg_bind_boundline(2)
  if (periodic) then
    call pg_init_boundline_gu_empty	
  else
    call pg_get_array_real(2,4,yc,nmax)
	call pg_init_boundline_gu_val(1,1,yc)
    call pg_init_boundline_gu_val_const(2,3,d0)
  endif
  !���� ������� ����
  call pg_bind_boundline(3)
  call pg_init_boundline_gu_val_const(1,1,h)
  call pg_init_boundline_gu_val_const(2,3,d0)
  !���� ������� ����� �������
  call pg_bind_boundline(4)
  if (periodic) then
    call pg_init_boundline_gu_empty	
  else
    call pg_get_array_real(2,4,yc,nmax)
	call pg_init_boundline_gu_val(1,1,yc)
    call pg_init_boundline_gu_val_const(2,3,d0)
  endif
  k=6
else
  !���� ����
  call pg_bind_boundline(2)
  call pg_get_array_real(2,4,yc,nmax)
  if (for_konform) then
    n=pg_get_int(2,1)
	yc(1:n)=h*dsin(yc(1:n))
  endif
  call pg_init_boundline_gu_val(1,1,yc)
  call pg_init_boundline_gu_val_const(2,3,d0)
  k=4
endif
!����� ����� ���������
call pg_bind_boundline(k-1)
call pg_init_boundline_gu_val_const(1,1,d0)
call pg_init_boundline_gu_val_const(2,3,d0)
!��� ������� ��������� ������ ����� ����������
call pg_allocate_boundline_gu_add
call pg_init_boundline_gu_symmetr(3,2,.false.,1,1,1,d0,d1)
call pg_init_boundline_gu_symmetr(4,4,.false.,1,1,1,d0,d1)
!����� ����
call pg_bind_boundline(k)
if (need_gu_circle) then
  call pg_init_boundline_gu_val_const(1,1,d0)
  if (slip_gu_circle) then
    call pg_init_boundline_gu_slip(2,slip_k)
  else
    call pg_init_boundline_gu_val_const(2,2,d0)
  endif
else
  call pg_init_boundline_gu_empty	
endif
end

subroutine init_gu3
!������������� ��������� ������� ��� ������ 2-�� ������� � ��������
use pgmod2
               !��� ���������   
               !1 - ������������� �������
               !4 - ��������� ��������
               !5 - ���������� ��������� �����������
               !6 - ������������ ��������� �����������
			   !16,17 - ��������� ��������
integer(4) i
call pg_allocate_bound_gu
do i=1,gsbnd%nline
  call pg_bind_boundline(i)
  call init_gu3_
enddo
end

subroutine init_gu3_
!������������� ��������� ������� ��� ������ 2-�� ������� � ��������
use pgmod2
real(8) get_funb,get_funb1
integer(4) typegu,bndu,j,n,pg_get_int
real(8) bndval(nmax),xc(nmax),yc(nmax)
typegu=1 !1,2
call pg_get_array_real(2,3,xc,nmax)
call pg_get_array_real(2,4,yc,nmax)
n=pg_get_int(2,1)
do j=1,n
  if (typegu==1) then
    bndval(j)=get_funb(xc(j),yc(j))
  else
    if (gsbndl%i==1) then
      bndval(j)=get_funb(xc(j),yc(j))
    else
      bndval(j)=get_funb1(xc(j),yc(j),gsbndl%i)
    endif
  endif
enddo
bndu=1
if (typegu==2.and.gsbndl%i>1) bndu=2
call pg_init_boundline_gu_val(1,bndu,bndval)
end

subroutine init_gu3_1
!������������� ��������� ������� ��� ������ 2-�� ������� � ��������
!� ���������� �� ����������
use pgmod2
               !��� ���������   
               !1 - ������������� �������
               !4 - ��������� ��������
               !5 - ���������� ��������� �����������
               !6 - ������������ ��������� �����������
			   !16,17 - ��������� ��������
integer(4) i,k,ix,iy
k=0
do ix=1,nsx
  do iy=1,nsy
    k=k+1
    call pg_bind_domain(k)
    call pg_bind_bound(1)
    call pg_allocate_bound_gu
    do i=1,gsbnd%nline
      call pg_bind_boundline(i)
      if (i==1.and.iy>1) then  !���
        call pg_init_boundline_gu_subdomain_continuity(k,1,1,k-1,1,3)
      elseif (i==2.and.ix<nsx) then  !�����
        call pg_init_boundline_gu_empty
      elseif (i==3.and.iy<nsy) then  !����
        call pg_init_boundline_gu_empty
      elseif (i==4.and.ix>1) then  !����
        call pg_init_boundline_gu_subdomain_continuity(k,1,4,k-nsy,1,2)
      else
        call init_gu3_
      endif
    enddo
  enddo
enddo
end

subroutine init_gu3_3
!������������� ��������� ������� ��� ������ 2-�� ������� � ��������
!� ���������� �� ����������
!������������ pg_init_gu_multiarea
use pgmod2
external set_gu_internal3_3,init_gu3_
call pg_init_subdmain_connectivity(.true.)
call pg_init_gu_multiarea(set_gu_internal3_3,init_gu3_)
end

subroutine set_gu_internal3_3(ia1,ib1,ibl1,ia2,ib2,ibl2)
use pgmod2
integer(4) ia1,ib1,ibl1,ia2,ib2,ibl2
!call pg_init_boundline_gu_subdomain_continuity(ia1,ib1,ibl1,ia2,ib2,ibl2)
!call pg_init_boundline_gu_subdomain_continuity(ia2,ib2,ibl2,ia1,ib1,ibl1)
call pg_init_boundline_gu_subdomain_continuity_with_direction(ia1,ib1,ibl1,ia2,ib2,ibl2,d0,d1)
end

subroutine init_gu3_2
!������������� ��������� ������� ��� ������ 4-�� ������� � �������� (psi, om)
use pgmod2
integer(4) i
call pg_allocate_bound_gu
do i=1,gsbnd%nline
  call pg_bind_boundline(i)
  call init_gu3_2_
enddo
end

subroutine init_gu3_2_
!������������� ��������� ������� ��� ������ 4-�� ������� � �������� (psi, om)
use pgmod2
integer(4) j,n,pg_get_int
real(8) get_fund,get_fund2
real(8) bndval(nmax),xc(nmax),yc(nmax)
call pg_get_array_real(2,3,xc,nmax)
call pg_get_array_real(2,4,yc,nmax)
n=pg_get_int(2,1)
do j=1,n
  bndval(j)=get_fund(xc(j),yc(j))
enddo
call pg_init_boundline_gu_val(1,1,bndval)
do j=1,n
  bndval(j)=get_fund2(xc(j),yc(j))
enddo
call pg_init_boundline_gu_val(2,3,bndval)
end

subroutine init_gu3_4
!������������� ��������� ������� ��� ������ 4-�� ������� � ��������
!� ���������� �� ����������
!������������ pg_init_gu_multiarea
use pgmod2
external set_gu_internal3_3,init_gu3_2_
call pg_init_subdmain_connectivity(.true.)
call pg_init_gu_multiarea(set_gu_internal3_3,init_gu3_2_)
end

subroutine init_gu4_1(square_cell)
!������������� ��������� ������� ��� ������ ��������� �������� � �������� ����� - ��� ������� �������� �����
use pgmod2
logical square_cell
integer(4) i,n,pg_get_int,j
real(8) get_func
real(8) bndval(nmax),xc(nmax),yc(nmax)
call pg_allocate_bound_gu
!������ ����� ���������
call pg_bind_boundline(1)
call pg_init_boundline_gu_val_const(1,1,d0)
!���� ����
if (square_cell) then
  do j=2,4
    call pg_bind_boundline(j)
    call pg_get_array_real(2,4,yc,nmax)
    call pg_init_boundline_gu_val(1,1,yc)
  enddo
  j=5
else
  call pg_bind_boundline(2)
  call pg_get_array_real(2,3,xc,nmax)
  call pg_get_array_real(2,4,yc,nmax)
  n=pg_get_int(2,1)
  do i=1,n
    bndval(i)=get_func(dcmplx(xc(i),yc(i))) 
  enddo
  call pg_init_boundline_gu_val(1,1,bndval)
  j=3
endif
!����� ����� ���������
call pg_bind_boundline(j)
call pg_init_boundline_gu_val_const(1,1,d0)
!����� ����
call pg_bind_boundline(j+1)
call pg_init_boundline_gu_empty	
end

subroutine init_gu4_2
!������������� ��������� ������� ��� ������ ��������� �������� � �������� ����� - ��� ���������� �������� �����
use pgmod2
call pg_allocate_bound_gu
!����
call pg_bind_boundline(1)
call pg_init_boundline_gu_empty	
!����� ���������
call pg_bind_boundline(2)
call pg_init_boundline_gu_val_const(1,1,d0)
end

subroutine init_gu4_3
!������������� ��������� ������� ��� ������ ��������� �������� (������ ���������)
use pgmod2
call pg_allocate_bound_gu
!����
call pg_bind_boundline(1)
call pg_init_boundline_gu_empty	
!����� ���������
call pg_bind_boundline(2)
call pg_init_boundline_gu_val_const(1,1,d0)
call pg_init_boundline_gu_val_const(2,3,d0)
end

subroutine init_gu3_per
!������������� ��������� ������� ��� ������ ��������� ������ � ������
use pgmod2
call pg_allocate_bound_gu
!������ ������
call pg_bind_boundline(1)
call pg_init_boundline_gu_val_const(1,1,d0)
!������ ������ (�����)
call pg_bind_boundline(2)
call pg_init_boundline_gu_val_const(1,2,d0)
!������� ������
call pg_bind_boundline(3)
call pg_init_boundline_gu_val_const(1,1,d0)
!����� ������ (����)
call pg_bind_boundline(4)
call pg_init_boundline_gu_val_const(1,1,d1)
end

subroutine init_gu2_per
!������������� ��������� ������� ��� ������ ��������� ������ � ������
use pgmod2
call pg_allocate_bound_gu
!������ ������ ����� ���������
call pg_bind_boundline(1)
call pg_init_boundline_gu_val_const(1,2,d0)
!������ ������ (�����)
call pg_bind_boundline(2)
call pg_init_boundline_gu_val_const(1,2,d0)
!������� ������
call pg_bind_boundline(3)
call pg_init_boundline_gu_val_const(1,2,d0)
!����� ������ (����)
call pg_bind_boundline(4)
call pg_init_boundline_gu_val_const(1,1,d1)
!������ ����� ����� ���������
call pg_bind_boundline(5)
call pg_init_boundline_gu_val_const(1,2,d0)
!�������
call pg_bind_boundline(6)
call pg_init_boundline_gu_val_const(1,1,d0)
end

subroutine init_gu4_4
!������������� ��������� ������� ��� ������ ��������� ��������� ��������� � ���������� ��������������
use pgmod2
real(8) yc(nmax)
call pg_allocate_bound_gu
!����
call pg_bind_boundline(1)
call pg_get_array_real(2,4,yc,nmax)
call pg_init_boundline_gu_val(1,1,yc)
!����� ��������� 
call pg_bind_boundline(2)
call pg_init_boundline_gu_val_const(1,1,d0)
end

subroutine init_gu_per(periodic)
!�������� ������������� ��������� �������
use pgmod2
integer(4) periodic
if (periodic==3) then
  call pg_init_boundline_gu_periodic(1,1,2,1,1,4,.true.,d0)
else
  call pg_init_boundline_gu_periodic(1,1,2,1,1,4,.false.,h)
endif
end
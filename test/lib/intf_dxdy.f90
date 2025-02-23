module ifdxdy_mod
  real(8), parameter :: ifdxdy_emptyval=Z'0FFFFFFFFFFFFFFF'
end module
  
subroutine intf_dxdy(zm,tr,n,ntr,npe,ff,ffknow,dfdl,ds)
!����������� ����������� ������� �� ��������� ����������� �� ������������� ����������� df/dl
use gen_mod
use slau_block
use ifdxdy_mod
integer(4) n !���������� ����� �����
integer(4) ntr !���������� ��������� �����
integer(4) npe !3,4 - ���������� ����� � ��������
complex(8) zm(n) !���������� ����� �����
integer(4) tr(npe,ntr) !������� ����� ��������� �����
real(8) ff(n) !������� ������� � ����� �����
logical ffknow(n) !������� ������� � ����� (true - ������)
real(8) dfdl !������� ��� ���������� ����������� � ������������ ����� �� ����������� dfdl(x,y,xdir,ydir)
real(8) ds !��� �������������� ����� ��������� ������� �����
           !=0, ���� ���������� ���������� �����, ���������� �� ���������� ������������� ����������� ��� ��������� �����
external dfdl
call intf_dxdy_(zm,tr,n,ntr,npe,ff,ffknow,dfdl,ds,0)
!call intf_dxdy_old(zm,tr,n,ntr,npe,ff,ffknow,dfdl,ds,0)
end

subroutine intf_dxdy_sd(zm,tr,n,ntr,npe,ff,ffknow,dfdl,ds,sd)
!����������� ����������� ������� �� ��������� ����������� �� ������������� ����������� df/dl
use gen_mod
use slau_block
use ifdxdy_mod
integer(4) n !���������� ����� �����
integer(4) ntr !���������� ��������� �����
integer(4) npe !3,4 - ���������� ����� � ��������
complex(8) zm(n) !���������� ����� �����
integer(4) tr(npe,ntr) !������� ����� ��������� �����
real(8) ff(n) !������� ������� � ����� �����
logical ffknow(n) !������� ������� � ����� (true - ������)
real(8) dfdl !������� ��� ���������� ����������� � ������������ ����� �� ����������� dfdl(x,y,xdir,ydir)
real(8) ds !��� �������������� ����� ��������� ������� �����
           !=0, ���� ���������� ���������� �����, ���������� �� ���������� ������������� ����������� ��� ��������� �����
integer(2) sd(ntr) !������� �����������. ������ ���� >0. ���� ������=0, �� ��� ���������� �������������� ������ � ���������� �������
external dfdl
call intf_dxdy_sd_(zm,tr,n,ntr,npe,ff,ffknow,dfdl,ds,0,sd)
end

subroutine intf_dxdy2(zm,tr,n,ntr,npe,ff,ffknow,dfdl,ds)
!����������� ����������� ������� �� ��������� ����������� �� ������������� ����������� df/dl
!��� ������� dfdl ���������� ����� ����� �����
use gen_mod
use slau_block
use ifdxdy_mod
integer(4) n !���������� ����� �����
integer(4) ntr !���������� ��������� �����
integer(4) npe !3,4 - ���������� ����� � ��������
complex(8) zm(n) !���������� ����� �����
integer(4) tr(npe,ntr) !������� ����� ��������� �����
real(8) ff(n) !������� ������� � ����� �����
logical ffknow(n) !������� ������� � ����� (true - ������)
real(8) dfdl !������� ��� ���������� ����������� � ������������ ����� �� ����������� dfdl(x,y,xdir,ydir,itr,ib,k1,k2,is_bound)
             !itr - ����� �������� (������ ������ � ������� tr)
             !ib - ����� ����� � �������� (������ ������ � ������� tr)
             !k1,k2 - ������� ������ � ������� zm
             !is_bound - ����� ���������
real(8) ds !��� �������������� ����� ��������� ������� �����
           !=0, ���� ���������� ���������� �����, ���������� �� ���������� ������������� ����������� ��� ��������� �����
external dfdl
call intf_dxdy_(zm,tr,n,ntr,npe,ff,ffknow,dfdl,ds,1)
!call intf_dxdy_old(zm,tr,n,ntr,npe,ff,ffknow,dfdl,ds,1)
end

subroutine intf_dxdy2_sd(zm,tr,n,ntr,npe,ff,ffknow,dfdl,ds,sd)
!����������� ����������� ������� �� ��������� ����������� �� ������������� ����������� df/dl
!��� ������� dfdl ���������� ����� ����� �����
use gen_mod
use slau_block
use ifdxdy_mod
integer(4) n !���������� ����� �����
integer(4) ntr !���������� ��������� �����
integer(4) npe !3,4 - ���������� ����� � ��������
complex(8) zm(n) !���������� ����� �����
integer(4) tr(npe,ntr) !������� ����� ��������� �����
real(8) ff(n) !������� ������� � ����� �����
logical ffknow(n) !������� ������� � ����� (true - ������)
real(8) dfdl !������� ��� ���������� ����������� � ������������ ����� �� ����������� dfdl(x,y,xdir,ydir,itr,ib,k1,k2,is_bound)
             !itr - ����� �������� (������ ������ � ������� tr)
             !ib - ����� ����� � �������� (������ ������ � ������� tr)
             !k1,k2 - ������� ������ � ������� zm
             !is_bound - ����� ���������
real(8) ds !��� �������������� ����� ��������� ������� �����
           !=0, ���� ���������� ���������� �����, ���������� �� ���������� ������������� ����������� ��� ��������� �����
integer(2) sd(ntr) !������� �����������. ������ ���� >0. ���� ������=0, �� ��� ���������� �������������� ������ � ���������� �������
external dfdl
call intf_dxdy_sd_(zm,tr,n,ntr,npe,ff,ffknow,dfdl,ds,1,sd)
end

subroutine intf_dxdy_old(zm,tr,n,ntr,npe,ff,ffknow,dfdl,ds,mode)
!����������� ����������� ������� �� ��������� ����������� �� ������������� ����������� df/dl
!��� ������� dfdl ���������� ����� ����� �����
use gen_mod
use slau_block
use ifdxdy_mod
integer(4) n !���������� ����� �����
integer(4) ntr !���������� ��������� �����
integer(4) npe !3,4 - ���������� ����� � ��������
complex(8) zm(n) !���������� ����� �����
integer(4) tr(npe,ntr) !������� ����� ��������� �����
real(8) ff(n) !������� ������� � ����� �����
logical ffknow(n) !������� ������� � ����� (true - ������)
real(8) dfdl !������� ��� ���������� ����������� � ������������ ����� �� ����������� dfdl(x,y,xdir,ydir,itr,ib,k1,k2,is_bound)
               !itr - ����� ��������
               !ib - ����� ����� � ��������
               !k1,k2 - ������� ������ � ������� zm
               !is_bound - ����� ���������
             !���� �������� ����������� ��� ��������� ����� �� ����������, ������� ������ ������� �������� ifdxdy_emptyval
real(8) ds !��� �������������� ����� ��������� ������� �����
           !=0, ���� ���������� ���������� �����, ���������� �� ���������� ������������� ����������� ��� ��������� �����
integer(4) mode !0 - intf_dxdy
                !1 - intf_dxdy2
external dfdl
logical, allocatable :: trbound(:,:)
integer(4), allocatable :: ind(:)
real(8), allocatable :: m(:,:),b(:),x(:)
integer(4) i,j,j1,nn,nb,nt,k,k1,nj,ni
complex(8) dz,z
real(8) eps,df,dl,norm,find_df
type(sparse_matrix) sparse
integer(4) matrix_type !0 - ������� �������
                       !1 - sparse
logical have_empty
logical, allocatable :: ff_notempty(:) !������� ������� � ����� ��������� � ����� ���� ���������
!-------------
have_empty=.false.
eps=1.0d-8
matrix_type=0
allocate(trbound(npe,ntr),ind(n),ff_notempty(n))
call get_trbound(tr,n,ntr,npe,trbound)
!������ �������� �����������
ind=0
ff_notempty=.false.
nj=0 !����� �����������
do i=1,n
  if (ffknow(i)) cycle
  nj=nj+1
  ind(i)=nj
enddo
!����� ��������� �����
nb=0
do i=1,ntr
  do j=1,npe
    if (trbound(j,i)) nb=nb+1
  enddo
enddo
!������� ����
nn=(ntr*npe-nb)/2+nb !��������� ����� ���������
if (matrix_type==0) then
  allocate(m(nn,nj))
  m=d0
else
  allocate(sparse%m(nn*2))
  allocate(sparse%c(nn*2))
  allocate(sparse%r(nn+1))
  ni=0
endif
allocate(b(nn))
b=d0
nt=0
do i=1,ntr
  do j=1,npe
    j1=j+1
    if (j1>npe) j1=1
    k=tr(j,i)
    k1=tr(j1,i)
    if (ffknow(k).and.ffknow(k1)) cycle
    if ((.not.trbound(j,i)).and.k1<k) cycle
    dz=zm(k1)-zm(k)
    z=(zm(k1)+zm(k))*d5
    dl=cdabs(dz)
    if (dl<eps) cycle
    if (ds==d0) then
      if (mode==0) then
        df=dfdl(dreal(z),dimag(z),dreal(dz),dimag(dz))
      else
        df=dfdl(dreal(z),dimag(z),dreal(dz),dimag(dz),i,j,k,k1,trbound(j,i))
      endif
    else
      df=find_df(zm(k),zm(k1),ds,dfdl,i,j,k,k1,trbound(j,i),mode)
      dl=d1
    endif
    if (df==ifdxdy_emptyval) then
      have_empty=.true.
      cycle
    endif
    nt=nt+1
    b(nt)=df*dl
    if (ffknow(k)) then
      b(nt)=b(nt)+ff(k)
      if (matrix_type==0) then
        m(nt,ind(k1))=d1
      else
        ni=ni+1
        sparse%r(nt)=ni
        sparse%m(ni)=d1
        sparse%c(ni)=ind(k1)
      endif
      ff_notempty(k1)=.true.
    elseif (ffknow(k1)) then
      b(nt)=b(nt)-ff(k1)
      if (matrix_type==0) then
        m(nt,ind(k))=-d1
      else
        ni=ni+1
        sparse%r(nt)=ni
        sparse%m(ni)=-d1
        sparse%c(ni)=ind(k)
      endif
      ff_notempty(k)=.true.
    else
      if (matrix_type==0) then
        m(nt,ind(k))=-d1
        m(nt,ind(k1))=d1
      else
        ni=ni+1
        sparse%r(nt)=ni
        sparse%m(ni)=d1
        sparse%c(ni)=ind(k1)
        ni=ni+1
        sparse%m(ni)=-d1
        sparse%c(ni)=ind(k)
      endif
      ff_notempty(k)=.true.
      ff_notempty(k1)=.true.
    endif
  enddo
enddo
if (matrix_type==1) then
  sparse%r(nt+1)=ni+1
  sparse%nc=ni
endif
deallocate(trbound)
if (have_empty) then
  !���� ���� ����������� ��� ������� ��� �� ������ ��������� - �������� ��� ��� ��������� ���������
  do i=1,n
    if (ind(i)>0.and.(.not.ff_notempty(i))) then
      nt=nt+1
      if (matrix_type==0) then
        m(nt,ind(i))=d1
      else
        !!!�� ��������!!!
        ni=ni+1
        sparse%r(nt)=ni
        sparse%m(ni)=d1
        sparse%c(ni)=ind(k)
      endif
    endif
  enddo
endif
allocate(x(nj))
if (matrix_type==0) then
  call calc_slau_ls(m, b, nt, nj, nn, x, 1,.true.,norm)
else
  call calc_slau_dss(sparse%m, sparse%c, sparse%r, b, nt, nj, sparse%nc, x)
  !call calc_slau_pardiso(sparse%m, sparse%c, sparse%r, b, nt, nj, sparse%nc, x, 6)
endif
do i=1,n
  if (ind(i)>0) then
    if (ff_notempty(i)) then
      ff(i)=x(ind(i))
    else
      ff(i)=ifdxdy_emptyval
    endif
  endif
enddo
deallocate(b,x,ind,ff_notempty)
if (matrix_type==0) then
  deallocate(m)
else
  call bm_deallocate_sparse_m(sparse)
endif
end

subroutine intf_dxdy_(zm,tr,n,ntr,npe,ff,ffknow,dfdl,ds,mode)
!����������� ����������� ������� �� ��������� ����������� �� ������������� ����������� df/dl
!��� ������� dfdl ���������� ����� ����� �����
use gen_mod
use slau_block
use ifdxdy_mod
integer(4) n !���������� ����� �����
integer(4) ntr !���������� ��������� �����
integer(4) npe !3,4 - ���������� ����� � ��������
complex(8) zm(n) !���������� ����� �����
integer(4) tr(npe,ntr) !������� ����� ��������� �����
real(8) ff(n) !������� ������� � ����� �����
logical ffknow(n) !������� ������� � ����� (true - ������)
real(8) dfdl !������� ��� ���������� ����������� � ������������ ����� �� ����������� dfdl(x,y,xdir,ydir,itr,ib,k1,k2,is_bound)
               !itr - ����� ��������
               !ib - ����� ����� � ��������
               !k1,k2 - ������� ������ � ������� zm
               !is_bound - ����� ���������
             !���� �������� ����������� ��� ��������� ����� �� ����������, ������� ������ ������� �������� ifdxdy_emptyval
real(8) ds !��� �������������� ����� ��������� ������� �����
           !=0, ���� ���������� ���������� �����, ���������� �� ���������� ������������� ����������� ��� ��������� �����
integer(4) mode !0 - intf_dxdy
                !1 - intf_dxdy2
integer(2) sd(1) 
external dfdl
integer(2), allocatable :: trbound(:,:) !������ ������� �������� ������
logical, allocatable :: usebound(:,:) !true - �� ������� ������������ ���������
integer(4), allocatable :: ind(:)
integer(4) i,j,nb,j1
logical add_nb

allocate(trbound(npe,ntr))
call get_trbound_sd(tr,sd,n,ntr,npe,trbound,.false.)
allocate(usebound(npe,ntr),ind(n))
usebound=.false.
ind=0
nb=0
do i=1,ntr
  do j=1,npe
    add_nb=trbound(j,i)<0
    if (.not.add_nb) then
      j1=j+1
      if (j1>npe) j1=1
      add_nb=tr(j1,i)>tr(j,i)
    endif
    if (add_nb) then
      nb=nb+1
      usebound(j,i)=.true.
    endif
  enddo
enddo
call intf_dxdy_one_domain(zm,tr,n,ntr,npe,ff,ffknow,dfdl,ds,mode,usebound,nb,ind,trbound)
deallocate(trbound,usebound,ind)
end

subroutine intf_dxdy_sd_(zm,tr,n,ntr,npe,ff,ffknow,dfdl,ds,mode,sd)
!����������� ����������� ������� �� ��������� ����������� �� ������������� ����������� df/dl
!��� ������� dfdl ���������� ����� ����� �����
use gen_mod
use slau_block
use ifdxdy_mod
integer(4) n !���������� ����� �����
integer(4) ntr !���������� ��������� �����
integer(4) npe !3,4 - ���������� ����� � ��������
complex(8) zm(n) !���������� ����� �����
integer(4) tr(npe,ntr) !������� ����� ��������� �����
real(8) ff(n) !������� ������� � ����� �����
logical ffknow(n) !������� ������� � ����� (true - ������)
real(8) dfdl !������� ��� ���������� ����������� � ������������ ����� �� ����������� dfdl(x,y,xdir,ydir,itr,ib,k1,k2,is_bound)
               !itr - ����� ��������
               !ib - ����� ����� � ��������
               !k1,k2 - ������� ������ � ������� zm
               !is_bound - ����� ���������
             !���� �������� ����������� ��� ��������� ����� �� ����������, ������� ������ ������� �������� ifdxdy_emptyval
real(8) ds !��� �������������� ����� ��������� ������� �����
           !=0, ���� ���������� ���������� �����, ���������� �� ���������� ������������� ����������� ��� ��������� �����
integer(4) mode !0 - intf_dxdy
                !1 - intf_dxdy2
integer(2) sd(ntr) !������� �����������                
external dfdl
integer(2), allocatable :: trbound(:,:) !������ ������� �������� ������
logical, allocatable :: usebound(:,:) !true - �� ������� ������������ ���������
integer(2) max_sd,isd
integer(4), allocatable :: sd_cell_count(:),ind(:)
integer(4) i0,i,j,nb,tb,j1,k,k1
logical ub,need_test_k_k1

allocate(trbound(npe,ntr))
call get_trbound_sd(tr,sd,n,ntr,npe,trbound,.true.)
max_sd=maxval(sd)
allocate(sd_cell_count(0:max_sd))
sd_cell_count=0
do i=1,ntr
  isd=max(sd(i),0)
  sd_cell_count(isd)=sd_cell_count(isd)+1
enddo
allocate(usebound(npe,ntr),ind(n))
do i0=0,max_sd
  if (i0>0.and.sd_cell_count(i0)==0) cycle
  nb=0
  ind=-1
  do i=1,ntr
    do j=1,npe
      ub=.false.
      tb=trbound(j,i)
      need_test_k_k1=.false.
      if (sd(i)/=tb) then
        !��������� ����� ����� ������������ ��� ����������
        if (i0==0) then
          need_test_k_k1=tb>=0
          ub=.not.need_test_k_k1
        endif
      else
        need_test_k_k1=sd(i)==i0
      endif
      if (need_test_k_k1) then
        !��������� ������� �����, ����� ��� ������ ����� �� ��������� ��������� ������
        j1=j+1
        if (j1>npe) j1=1
        k1=tr(j1,i)
        k=tr(j,i)
        ub=k1>k
      endif
      usebound(j,i)=ub
      if (ub) then 
        nb=nb+1
        if (.not.need_test_k_k1) then
          j1=j+1
          if (j1>npe) j1=1
          k1=tr(j1,i)
          k=tr(j,i)
        endif
        ind(k)=0
        ind(k1)=0
      endif
    enddo
  enddo
  if (nb==0) cycle
  call intf_dxdy_one_domain(zm,tr,n,ntr,npe,ff,ffknow,dfdl,ds,mode,usebound,nb,ind,trbound)
enddo
deallocate(trbound,sd_cell_count,usebound,ind)
end

subroutine intf_dxdy_one_domain(zm,tr,n,ntr,npe,ff,ffknow,dfdl,ds,mode,usebound,nn,ind,trbound)
!����������� ����������� ������� �� ��������� ����������� �� ������������� ����������� df/dl
!��� ������� dfdl ���������� ����� ����� �����
use gen_mod
use slau_block
use ifdxdy_mod
integer(4) n !���������� ����� �����
integer(4) ntr !���������� ��������� �����
integer(4) npe !3,4 - ���������� ����� � ��������
complex(8) zm(n) !���������� ����� �����
integer(4) tr(npe,ntr) !������� ����� ��������� �����
real(8) ff(n) !������� ������� � ����� �����
logical ffknow(n) !������� ������� � ����� (true - ������)
real(8) dfdl !������� ��� ���������� ����������� � ������������ ����� �� ����������� dfdl(x,y,xdir,ydir,itr,ib,k1,k2,is_bound)
               !itr - ����� ��������
               !ib - ����� ����� � ��������
               !k1,k2 - ������� ������ � ������� zm
               !is_bound - ����� ���������
             !���� �������� ����������� ��� ��������� ����� �� ����������, ������� ������ ������� �������� ifdxdy_emptyval
real(8) ds !��� �������������� ����� ��������� ������� �����
           !=0, ���� ���������� ���������� �����, ���������� �� ���������� ������������� ����������� ��� ��������� �����
integer(4) mode !0 - intf_dxdy
                !1 - intf_dxdy2
logical usebound(npe,ntr) !true - �� ������� ������������ ���������
integer(4) nn  !���������� �����, �� ������� ����� ������������ ��������� - ��������� ����� ���������
external dfdl
integer(4) ind(n)  !�� �����: -1 - ���� �� ��������� � �������, 0 - ���������
                   !� ���������� ������ �������� �����������
integer(2) trbound(npe,ntr) !������ ������� �������� ������, -1 - ��������� �����
real(8), allocatable :: m(:,:),b(:),x(:)
integer(4) i,j,j1,nt,k,k1,nj,ni
complex(8) dz,z
real(8) eps,df,dl,norm,find_df
type(sparse_matrix) sparse
integer(4) matrix_type !0 - ������� �������
                       !1 - sparse
logical have_empty
logical, allocatable :: ff_notempty(:) !������� ������� � ����� ��������� � ����� ���� ���������
!-------------
have_empty=.false.
eps=1.0d-8
matrix_type=0
allocate(ff_notempty(n))
!������ �������� �����������
ff_notempty=.false.
nj=0 !����� �����������
do i=1,n
  if (ffknow(i).or.ind(i)<0) cycle
  nj=nj+1
  ind(i)=nj
enddo
!������� ����
if (matrix_type==0) then
  allocate(m(nn,nj))
  m=d0
else
  allocate(sparse%m(nn*2))
  allocate(sparse%c(nn*2))
  allocate(sparse%r(nn+1))
  ni=0
endif
allocate(b(nn))
b=d0
nt=0
do i=1,ntr
  do j=1,npe
    j1=j+1
    if (j1>npe) j1=1
    k=tr(j,i)
    k1=tr(j1,i)
    if (ffknow(k).and.ffknow(k1)) cycle
    if (.not.usebound(j,i)) cycle
    dz=zm(k1)-zm(k)
    z=(zm(k1)+zm(k))*d5
    dl=cdabs(dz)
    if (dl<eps) cycle
    if (ds==d0) then
      if (mode==0) then
        df=dfdl(dreal(z),dimag(z),dreal(dz),dimag(dz))
      else
        df=dfdl(dreal(z),dimag(z),dreal(dz),dimag(dz),i,j,k,k1,trbound(j,i)<0)
      endif
    else
      df=find_df(zm(k),zm(k1),ds,dfdl,i,j,k,k1,trbound(j,i)<0,mode)
      dl=d1
    endif
    if (df==ifdxdy_emptyval) then
      have_empty=.true.
      cycle
    endif
    nt=nt+1
    b(nt)=df*dl
    if (ffknow(k)) then
      b(nt)=b(nt)+ff(k)
      if (matrix_type==0) then
        m(nt,ind(k1))=d1
      else
        ni=ni+1
        sparse%r(nt)=ni
        sparse%m(ni)=d1
        sparse%c(ni)=ind(k1)
      endif
      ff_notempty(k1)=.true.
    elseif (ffknow(k1)) then
      b(nt)=b(nt)-ff(k1)
      if (matrix_type==0) then
        m(nt,ind(k))=-d1
      else
        ni=ni+1
        sparse%r(nt)=ni
        sparse%m(ni)=-d1
        sparse%c(ni)=ind(k)
      endif
      ff_notempty(k)=.true.
    else
      if (matrix_type==0) then
        m(nt,ind(k))=-d1
        m(nt,ind(k1))=d1
      else
        ni=ni+1
        sparse%r(nt)=ni
        sparse%m(ni)=d1
        sparse%c(ni)=ind(k1)
        ni=ni+1
        sparse%m(ni)=-d1
        sparse%c(ni)=ind(k)
      endif
      ff_notempty(k)=.true.
      ff_notempty(k1)=.true.
    endif
  enddo
enddo
if (matrix_type==1) then
  sparse%r(nt+1)=ni+1
  sparse%nc=ni
endif
if (have_empty) then
  !���� ���� ����������� ��� ������� ��� �� ������ ��������� - �������� ��� ��� ��������� ���������
  do i=1,n
    if (ind(i)>0.and.(.not.ff_notempty(i))) then
      nt=nt+1
      if (matrix_type==0) then
        m(nt,ind(i))=d1
      else
        !!!�� ��������!!!
        ni=ni+1
        sparse%r(nt)=ni
        sparse%m(ni)=d1
        sparse%c(ni)=ind(k)
      endif
    endif
  enddo
endif
allocate(x(nj))
if (matrix_type==0) then
  call calc_slau_ls(m, b, nt, nj, nn, x, 1,.true.,norm)
else
  call calc_slau_dss(sparse%m, sparse%c, sparse%r, b, nt, nj, sparse%nc, x)
  !call calc_slau_pardiso(sparse%m, sparse%c, sparse%r, b, nt, nj, sparse%nc, x, 6)
endif
do i=1,n
  if (ind(i)>0) then
    if (ff_notempty(i)) then
      ff(i)=x(ind(i))
      ffknow(i)=.true.
    else
      ff(i)=ifdxdy_emptyval
    endif
  endif
enddo
deallocate(b,x,ff_notempty)
if (matrix_type==0) then
  deallocate(m)
else
  call bm_deallocate_sparse_m(sparse)
endif
end

function find_df(z1,z2,ds,dfdl,itr,ib,k1,k2,is_bound,mode)
use gen_mod
use ifdxdy_mod
use func_mod
integer(4) mode !0 - �������� ������ ����������
                !1 - �������� ���������� � ������� �����
complex(8) z1,z2,dz
real(8) dfdl,ds,find_df,d_int
external dfdl
integer(4) n,i
integer(4) itr,ib,k1,k2
logical is_bound
real(8), allocatable :: x(:),y(:),s(:),f(:)
call line_array_ds(dreal(z1),dimag(z1),dreal(z2),dimag(z2),ds,x,y,s,n,4)
allocate(f(n))
dz=z2-z1
find_df=ifdxdy_emptyval
do i=1,n
  if (mode==0) then
    f(i)=dfdl(x(i),y(i),dreal(dz),dimag(dz))
  else
    f(i)=dfdl(x(i),y(i),dreal(dz),dimag(dz),itr,ib,k1,k2,is_bound)
  endif
  if (f(i)==ifdxdy_emptyval) return
enddo
find_df=d_int(n,s,f,s(1),s(n))
deallocate(x,y,s,f)
end

subroutine get_trbound(tr,n,ntr,npe,trbound)
!����������� ����� ���������, ���������� ����������
use gen_mod
integer(4) n !���������� �����
integer(4) ntr !���������� ���������
integer(4) npe !3,4 - ���������� ����� � ��������
integer(4) tr(npe,ntr) !������� ����� ���������
logical trbound(npe,ntr) !������� ���������� �����
logical, allocatable :: m(:,:)
integer(4) i,j,j1,k,k1
allocate(m(n,n))
m=.false.
do i=1,ntr
  do j=1,npe
    j1=j+1
    if (j1>npe) j1=1
    k=min(tr(j,i),tr(j1,i))
    k1=max(tr(j,i),tr(j1,i))
    m(k,k1)=.not.m(k,k1)
  enddo
enddo
do i=2,n
  do j=1,i-1
    m(i,j)=m(j,i)
  enddo
enddo
do i=1,ntr
  do j=1,npe
    j1=j+1
    if (j1>npe) j1=1
    trbound(j,i)=m(tr(j,i),tr(j1,i))
  enddo
enddo
deallocate(m)
end

subroutine get_trbound_sd(tr,sd,n,ntr,npe,trbound,use_sd)
!����������� ����� ���������, ���������� ����������
use gen_mod
integer(4) n !���������� �����
integer(4) ntr !���������� ���������
integer(4) npe !3,4 - ���������� ����� � ��������
integer(4) tr(npe,ntr) !������� ����� ���������
integer(2) sd(ntr) !������� �����������
integer(2) trbound(npe,ntr) !������ ������� �������� ������
logical use_sd    !������ sd �� ��������� � ��������� ������� �����������
integer(2), allocatable :: m(:,:)
integer(4) i,j,j1,k
allocate(m(n,n))
m=-1
do i=1,ntr
  do j=1,npe
    j1=j+1
    if (j1>npe) j1=1
    m(tr(j,i),tr(j1,i))=i
  enddo
enddo
trbound=-1
do i=1,ntr
  do j=1,npe
    j1=j+1
    if (j1>npe) j1=1
    k=m(tr(j1,i),tr(j,i)) !������ �������� ������ (-1 - ��� ������)
    if (k>0) then
      if (use_sd) then
        trbound(j,i)=sd(k)
      else
        trbound(j,i)=0
      endif
    endif
  enddo
enddo
deallocate(m)
end

!Sparse QR Routines
!https://software.intel.com/en-us/mkl-developer-reference-fortran-sparse-qr-routines

!�������� ����������� *.lib
!https://software.intel.com/content/www/us/en/develop/tools/oneapi/components/onemkl/link-line-advisor.html
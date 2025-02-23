!��� ������ ������� �������

subroutine pg_init_boundline_gu_val(igu,bndf,bndval)
!dec$ attributes dllexport:: pg_init_boundline_gu_val
!��������������� ��������� ������� ��� ������� ������� �� ����������
use pgmod
integer(4) igu !����� ���������� �������
integer(4) bndf !��� ���������� �������
real(8) bndval(gsbndl%npanel)
type(TBoundline_GU), pointer :: gu
call allocate_gu(igu,1)
if (igu>gsarea%nu) then
  gu=>gsbndl%gu_add(igu-gsarea%nu)
else
  gu=>gsbndl%gu(igu)
endif
gu%bndu=1
gu%bndf=bndf
gu%bndval(:,1)=bndval
end

subroutine pg_init_boundline_gu_val_const(igu,bndf,bndval)
!dec$ attributes dllexport:: pg_init_boundline_gu_val_const
!��������������� ��������� ������� ��� ������� ������� �� ����������
use pgmod
integer(4) igu !����� ���������� �������
integer(4) bndf !��� ���������� �������
real(8) bndval
type(TBoundline_GU), pointer :: gu
if (igu>gsarea%nu) then
  gu=>gsbndl%gu_add(igu-gsarea%nu)
else
  gu=>gsbndl%gu(igu)
endif
gu%bndu=4
gu%bndf=bndf
gu%constval=bndval
end

subroutine pg_init_boundline_gu_empty
!dec$ attributes dllexport:: pg_init_boundline_gu_empty
!��������������� ��������� ������� ��� ������� �������, ����� ��� �� ������
use pgmod
integer(4) i
do i=1,gsarea%nu
  call pg_init_boundline_gu_empty_one(i)
enddo
gsbndl%gu_inited=.true.
end

subroutine pg_init_boundline_gu_empty_one(igu)
!dec$ attributes dllexport:: pg_init_boundline_gu_empty_one
!��������������� ���� ��������� ������� ��� ������� �������, ����� ��� �� ������
use pgmod
integer(4) igu !����� ���������� �������
type(TBoundline_GU), pointer :: gu
gu=>gsbndl%gu(igu)
gu%bndu=0
end

subroutine pg_init_boundline_gu_constval(igu,bndf,constvali)
!dec$ attributes dllexport:: pg_init_boundline_gu_constval
!��������������� ��������� ������� ��� ������� �������, ����� ������� ���������, �� �� ������
use pgmod
integer(4) igu !����� ���������� �������
integer(4) bndf !��� ���������� �������
integer(4) constvali
type(TBoundline_GU), pointer :: gu
if (igu>gsarea%nu) then
  gu=>gsbndl%gu_add(igu-gsarea%nu)
else
  gu=>gsbndl%gu(igu)
endif
gu%bndu=2
gu%bndf=bndf
gu%constvalind=constvali
end

subroutine pg_init_boundline_gu_gen(igu,bndf,is_direct,ia,ibnd,ibndl,n,c0_,second_bnd,bndf2,c)
!dec$ attributes dllexport:: pg_init_boundline_gu_gen
!��������������� ��������� ������� ��� ������� �������, � ���� ����� �����������
use pgmod
integer(4) igu !����� ���������� �������
integer(4) bndf !�������
logical is_direct  !������ ��� �������� ������������ ������� �� �������� �������
integer(4) ia      !������ ������ �������
integer(4) ibnd    !������ ������ �������
integer(4) ibndl   !������ ������� ������� �������
integer(4) n       !���������� ����������� � ������ �����
real(8) c0_          !��������� ���� � ������ �����
logical second_bnd(n) !������� (false) ��� ������ (true) �������
integer(4) bndf2(n)   !��� �������
real(8) c(n) !������������ ��� �����������
type(TBoundline_GU), pointer :: gu
if (igu>gsarea%nu) then
  gu=>gsbndl%gu_add(igu-gsarea%nu)
else
  gu=>gsbndl%gu(igu)
endif
call init_boundline_gu_gen(gu,bndf,is_direct,ia,ibnd,ibndl,n,c0_,second_bnd,bndf2,c)
end

subroutine pg_init_boundline_gu_gen_var(igu,i,var)
!dec$ attributes dllexport:: pg_init_boundline_gu_gen_var
!��������������� ��������� ������� ��� ������� �������, � ���� ����� �����������
!������������� �������� �������������
use pgmod
integer(4) igu !����� ���������� �������
integer(4) i    !����� �������, =0 ��� c0
real(8) var(*)   !������ �������������
type(TBoundline_GU), pointer :: gu
if (igu>gsarea%nu) then
  gu=>gsbndl%gu_add(igu-gsarea%nu)
else
  gu=>gsbndl%gu(igu)
endif
call init_boundline_gu_gen_var(gu,i,var)
end

subroutine pg_init_boundline_gu_genGlobal(igu,bndf,n)
!dec$ attributes dllexport:: pg_init_boundline_gu_genGlobal
!��������������� ��������� ������� ��� ������� �������, ������� ������ ���� �� ���������� ����������
use pgmod
integer(4) igu !����� ���������� �������
integer(4) bndf !�������
integer(4) n       !���������� ����������� � ������ �����
type(TBoundline_GU), pointer :: gu
if (igu>gsarea%nu) then
  gu=>gsbndl%gu_add(igu-gsarea%nu)
else
  gu=>gsbndl%gu(igu)
endif
call init_boundline_gu_genGlobal(gu,bndf,n)
end

subroutine pg_init_boundline_gu_genGlobal_var(igu,i,ind,var)
!dec$ attributes dllexport:: pg_init_boundline_gu_genGlobal_var
!��������������� ��������� ������� ��� ������� �������, ������� ������ ���� �� ���������� ����������
!������������� �������� �������������
use pgmod
integer(4) igu !����� ���������� �������
integer(4) i    !����� �������, =0 ��� c0
integer(4) ind  !������ ���������� ���������� � gs%constvalinds
real(8) var(*)   !������ �������������
type(TBoundline_GU), pointer :: gu
if (igu>gsarea%nu) then
  gu=>gsbndl%gu_add(igu-gsarea%nu)
else
  gu=>gsbndl%gu(igu)
endif
call init_boundline_gu_genGlobal_var(gu,i,ind,var)
end

subroutine pg_init_boundline_gu_genGlobal_const(igu,i,ind,varconst)
!dec$ attributes dllexport:: pg_init_boundline_gu_genGlobal_const
!��������������� ��������� ������� ��� ������� �������, ������� ������ ���� �� ���������� ����������
!������������� �������� �������������, ���� ������ ���������
use pgmod
integer(4) igu !����� ���������� �������
integer(4) i    !����� �������, =0 ��� c0
integer(4) ind  !������ ���������� ���������� � gs%constvalinds
real(8) varconst     !��������� ��� ������ �������������
type(TBoundline_GU), pointer :: gu
if (igu>gsarea%nu) then
  gu=>gsbndl%gu_add(igu-gsarea%nu)
else
  gu=>gsbndl%gu(igu)
endif
call init_boundline_gu_genGlobal_const(gu,i,ind,varconst)
end

subroutine init_boundline_gu_symmetr(gu,bndf,is_direct,ia,ibnd,ibndl,c0_,c1_)
!��������������� ��������� ������� ��� ������� �������, ������� ���������
use pgmod
integer(4) bndf !�������
logical is_direct  !������ ��� �������� ������������ ������� �� �������� �������
integer(4) ia      !������ ������ �������
integer(4) ibnd    !������ ������ �������
integer(4) ibndl   !������ ������� ������� �������
real(8) c0_          !��������� ���� � ������ �����
real(8) c1_          !����������� ��� �����������
logical second_bnd(1) !������� (false) ��� ������ (true) �������
integer(4) bndf2(1)   !��� �������
real(8) c(1) !������������ ��� �����������
type(TBoundline_GU) :: gu
second_bnd(1)=.true.
bndf2(1)=bndf
c(1)=c1_
call init_boundline_gu_gen(gu,bndf,is_direct,ia,ibnd,ibndl,1,c0_,second_bnd,bndf2,c)
end

subroutine pg_init_boundline_gu_symmetr(igu,bndf,is_direct,ia,ibnd,ibndl,c0_,c1_)
!dec$ attributes dllexport:: pg_init_boundline_gu_add_symmetr
!��������������� ��������� ������� ��� ������� �������, ������� ���������
use pgmod
integer(4) igu !����� ���������� �������
integer(4) bndf !�������
logical is_direct  !������ ��� �������� ������������ ������� �� �������� �������
integer(4) ia      !������ ������ �������
integer(4) ibnd    !������ ������ �������
integer(4) ibndl   !������ ������� ������� �������
real(8) c0_          !��������� ���� � ������ �����
real(8) c1_          !����������� ��� �����������
type(TBoundline_GU), pointer :: gu
if (igu>gsarea%nu) then
  gu=>gsbndl%gu_add(igu-gsarea%nu)
else
  gu=>gsbndl%gu(igu)
endif
call init_boundline_gu_symmetr(gu,bndf,is_direct,ia,ibnd,ibndl,c0_,c1_)
end

subroutine init_boundline_gu_gen(gu,bndf,is_direct,ia,ibnd,ibndl,n,c0_,second_bnd,bndf2,c)
!��������������� ��������� ������� ��� ������� �������, ������� ������ ���� �������� �����������
use pgmod
integer(4) bndf !�������
logical is_direct  !������ ��� �������� ������������ ������� �� �������� �������
integer(4) ia      !������ ������ �������
integer(4) ibnd    !������ ������ �������
integer(4) ibndl   !������ ������� ������� �������
integer(4) n       !���������� ����������� � ������ �����
real(8) c0_          !��������� ���� � ������ �����
logical second_bnd(n) !������� (false) ��� ������ (true) �������
integer(4) bndf2(n)   !��� �������
real(8) c(n) !������������ ��� �����������
type(TBoundline_GU) gu
integer(4) i
if (gu%bndl%npanel/=gs%a(ia)%bnd(ibnd)%line(ibndl)%npanel) call gs_print_stop('Error init_boundline_gu_gen!')
call allocate_gu_gen(gu,n)
gu%bndu=3
gu%bndf=bndf
gu%bndg%ia=ia
gu%bndg%ibnd=ibnd
gu%bndg%ibndl=ibndl
gu%bndg%is_direct=is_direct
gu%bndg%carr(0)%c=c0_
gu%bndg%second_bnd=second_bnd
gu%bndg%bndf=bndf2
do i=1,n
  gu%bndg%carr(i)%c=c(i)
enddo
end

subroutine init_boundline_gu_gen_var(gu,i,var)
!��������������� ��������� ������� ��� ������� �������, � ���� ����� �����������
!������������� �������� �������������
use pgmod
type(TBoundline_GU), target :: gu
integer(4) i    !����� �������, =0 ��� c0
real(8) var(gu%bndl%npanel)   !������ �������������
type(TAreaValue_c), pointer :: av
av=>gu%bndg%carr(i)
if (allocated(av%v)) deallocate(av%v)
allocate(av%v(gu%bndl%npanel))
av%v=var
end

subroutine init_boundline_gu_genGlobal(gu,bndf,n)
!��������������� ��������� ������� ��� ������� �������, ������� ������ ���� �� ���������� ����������
use pgmod
integer(4) bndf !�������
integer(4) n       !���������� ����������� � ������ �����
type(TBoundline_GU) :: gu
call allocate_gu_genGlobal(gu,n)
gu%bndu=5
gu%bndf=bndf
end

subroutine init_boundline_gu_genGlobal_var(gu,i,ind,var)
!��������������� ��������� ������� ��� ������� �������, ������� ������ ���� �� ���������� ����������
!������������� �������� �������������
use pgmod
type(TBoundline_GU), target :: gu
integer(4) i    !����� �������, =0 ��� c0
integer(4) ind  !������ ���������� ���������� � gs%constvalinds
real(8) var(gu%bndl%npanel)   !������ �������������
type(TAreaValue_c), pointer :: av
av=>gu%guglob%carr(i)
if (allocated(av%v)) deallocate(av%v)
allocate(av%v(gu%bndl%npanel))
av%v=var
if (i>0) gu%guglob%gu_indsi(i)=ind
end

subroutine init_boundline_gu_genGlobal_const(gu,i,ind,varconst)
!��������������� ��������� ������� ��� ������� �������, ������� ������ ���� �� ���������� ����������
!������������� �������� �������������, ���� ������ ���������
use pgmod
type(TBoundline_GU), target :: gu
integer(4) i    !����� �������, =0 ��� c0
integer(4) ind  !������ ���������� ���������� � gs%constvalinds
real(8) varconst     !��������� ��� ������ �������������
type(TAreaValue_c), pointer :: av
av=>gu%guglob%carr(i)
av%c=varconst
if (i>0) gu%guglob%gu_indsi(i)=ind
end

subroutine pg_set_boundline_gu_panel(igu,iuse)
!dec$ attributes dllexport:: pg_set_boundline_gu_panel
!���������������� �������� ������� �� ������� �������, �� ������� ����� ������� ��������� �� ���������� �������
use pgmod
integer(4) igu !����� ���������� �������
logical iuse(gsbndl%npanel)
if (.not.allocated(gsbndl%use_gu)) call allocate_gu_use
gsbndl%use_gu(:,igu)=iuse
end

subroutine pg_init_bound_gu_val(igu,bndf,bndval)
!dec$ attributes dllexport:: pg_init_bound_gu_val
!��������������� ��������� ������� ��� ������� �� ����������
use pgmod
integer(4) igu !����� ���������� �������
integer(4) bndf !��� ���������� �������
real(8) bndval(gsbnd%npanel)
integer(4) i
do i=1,gsbnd%nline
  call pg_bind_boundline(i)
  call pg_init_boundline_gu_val(igu,bndf,bndval(gsbndl%i_begin:gsbndl%i_end))
enddo
end

subroutine pg_init_bound_gu_val_const(igu,bndf,bndval)
!dec$ attributes dllexport:: pg_init_bound_gu_val_const
!��������������� ��������� ������� ��� ������� �� ����������
use pgmod
integer(4) igu !����� ���������� �������
integer(4) bndf !��� ���������� �������
real(8) bndval
integer(4) i
do i=1,gsbnd%nline
  call pg_bind_boundline(i)
  call pg_init_boundline_gu_val_const(igu,bndf,bndval)
enddo
end

subroutine pg_init_bound_gu_empty
!dec$ attributes dllexport:: pg_init_bound_gu_empty
!��������������� ��������� ������� ��� �������, ����� ��� �� ������
use pgmod
integer(4) i
do i=1,gsbnd%nline
  call pg_bind_boundline(i)
  call pg_init_boundline_gu_empty
enddo
end

subroutine pg_init_bound_gu_constval(igu,bndf,constvali)
!dec$ attributes dllexport:: pg_init_bound_gu_constval
!��������������� ��������� ������� ��� �������, ����� ������� ���������, �� �� ������
use pgmod
integer(4) igu !����� ���������� �������
integer(4) bndf !��� ���������� �������
integer(4) constvali
integer(4) i
do i=1,gsbnd%nline
  call pg_bind_boundline(i)
  call pg_init_boundline_gu_constval(igu,bndf,constvali)
enddo
end

subroutine init_bound_gu_cylider_particle_base(kk_om)
!���������������� ����������� �� �� �������� ������� � 6 ������������
!��� boundLineType=2
!����� ������� ���������� ��� ���������������
use pgmod
type(TBounLineFuncApprox), pointer :: ga
real(8) kk_om
logical use_kk_r
use_kk_r=.false.
!***����������������� �������
!1 - psi
ga=>gsbndl2%ga(1)
ga%typea=1
ga%bnda=2
call pg_allocate_and_set_ga2(1,3)
!3 - omega
ga=>gsbndl2%ga(3)
ga%typea=1
ga%bnda=4 
call pg_allocate_and_set_ga2(3,3)
kk_om=d1
if (use_kk_r) then
  kk_om=d1/gsbndl2%r
  ga%kk(2:3)=kk_om
endif
!4 - domega/dn
ga=>gsbndl2%ga(4)
ga%typea=1
ga%bnda=3
call pg_allocate_and_set_ga2(4,3)
if (gsbndl2%type_ga==0) then
  if (use_kk_r) ga%kk(2:3)=gsbndl2%r**(-2)
else
  ga%kk(2:3)=d1/gsbndl2%r*kk_om
  ga%ga_ref=3
endif
end

subroutine pg_init_bound_gu_cylider_particle
!dec$ attributes dllexport:: pg_init_bound_gu_cylider_particle
!���������������� ����������� �� �� �������� ������� � 6 ������������
!��� boundLineType=2
use pgmod
type(TBounLineFuncApprox), pointer :: ga
real(8) kk_om
call pg_bind_boundline(1)
call init_bound_gu_cylider_particle_base(kk_om)
!2 - dpsi/dn
ga=>gsbndl2%ga(2)
ga%typea=1
ga%bnda=1
call pg_allocate_and_set_ga2(2,3)
ga%val(1)=d0
end

subroutine pg_init_bound_gu_cylider_particle_slip(k_slip)
!dec$ attributes dllexport:: pg_init_bound_gu_cylider_particle_slip
!���������������� ����������� �� �� �������� ������� � 6 ������������
!� �������� ��������������� ��� �����������
!��� boundLineType=2
use pgmod
real(8) k_slip
type(TBounLineFuncApprox), pointer :: ga
real(8) kk_om
call pg_bind_boundline(1)
call init_bound_gu_cylider_particle_base(kk_om)
!2 - dpsi/dn
ga=>gsbndl2%ga(2)
ga%typea=1
ga%bnda=4
call pg_allocate_and_set_ga2(2,3)
ga%kk(1:3)=-k_slip*kk_om
ga%ga_ref=3
end

subroutine pg_update_bound_gu_cylider_particle_slip(k_slip)
!dec$ attributes dllexport:: pg_update_bound_gu_cylider_particle_slip
!���������������� ����������� �� �� �������� ������� � 6 ������������
!� �������� ��������������� ��� �����������
!��� boundLineType=2
use pgmod
real(8) k_slip
type(TBounLineFuncApprox), pointer :: ga,ga_om
call pg_bind_boundline(1)
!3 - om
ga_om=>gsbndl2%ga(3)
!2 - dpsi/dn
ga=>gsbndl2%ga(2)
ga%kk=-k_slip*ga_om%kk
end

subroutine pg_init_bound_gu_cylider_particle_all
!dec$ attributes dllexport:: pg_init_bound_gu_cylider_particle_all
!���������������� ����������� �� �� ���� �������� �������� ������� �������
!��� boundLineType=2
use pgmod
call pg_init_bound_gu_cylider_particle_all2(0,d0)
end

subroutine pg_init_bound_gu_cylider_particle_all2(mode,param)
!dec$ attributes dllexport:: pg_init_bound_gu_cylider_particle_all2
!���������������� ����������� �� �� ���� �������� �������� ������� �������
!��� boundLineType=2
!� ������� ���� �������
use pgmod
integer(4) mode !0 - ����������, 1 - ���������������
real(8) param !k_slip (��� mode=1)
integer(4) i
do i=1,gsarea%nb
  call pg_bind_bound(i)
  if (gsbnd%boundLineType.ne.2) cycle
  call pg_allocate_bound_gu2
  select case (mode)
  case (0)
    call pg_init_bound_gu_cylider_particle
  case (1)
    call pg_init_bound_gu_cylider_particle_slip(param)
  case default
    call gs_print_stop("Error pg_init_bound_gu_cylider_particle_all2!!!")
  endselect
enddo
end

subroutine pg_init_bound_gu_cylider_particle_all_domains
!dec$ attributes dllexport:: pg_init_bound_gu_cylider_particle_all_domains
!���������������� ����������� �� �� ���� �������� �������� ���� ��������
!��� boundLineType=2
use pgmod
call pg_init_bound_gu_cylider_particle_all_domains2(0,d0)
end

subroutine pg_init_bound_gu_cylider_particle_all_domains2(mode,param)
!dec$ attributes dllexport:: pg_init_bound_gu_cylider_particle_all_domains2
!���������������� ����������� �� �� ���� �������� �������� ���� ��������
!��� boundLineType=2
!� ������� ���� �������
use pgmod
integer(4) mode !0 - ����������, 1 - ���������������
real(8) param !k_slip (��� mode=1)
integer(4) i
do i=1,gs%na
  call pg_bind_domain(i)
  call pg_init_bound_gu_cylider_particle_all2(mode,param)
enddo
end

subroutine pg_init_boundline_gu_periodic(ia1,ibnd1,ibndl1,ia2,ibnd2,ibndl2,is_direct,dpsi1)
!dec$ attributes dllexport:: pg_init_boundline_gu_periodic
!���������������� ������������� ��������� ������� ��� ��������� �������, ���������� ��������
use pgmod
integer(4) ia1    !������� 1
integer(4) ibnd1  !������� 1
integer(4) ibndl1  !������� 1 - ����� psi',om' - �������� �����������, psi,om - ���������
integer(4) ia2    !������� 2
integer(4) ibnd2  !������� 2
integer(4) ibndl2  !������� 2 - ����� psi,om - �������� �����������, psi',om' - ��������� 
logical is_direct  !������ ��� �������� (��������� �������) �������������
real(8) dpsi1   !�������� ��� ������� ����
                  !=0 ��� ������ ������������� (���� ���� ��� ������)
                !=h ��� �������� ������������� (���� � ��������� �������)
type(TBoundline), pointer :: bl1,bl2
type(TArea),pointer :: a1,a2
logical second_bnd(1),dir
integer(4) bndf2(1)
real(8) c(1)
real(8) k_psi
dir=.not.is_direct
a1=>gs%a(ia1)
a2=>gs%a(ia2)
bl1=>a1%bnd(ibnd1)%line(ibndl1)
bl2=>a2%bnd(ibnd2)%line(ibndl2)
if (a1%nu.ne.a2%nu) call gs_print_stop('Error in pg_init_boundline_gu_periodic!')
second_bnd(1)=.true.
if (is_direct) then
  k_psi=d1
else
  k_psi=-d1
endif
!psi=dpsi+psi or psi=dpsi-psi 
bndf2(1)=1
c(1)=k_psi
call init_boundline_gu_gen(bl1%gu(1),1,dir,ia2,ibnd2,ibndl2,1,dpsi1,second_bnd,bndf2,c)
!dpsi/dn=-dpsi/dn or dpsi/dn=dpsi/dn
bndf2(1)=2
c(1)=-k_psi
call init_boundline_gu_gen(bl2%gu(1),2,dir,ia1,ibnd1,ibndl1,1,d0,second_bnd,bndf2,c)
if (a1%nu==2) then
  !om=om or om=-om
  bndf2(1)=3
  c(1)=k_psi
  call init_boundline_gu_gen(bl1%gu(2),3,dir,ia2,ibnd2,ibndl2,1,d0,second_bnd,bndf2,c)
  !dom/dn=-dom/dn or dom/dn=dom/dn
  bndf2(1)=4
  c(1)=-k_psi
  call init_boundline_gu_gen(bl2%gu(2),4,dir,ia1,ibnd1,ibndl1,1,d0,second_bnd,bndf2,c) 
endif
bl1%gu_inited=.true.
bl2%gu_inited=.true.
end

subroutine pg_init_boundline_gu_subdomain_continuity(ia1,ibnd1,ibndl1,ia2,ibnd2,ibndl2)
!dec$ attributes dllexport:: pg_init_boundline_gu_subdomain_continuity
!���������������� �� ������������� �� ������� ����������
use gen_mod
integer(4) ia1    !������� 1
integer(4) ibnd1  !������� 1
integer(4) ibndl1  !������� 1
integer(4) ia2    !������� 2
integer(4) ibnd2  !������� 2
integer(4) ibndl2  !������� 2
call pg_init_boundline_gu_periodic(ia1,ibnd1,ibndl1,ia2,ibnd2,ibndl2,.true.,d0)
end

subroutine pg_init_boundline_gu_subdomain_continuity_with_direction_f(ia1,ibnd1,ibndl1,ia2,ibnd2,ibndl2,f_dir)
!dec$ attributes dllexport:: pg_init_boundline_gu_subdomain_continuity_with_direction_f
!���������������� �� ������������� �� ������� ����������
!� ������ ����������� ��������������� ������� �������
use pgmod
integer(4) ia1    !������� 1
integer(4) ibnd1  !������� 1
integer(4) ibndl1  !������� 1
integer(4) ia2    !������� 2
integer(4) ibnd2  !������� 2
integer(4) ibndl2  !������� 2
external f_dir  !������� ����������� (x,y,dx,dy)
type(TArea), pointer :: a1
type(TBound), pointer :: b1
type(TBoundline), pointer :: bl1
integer(4) i !������ ������� ������
real(8) dx,dy
complex(8) dir,n
a1=>gs%a(ia1)
b1=>a1%bnd(ibnd1)
bl1=>b1%line(ibndl1)
i=bl1%i_begin+(bl1%i_end-bl1%i_begin)/2
call f_dir(b1%xc(i),b1%yc(i),dx,dy)
dir=dcmplx(dx,dy) !����������� ��������������� �������
n=-b1%ett(i)*ii !������� ������� � ������� 1-� �������
call init_boundline_gu_subdomain_continuity_with_direction(ia1,ibnd1,ibndl1,ia2,ibnd2,ibndl2,dir,n)
end

subroutine pg_init_boundline_gu_subdomain_continuity_with_direction(ia1,ibnd1,ibndl1,ia2,ibnd2,ibndl2,dx,dy)
!dec$ attributes dllexport:: pg_init_boundline_gu_subdomain_continuity_with_direction
!���������������� �� ������������� �� ������� ����������
!� ������ ����������� ��������������� ������� �������
use pgmod
integer(4) ia1    !������� 1
integer(4) ibnd1  !������� 1
integer(4) ibndl1  !������� 1
integer(4) ia2    !������� 2
integer(4) ibnd2  !������� 2
integer(4) ibndl2  !������� 2
type(TArea), pointer :: a1
type(TBound), pointer :: b1
type(TBoundline), pointer :: bl1
integer(4) i !������ ������� ������
real(8) dx,dy
complex(8) dir,n
a1=>gs%a(ia1)
b1=>a1%bnd(ibnd1)
bl1=>b1%line(ibndl1)
i=bl1%i_begin+(bl1%i_end-bl1%i_begin)/2
dir=dcmplx(dx,dy) !����������� ��������������� �������
n=-b1%ett(i)*ii !������� ������� � ������� 1-� �������
call init_boundline_gu_subdomain_continuity_with_direction(ia1,ibnd1,ibndl1,ia2,ibnd2,ibndl2,dir,n)
end

subroutine init_boundline_gu_subdomain_continuity_with_direction(ia1,ibnd1,ibndl1,ia2,ibnd2,ibndl2,dir,n)
!���������������� �� ������������� �� ������� ����������
!� ������ ����������� ��������������� ������� �������
use pgmod
integer(4) ia1    !������� 1
integer(4) ibnd1  !������� 1
integer(4) ibndl1  !������� 1
integer(4) ia2    !������� 2
integer(4) ibnd2  !������� 2
integer(4) ibndl2  !������� 2
real(8) sp,scal_p
complex(8) dir,n
sp=scal_p(dir,n)
if (sp>d0) then
  call pg_init_boundline_gu_subdomain_continuity(ia2,ibnd2,ibndl2,ia1,ibnd1,ibndl1)
else
  call pg_init_boundline_gu_subdomain_continuity(ia1,ibnd1,ibndl1,ia2,ibnd2,ibndl2)
endif
end

subroutine pg_init_boundline_gu_periodic_perenos(ia1,ibnd1,ibndl1,ia2,ibnd2,ibndl2,is_direct,kC)
!dec$ attributes dllexport:: pg_init_boundline_gu_periodic_perenos
!���������������� ������������� ��������� ������� ��� ��������� �������� � ������������� kC
!C1=kC*C2
!dC1/dn=-kC*dC2/dn
!eta1=kC*eta2
!deta1/dn=-kC*deta2/dn
use pgmod
integer(4) ia1    !������� 1
integer(4) ibnd1  !������� 1
integer(4) ibndl1  !������� 1
integer(4) ia2    !������� 2
integer(4) ibnd2  !������� 2
integer(4) ibndl2  !������� 2
logical is_direct  !������ ��� �������� (��������� �������) �������������: ������ ������������� - ��������� ������������ �������, �������� ������������� - ������ ������������ �������
real(8) kC        !����������� �������������� �������
type(TBoundline), pointer :: bl1,bl2
type(TArea),pointer :: a1,a2
logical second_bnd(1),dir
integer(4) bndf2(1)
real(8) c(1)
dir=.not.is_direct
a1=>gs%a(ia1)
a2=>gs%a(ia2)
if (a1%nu/=a2.nu) call gs_print_stop('Error pg_init_boundline_gu_periodic_perenos!')
bl1=>a1%bnd(ibnd1)%line(ibndl1)
bl2=>a2%bnd(ibnd2)%line(ibndl2)
second_bnd(1)=.true.
!C1=kC*C2
bndf2(1)=1
c(1)=kC
call init_boundline_gu_gen(bl1%gu(1),1,dir,ia2,ibnd2,ibndl2,1,d0,second_bnd,bndf2,c)
!dC1/dn=-kC*dC2/dn
bndf2(1)=2
c(1)=-d1/kC
call init_boundline_gu_gen(bl2%gu(1),2,dir,ia1,ibnd1,ibndl1,1,d0,second_bnd,bndf2,c)
if (a1%nu==2) then
  !eta1=kC*eta2
  bndf2(1)=3
  c(1)=kC
  call init_boundline_gu_gen(bl1%gu(2),3,dir,ia2,ibnd2,ibndl2,1,d0,second_bnd,bndf2,c)
  !deta1/dn=-kC*deta2/dn
  bndf2(1)=4
  c(1)=-d1/kC
  call init_boundline_gu_gen(bl2%gu(2),4,dir,ia1,ibnd1,ibndl1,1,d0,second_bnd,bndf2,c)
endif
bl1%gu_inited=.true.
bl2%gu_inited=.true.
end

subroutine pg_init_boundline_gu_periodic_perenos_symmetr1(ia1,ibnd1,ibndl1,ia2,ibnd2,ibndl2,ia3,ibnd3,ibndl3,is_direct,kC)
!dec$ attributes dllexport:: pg_init_boundline_gu_periodic_perenos_symmetr1
!���������������� ������� ��������� (��������� �� ������� 3 ��������� � ����������� �� ������� 1)
!��� �������������� ���������� ������� ��� ��������� �������� � ������������� kC
!��� ��������� ����� ia3,ibnd3,ibndl3 ������ ���� �������� ����� ��, ��� ��� ������ pg_init_boundline_gu_periodic_perenos
!C1=kC*C2
!dC1/dn=-kC*dC2/dn
!eta1=kC*eta2
!deta1/dn=-kC*deta2/dn
use pgmod
integer(4) ia1    !������� 1
integer(4) ibnd1  !������� 1
integer(4) ibndl1  !������� 1
integer(4) ia2    !������� 2
integer(4) ibnd2  !������� 2
integer(4) ibndl2  !������� 2
integer(4) ia3    !������� 3
integer(4) ibnd3  !������� 3
integer(4) ibndl3  !������� 3
logical is_direct  !������ ��� �������� (��������� �������) �������������: ������ ������������� - ��������� ������������ �������, �������� ������������� - ������ ������������ �������
real(8) kC        !����������� �������������� �������
type(TBoundline), pointer :: bl3
type(TArea),pointer :: a1,a2,a3
logical dir
dir=.not.is_direct
a1=>gs%a(ia1)
a2=>gs%a(ia2)
a3=>gs%a(ia3)
if (a1%nu/=a2.nu.and.a1%nu/=a3.nu) call gs_print_stop('Error pg_init_boundline_gu_periodic_perenos_symmetr1!')
bl3=>a3%bnd(ibnd3)%line(ibndl3)
!C1=kC*C2
call init_boundline_gu_symmetr(bl3%gu(1),1,dir,ia2,ibnd2,ibndl2,d0,d1)
!dC1/dn=-kC*dC2/dn
call init_boundline_gu_gen(bl3%gu(2),2,.not.dir,ia1,ibnd1,ibndl1,1,d0,[.true.],[2],[-d1/kC])
if (a1%nu==2) then
  !pg_allocate_boundline_gu_add ������ ���� ������� �������
  !eta1=kC*eta2
  call init_boundline_gu_symmetr(bl3%gu_add(1),3,dir,ia2,ibnd2,ibndl2,d0,d1)
  !deta1/dn=-kC*deta2/dn
  call init_boundline_gu_gen(bl3%gu_add(2),4,.not.dir,ia1,ibnd1,ibndl1,1,d0,[.true.],[4],[-d1/kC])
endif
bl3%gu_inited=.true.
end

subroutine pg_init_boundline_gu_periodic_perenos_symmetr2(ia1,ibnd1,ibndl1,ia2,ibnd2,ibndl2,ia3,ibnd3,ibndl3,is_direct,kC)
!dec$ attributes dllexport:: pg_init_boundline_gu_periodic_perenos_symmetr2
!���������������� ������� ��������� (��������� �� ������� 3 ��������� � ����������� �� ������� 1)
!��� �������������� ���������� ������� ��� ��������� �������� � ������������� kC
!��� ��������� ����� ia3,ibnd3,ibndl3 ������ ���� �������� ����� ��, ��� ��� ������ pg_init_boundline_gu_periodic_perenos
!C1=kC*C2
!dC1/dn=-kC*dC2/dn
!eta1=kC*eta2
!deta1/dn=-kC*deta2/dn
use pgmod
integer(4) ia1    !������� 1
integer(4) ibnd1  !������� 1
integer(4) ibndl1  !������� 1
integer(4) ia2    !������� 2
integer(4) ibnd2  !������� 2
integer(4) ibndl2  !������� 2
integer(4) ia3    !������� 3
integer(4) ibnd3  !������� 3
integer(4) ibndl3  !������� 3
logical is_direct  !������ ��� �������� (��������� �������) �������������: ������ ������������� - ��������� ������������ �������, �������� ������������� - ������ ������������ �������
real(8) kC        !����������� �������������� �������
type(TBoundline), pointer :: bl3
type(TArea),pointer :: a1,a2,a3
logical dir
dir=.not.is_direct
a1=>gs%a(ia1)
a2=>gs%a(ia2)
a3=>gs%a(ia3)
if (a1%nu/=a2.nu.and.a1%nu/=a3.nu) call gs_print_stop('Error pg_init_boundline_gu_periodic_perenos_symmetr1!')
bl3=>a3%bnd(ibnd3)%line(ibndl3)
!C1=kC*C2
call init_boundline_gu_gen(bl3%gu(1),1,.not.dir,ia2,ibnd2,ibndl2,1,d0,[.true.],[1],[kC])
!dC1/dn=-kC*dC2/dn
call init_boundline_gu_gen(bl3%gu(2),2,dir,ia1,ibnd1,ibndl1,1,d0,[.true.],[2],[d1])
if (a1%nu==2) then
  !pg_allocate_boundline_gu_add ������ ���� ������� �������
  !eta1=kC*eta2
  call init_boundline_gu_gen(bl3%gu_add(1),3,.not.dir,ia2,ibnd2,ibndl2,1,d0,[.true.],[3],[kC])
  !deta1/dn=-kC*deta2/dn
  call init_boundline_gu_gen(bl3%gu_add(2),4,dir,ia1,ibnd1,ibndl1,1,d0,[.true.],[4],[d1])
endif
bl3%gu_inited=.true.
end

subroutine pg_init_boundline_gu_darcy_darcy(ia1,ibnd1,ibndl1,ia2,ibnd2,ibndl2,k1_div_k2)
!dec$ attributes dllexport:: pg_init_boundline_gu_darcy_darcy
!���������������� ��������� ������� ��� ������� ������� ���� �����-����� � ������� ���������������
use pgmod
integer(4) ia1    !������� 1
integer(4) ibnd1  !������� 1
integer(4) ibndl1  !������� 1
integer(4) ia2    !������� 2
integer(4) ibnd2  !������� 2
integer(4) ibndl2  !������� 2
real(8) k1_div_k2  !��������� �������������� � ������ � ������ ��������
type(TBoundline), pointer :: bl1,bl2
type(TArea),pointer :: a1,a2
logical second_bnd(1)
integer(4) bndf2(1)
real(8) c(1)
a1=>gs%a(ia1)
a2=>gs%a(ia2)
bl1=>a1%bnd(ibnd1)%line(ibndl1)
bl2=>a2%bnd(ibnd2)%line(ibndl2)
second_bnd(1)=.true.
!psi_e=psi_i
bndf2(1)=1
c(1)=d1
call init_boundline_gu_gen(bl1%gu(1),1,.false.,ia2,ibnd2,ibndl2,1,d0,second_bnd,bndf2,c)
!dpsi_i/dn=-k1_div_k2*dpsi_e/dn
bndf2(1)=2
c(1)=-k1_div_k2
call init_boundline_gu_gen(bl2%gu(1),2,.false.,ia1,ibnd1,ibndl1,1,d0,second_bnd,bndf2,c)
bl1%gu_inited=.true.
bl2%gu_inited=.true.
end

subroutine pg_init_boundline_gu_stoks_darcy(ia1,ibnd1,ibndl1,ia2,ibnd2,ibndl2,s_darcy,alf_darcy,bet_darcy)
!dec$ attributes dllexport:: pg_init_boundline_gu_stoks_darcy
!���������������� ��������� ������� ��� ������� ������� ���� �����-�����
!������ ������� �����, ������ ������� �����
use pgmod
integer(4) ia1    !������� 1 (�����)
integer(4) ibnd1  !������� 1
integer(4) ibndl1  !������� 1
integer(4) ia2    !������� 2 (�����)
integer(4) ibnd2  !������� 2
integer(4) ibndl2  !������� 2
real(8) s_darcy  !S=1/sqrt(k)
real(8) alf_darcy,bet_darcy  !������������ ��� ��������� V_stoks � V_darcy � ������� ��������� ������� ��������������
type(TBoundline), pointer :: bl1,bl2
type(TArea),pointer :: a1,a2
logical second_bnd(2)
integer(4) bndf2(2)
real(8) c(2)
a1=>gs%a(ia1)
a2=>gs%a(ia2)
bl1=>a1%bnd(ibnd1)%line(ibndl1)
bl2=>a2%bnd(ibnd2)%line(ibndl2)
second_bnd(1)=.true.
!�� ��������� ��������
bndf2(1)=2
c(1)=s_darcy**2
call init_boundline_gu_gen(bl1%gu(1),4,.false.,ia2,ibnd2,ibndl2,1,d0,second_bnd,bndf2,c)


!om_e=S(alf*v_e_tau-bet*v_i_tau)  ������ �������
bndf2(1)=2
bndf2(2)=2
c(1)=-s_darcy*bet_darcy
c(2)=-s_darcy*alf_darcy
second_bnd(1)=.true.
second_bnd(2)=.false.
call init_boundline_gu_gen(bl1%gu(2),3,.false.,ia2,ibnd2,ibndl2,2,d0,second_bnd,bndf2,c)

!!! ������ ��� ������! dome/dn - ������� ��� ��������� ���������� � ������ ��!
!om_e=alf*S*v_e_tau-bet/S*dome/dn  ������ �������, ���������� ����� dome/dn
!bndf2(1)=2
!bndf2(2)=4
!c(1)=-alf_darcy*s_darcy
!c(2)=-bet_darcy/s_darcy
!second_bnd(1)=.false.
!second_bnd(2)=.false.
!call init_boundline_gu_gen(bl1%gu(2),3,.false.,ia2,ibnd2,ibndl2,2,d0,second_bnd,bndf2,c)

!psi_e'=-psi_i' ���������� ����������� ��������� ���������
!bndf2(1)=2
!c(1)=-d1
!second_bnd(1)=.true.
!call init_boundline_gu_gen(bl1%gu(2),2,.false.,ia2,ibnd2,ibndl2,1,d0,second_bnd,bndf2,c)

!psi_e'=0  v_e_tau=0 ������� �������� ������ �� �������
!call init_boundline_gu_gen(bl1%gu(2),2,.false.,ia2,ibnd2,ibndl2,0,d0,second_bnd,bndf2,c)


!psi_e=psi_i
bndf2(1)=1
c(1)=d1
second_bnd(1)=.true.
call init_boundline_gu_gen(bl2%gu(1),1,.false.,ia1,ibnd1,ibndl1,1,d0,second_bnd,bndf2,c)
bl1%gu_inited=.true.
bl2%gu_inited=.true.
end

subroutine pg_init_boundline_gu_stoks_brinkman(ia1,ibnd1,ibndl1,ia2,ibnd2,ibndl2,bet,delt,k_bri,mu_bri)
!dec$ attributes dllexport:: pg_init_boundline_gu_stoks_brinkman
!���������������� ��������� ������� ��� ������� ������� ���� �����-��������
!������ ������� �����, ������ ������� �����
use pgmod
integer(4) ia1    !������� 1 (�����)
integer(4) ibnd1  !������� 1
integer(4) ibndl1  !������� 1
integer(4) ia2    !������� 2 (��������)
integer(4) ibnd2  !������� 2
integer(4) ibndl2  !������� 2
real(8) k_bri,mu_bri  !k � mu_b
real(8) bet,delt  !������������ � ��������� �������� alpha,beta,gamma,delta
type(TBoundline), pointer :: bl1,bl2
type(TArea),pointer :: a1,a2
real(8) sqk
sqk=dsqrt(k_bri)
a1=>gs%a(ia1)
a2=>gs%a(ia2)
bl1=>a1%bnd(ibnd1)%line(ibndl1)
bl2=>a2%bnd(ibnd2)%line(ibndl2)
!psi_e=psi_i
call init_boundline_gu_gen(bl1%gu(1),1,.false.,ia2,ibnd2,ibndl2,1,d0,[.true.],[1],[d1])
!bet*vtau_e-vtau_i=delt*sqrt(k)*omega_e
call init_boundline_gu_gen(bl1%gu(2),2,.false.,ia2,ibnd2,ibndl2,2,d0,[.true.,.false.],[2,3],[-d1/bet,-delt/bet*sqk])
!ome=mub*omi
call init_boundline_gu_gen(bl2%gu(1),3,.false.,ia1,ibnd1,ibndl1,1,d0,[.true.],[3],[d1/mu_bri])
!�� ��������� ��������
call init_boundline_gu_gen(bl2%gu(2),4,.false.,ia1,ibnd1,ibndl1,2,d0,[.false.,.true.],[2,4],[d1/(mu_bri*k_bri),-d1/mu_bri])
bl1%gu_inited=.true.
bl2%gu_inited=.true.
end

subroutine pg_init_boundline_gu_slip(igu,slip_k)
!dec$ attributes dllexport:: pg_init_boundline_gu_slip
!��������������� ��������� ������� ����������
!dpsi/dn=-slip_k*eta    (eta=-om)
!dpsi/dn=slip_k*om
use pgmod
integer(4) igu !����� ���������� �������
real(8) slip_k !����������� ����������
logical second_bnd(1)
integer(4) bndf2(1)
real(8) c(1)
second_bnd(1)=.false.
bndf2(1)=3
c(1)=-slip_k
call pg_init_boundline_gu_gen(igu,2,.true.,gsarea%i,gsbnd%i,gsbndl%i,1,d0,second_bnd,bndf2,c)
end

subroutine pg_init_boundline_gu_slip2(igu,slip_k)
!dec$ attributes dllexport:: pg_init_boundline_gu_slip2
!��������������� ��������� ������� ����������
!dpsi/dn=-slip_k*eta    (eta=-om)
!dpsi/dn=slip_k*om
use pgmod
integer(4) igu !����� ���������� �������
real(8) slip_k !����������� ����������
logical second_bnd(1)
integer(4) bndf2(1)
real(8) c(1)
second_bnd(1)=.false.
bndf2(1)=4
c(1)=-slip_k
call pg_init_boundline_gu_gen(igu,2,.true.,gsarea%i,gsbnd%i,gsbndl%i,1,d0,second_bnd,bndf2,c)
end

subroutine pg_init_gu_multiarea(set_gu_internal,set_gu_external)
!dec$ attributes dllexport:: pg_init_gu_multiarea
!����������������� ��������� ������� ��� ����� � ������������� ��������� ������� �� ��������� ��������
!��� boundLineType=1
use pgmod
external set_gu_internal !��������� ��� ������� ����������� ���������� �������
                         !��������� - (ia1,ib1,ibl1,ia2,ib2,ibl2)
external set_gu_external !��������� ��� ������� �������� ���������� ������� ��� �������� gsbndl                        
                         !��� ����������
integer(4) i1,j1,k1
type(TArea), pointer :: a1
type(TBound), pointer :: b1
type(TBoundline), pointer :: bl1
real(8), parameter :: eps=1.0d-8
do i1=1,gs%na
  call pg_bind_domain(i1)
  do j1=1,gsarea%nb
    call pg_bind_bound(j1)
    call pg_allocate_bound_gu
  enddo
enddo
do i1=1,gs%na
  a1=>gs%a(i1)
  do j1=1,a1%nb
    b1=>a1%bnd(j1)
    if (b1%boundLineType/=1) cycle
    do k1=1,b1%nline
      bl1=>b1%line(k1)
      if (bl1%gu_inited) cycle
      if (bl1%ci%is_internal) then
        call set_gu_internal(i1,j1,k1,bl1%ci%ia,bl1%ci%ibnd,bl1%ci%ibndl)
        bl1%gu_inited=.true.
        gs%a(bl1%ci%ia)%bnd(bl1%ci%ibnd)%line(bl1%ci%ibndl)%gu_inited=.true.
      endif
    enddo
  enddo
enddo
!���������� ������� ������� � �������� ��
do i1=1,gs%na
  a1=>gs%a(i1)
  do j1=1,a1%nb
    b1=>a1%bnd(j1)
    if (b1%boundLineType/=1) cycle
    do k1=1,b1%nline
      bl1=>b1%line(k1)
      if (bl1%gu_inited) cycle
      call pg_bind_domain(i1)
      call pg_bind_bound(j1)
      call pg_bind_boundline(k1)
      call set_gu_external
    enddo
  enddo
enddo
end 

subroutine pg_init_gu_multiarea2(set_gu_external)
!dec$ attributes dllexport:: pg_init_gu_multiarea2
!����������������� ��������� ������� ��� ����� � ������������� ��������� ������� �� ��������� ��������
!��� boundLineType=1
use pgmod
external set_gu_external !��������� ��� ������� �������� ���������� ������� ��� �������� gsbndl                        
                         !��� ����������
integer(4) i1,j1,k1
type(TArea), pointer :: a1,a2
type(TBound), pointer :: b1,b2
type(TBoundline), pointer :: bl1,bl2
real(8), parameter :: eps=1.0d-8
do i1=1,gs%na
  call pg_bind_domain(i1)
  do j1=1,gsarea%nb
    call pg_bind_bound(j1)
    if (gsbnd%boundLineType==1) call pg_allocate_bound_gu
  enddo
enddo
do i1=1,gs%na
  a1=>gs%a(i1)
  do j1=1,a1%nb
    b1=>a1%bnd(j1)
    if (b1%boundLineType/=1) cycle
    do k1=1,b1%nline
      bl1=>b1%line(k1)
      if (bl1%gu_inited) cycle
      if (bl1%ci%is_internal) then
        a2=>gs%a(bl1%ci%ia)
        if (a1%iorder<a2%iorder) then
          call pg_init_boundline_gu_subdomain_continuity(bl1%ci%ia,bl1%ci%ibnd,bl1%ci%ibndl,i1,j1,k1)
        elseif (a1%iorder==a2%iorder) then
          if (a1%i<a2%i) then
            call pg_init_boundline_gu_subdomain_continuity(bl1%ci%ia,bl1%ci%ibnd,bl1%ci%ibndl,i1,j1,k1)
          else
            call pg_init_boundline_gu_subdomain_continuity(i1,j1,k1,bl1%ci%ia,bl1%ci%ibnd,bl1%ci%ibndl)
          endif
        else
          call pg_init_boundline_gu_subdomain_continuity(i1,j1,k1,bl1%ci%ia,bl1%ci%ibnd,bl1%ci%ibndl)
        endif
        bl1%gu_inited=.true.
        b2=>a2%bnd(bl1%ci%ibnd)
        bl2=>b2%line(bl1%ci%ibndl)
        bl2%gu_inited=.true.
      endif
    enddo
  enddo
enddo
!���������� ������� ������� � �������� ��
do i1=1,gs%na
  a1=>gs%a(i1)
  do j1=1,a1%nb
    b1=>a1%bnd(j1)
    if (b1%boundLineType/=1) cycle
    do k1=1,b1%nline
      bl1=>b1%line(k1)
      if (bl1%gu_inited) cycle
      call pg_bind_domain(i1)
      call pg_bind_bound(j1)
      call pg_bind_boundline(k1)
      call set_gu_external
    enddo
  enddo
enddo
end
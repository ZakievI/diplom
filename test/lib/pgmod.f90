module pgmod
    use gen_mod
    use func_mod
    use slau_block

    real(8), parameter :: ds0numint=1.0d-3
    integer(4), parameter :: max_int=19 !������������ ���������� ����� ���������� ��� ����
    !real(8), parameter :: real8_inf=1.0d0/0.0d0
    real(8), parameter :: real8_inf=Z'0FFFFFFFFFFFFFFF'
    integer(4), parameter :: int4_empty=Z'0FFFFFFF'
    integer(4), parameter :: gs_log_file_i=999
    integer(4), parameter :: int_type(0:max_int) = [0,0,0,0,1,2,1,2,1,2,1,2,0,0,0,0,1,2,1,2] !0 - ������� ��������, 1 - d/dx, 2 - d/dy
    integer(4), parameter :: int_type_dn(0:max_int) = [0,1,0,1,0,0,1,1,0,0,1,1,0,1,0,1,0,0,1,1] !0 - ������� ��������, 1 - �������� � d/dn
    integer(4), parameter :: int_type_area(0:max_int) = [1,0,1,0,1,1,0,0,1,1,0,0,0,0,0,0,0,0,0,0] !1 - ��������, ������� ����� ��������� �� �������, 0 - ������
    integer(4), parameter :: int_type_dual(0:max_int) = [-1,-1,-1,-1,5,-2,7,-2,9,-2,11,-2,-1,-1,-1,-1,-1,-1,-1,-1] !����� ������� ���������, -1 - ��� �������, -2 ������� �������� �������� ������ � ������� ���������
    
    type TAreaValue
      real(8), allocatable :: v(:) !(ntr) �������� � ������������
    endtype
    
    type TAreaValue_c
      real(8) c !���������� ��������
      real(8), allocatable :: v(:) !���������� ��������
    endtype
    
    type intf_struct
      integer(4) n
      real(8), allocatable :: x0(:),y0(:),s0(:)
      real(8), pointer :: x(:),y(:),s(:)
      real(8) xder,yder
    endtype
    

  type TCashValue
    real(8), allocatable :: i(:,:) !������ ������ - ������ (�����������) ��������������
                                   !������ ������ - ����������� �����
  endtype
  
  type TCashIntegral
    logical, allocatable :: inited(:) !(0:max_int) true - ������ intvals(k)%i ���������������
                                      !������ - ��� ���������
    logical, allocatable :: solved(:) !(0:max_int) true - ������ intvals(k)%i �������� � ������ gs_cash_presolve
                                      !������ - ��� ���������
    type(TCashValue), allocatable :: intvals(:) !(0:max_int) ��� ����������
  endtype
  
  type TCash 
    !��� ����������, ����������� ������� ��� �������, �� �������� ������� ��������������
    logical bnd_inited
    logical area_inited
    type(TCashIntegral), allocatable :: bnd(:) !(nb) ��� ���������� � ������������ ������� �� ��������
    type(TCashIntegral) area !��� ���������� � ������������ ������� � �������������
    integer(4) ia            !������ �������, ������� ����������� ���
  endtype
  
  type TBoungline_GU_gen
    logical is_direct  !������ ��� �������� ������������ ������� �� �������� �������
    integer(4) ia      !������ ������ �������
    integer(4) ibnd    !������ ������ �������
    integer(4) ibndl   !������ ������� ������� �������
    integer(4) n       !���������� ����������� � ������ �����
    logical, allocatable :: second_bnd(:) !(n) ������� (false) ��� ������ (true) �������
    integer(4), allocatable :: bndf(:)   !(n) ��� �������
                                          !1 - psi
                                          !2 - dpsi/dn
                                          !3 - omega
                                          !4 - domega/dn
    type(TAreaValue_c), allocatable :: carr(:) !(0:n) ������� ������������� ��� �����������. 0 - ��� ���������� ����� �0
  endtype
  
  type TBoungline_GU_genGlobal
    !��������� ������� � ����������� �� ���������� �������� � ���������� ��������������
    integer(4) n       !���������� ����������� � ������ �����
    type(TAreaValue_c), allocatable :: carr(:) !(0:n) ������� ������������� ��� �����������. 0 - ��� ���������� ����� �0
    integer(4), allocatable :: gu_indsi(:)  !(n) ������� ���������� ����������� � gs%constvalinds
  endtype
  
  type TBound_near   !���������� � ������ ��� ���������� �������� ������� ������ �������
    real(8) x,y      !���������� ��������� ����� �� �������
    real(8) dist        !��������� �� �����
    real(8) panelL      !����� ��������� ������
    logical in_domain   !����� ��������� � ������� 
  endtype
  
  type TBound_info !���������� � ����� �� �������
    integer(4) ibnd !����� �������
    real(8) s       !������� ��������
    real(8) tt      !���� ������� ����������� � �������  
    integer(4) in_bound !0 - �������, 1 - �������, 2 - ����������� ������� (���������� ������ ����� ������)
    type(TBound_near) bn
  endtype

	type TBoundline_GU !��������� �������
	  integer(4) bndu  ! ��� ���������� �������
		  !0 - ������ �� ������ 
      !1 - ������ �������
      !2 - �������� ���������� �� ���� �������, �� ������� �����������
      !3 - ����� ������� - �������� ����������� �� ������ ����������� (�������� �� ������ �������, ������ �������, ������ �������)
      !4 - ������ ������� - const
      !5 - �� � ����������� �� ���������� �������� 
    integer(4) bndf  ! �������, ��� ������� �������� ��
      !1 - psi
      !2 - dpsi/dn
      !3 - omega
      !4 - domega/dn
	  integer(4) constvalind !��� bndu=2 - ������ ����������� � constvals, constvalinds
    type(TBoungline_GU_gen) bndg !��� ������ �� bndu=3
	  real(8), allocatable :: bndval(:,:) !(npanel,nval) �������� � ��������� �������. nval = 1
    real(8) constval !�������� ������� ��� bndu=4
    type(TBoundline), pointer :: bndl !������ �� ������� �������
    type(TBoungline_GU_genGlobal) guglob !��� �� � ����������� �� ���������� �������� bndu=5
	endtype

	type TBounLineFuncApprox !������������� ������� �� ������� �������
	  !���������������� ������� ������������ �������� � ������� TBoundline2%ga
	      !1 - psi
          !2 - dpsi/dn
          !3 - omega
          !4 - domega/dn
	  integer(4) typea ! ��� �������������  
	      !1 - ����� ���� ����� (g=2*pi*s, s - ��������� ������������� ������� ��������: 0 - ������, 1 - �����)
	  integer(4) bnda  !��������� �������������  
	      !1 - C0 (��������)
          !2 - �0
          !3 - �1*sin(g)+�2*cos(g)
          !4 - �0+�1*sin(g)+�2*cos(g)
	  real(8), allocatable :: valkk(:)  !�������� ������������ � ������ ��������� �� kk
	  real(8), allocatable :: val(:)  !�������� ������������
	  integer(4), allocatable :: ind(:) !������� ����������� ��� bnda>1
	  real(8), allocatable :: kk(:)  !���������� ��������� ��� �������������
	  integer(4) n !���������� ����������� �������������
    integer(4) ga_ref !������ �� ������ ������������ ������������� (����� dom/dn ������� �� om)
	  real(8) cc !��������� 
	  real(8) delta !�������� ��� (�������� ������� ��������)    
	endtype

	type TBoundline_Collocate !����� ���������� � ��������� ������� ��� ������������������ �������
	  !����� ��������� ������������ �������� � ������� TBoundline2%gc
	  real(8), pointer :: s_gu(:) !��������� ������������� ������� �������� �����, � ������� ���������� ��������� ���������
	  integer(4) n !����� ����� ���������� � s_gu
	  real(8), allocatable :: g(:) !s_gu*pi2
	  complex(8),allocatable :: z(:) !����������� ����� ����������
	  complex(8),allocatable :: ztc(:) !������� �� ������ ����� �� ����� ����������
  endtype
  
  type TBoundline_geomdetale
    integer(4) mode !0 - �� ������
                    !1 - ������ ���������
                    !2 - ������������� �������
                    !3 - ���� ����������
    real(8) s1,s2   !������� �������� ��������� � �������� �����
    real(8) x,y     !������ ����� (mode=1), ����� ���������� (2)
    real(8) x2,y2   !�������� ����� (1)
    real(8) r,gam1,gam2 !������ � ���� (2)
    integer(4) i_begin, i_end !������� ������ � ��������� ����� � �������� x,y � TBound (� ������ � TBoundline)
    integer(4) ibndl    !������ TBoundline
    real(8) panelLmax   !������������ ����� ������ �� �������
  endtype

  type Connect_info
    logical is_internal   !���������� ������� ����� �����������
    integer(4) ia         !����� �������
    integer(4) ibnd       !����� �������
    integer(4) ibndl      !����� �������
  endtype
    
  type TBoundline !������� ����� �������
	    integer(4) i !�����
        !�������������� ��������������
		integer(4) i_begin,i_end       !������ ������ � ��������� ������ ������� � �������� �������
    integer(4) npanel    !����� �������
    real(8), allocatable :: x(:), y(:)   !(npanel+1) ���������� ������ �������        
        !��������� �������
		type(TBoundline_GU), allocatable :: gu(:) !(nu) ��������� �������
    type(TBoundline_GU), allocatable :: gu_add(:) !(nu) �������������� ��������� �������, ���������� ����� ���������� � ��������� � ���
		logical, allocatable :: use_gu(:,:)   !(npanel,nu) - true � ��� �������, ��� ����� �������� ���������
    logical, allocatable :: cp_line(:)   !(nu) - true ��� ���� ���������, ��� �������� ����� ������� ����� ����������
    integer(4) gu_mode  !��� ���������� ������� - ������������ ����� �������������
    integer(4) igd_begin !������ ���������� geomdetale � bnd%geom_detale0
    integer(4) igd_end !������ ��������� geomdetale � bnd%geom_detale0
    logical gu_inited  !��������� ������� ������������������. ������������ ��� pg_init_gu_multiarea
    type(Connect_info) ci  !���������� � �������� �������
  endtype

	type TBoundline2 !������� ����� �������
	    integer(4) i !�����
        !�������������� ��������������
		complex(8) zc !����� �������
		real(8) r !������ �������
		real(8) g0 !��������� ���� (������������� �� ��� x ������ ������� �������)
		integer(4) dir !����������� ������ ����������: 1 - ������ ������� ������� (��� ������������), -1 - �� ������� ������� (��� ���������)
		type(TBoundline_Collocate), allocatable :: gc(:) !(nu) ����� ���������� ���������� ������� ��� ������������������ �������
		!��������� �������
		type(TBounLineFuncApprox), allocatable :: ga(:) !(umax) ������������� ��� �������
    integer(4) type_ga !0 - 6 �����������, 1 - 4 ����������� (dom/dn ���������� ����� om)
    real(8) m1 !��������� �������� ������������ = 1-m
    real(8) u !��������� �������� �������� ����������
    real(8) tet !���� ������� ��������� �������� ����������
  endtype

	type TBound_Collocate_points !����� ���������� � ��������� ������� ��� ������������������ �������
	  !����� ��������� ������������ �������� � ������� TArea%cpp
	  integer(4) ncp !����� ����� ����������
	  integer(4), allocatable :: i(:) !(ncp) ������ ����� (boundLineType=1 - ������ ������, boundLineType=2 - ������ �����)
	  integer(4), allocatable :: ibndl(:) !(ncp) ������ ������� ������� (��� boundLineType /= 1)
	  integer(4), allocatable :: ibnd(:) !(ncp) ������ ������� 
	  integer(4), allocatable :: iu(:) !(ncp) ��������� (1..nu)
	  integer(4), allocatable :: j(:) !(ncp) ������ ��������� � �������
    logical inited
	endtype

	type TArea_Collocate_points !����� ���������� � ��������� ������� ��� ������������������ �������
	  !����� ��������� ������������ �������� � ������� TArea%cppa
	  integer(4) ncp !����� ����� ����������
	  integer(4), allocatable :: iu(:) !(ncp) ����� ��������� � ������������
	  integer(4), allocatable :: itr(:) !(ncp) ������ ������������
	  integer(4), allocatable :: j(:) !(ncp) ������ ���������
    logical inited
	endtype

    type TBound  !����� �������
	    integer(4) i !�����
	    integer(4) nline !����� �������� �������
		integer(4) boundLineType !���

		!***boundLineType=1 - ������� � ������� ������,
	    type(TBoundline), allocatable :: line(:)  !������� �������
      type(TBoundline_geomdetale), allocatable :: geom_detale0(:) !���������� � ��������� �������� ������� (��������)
      type(TBoundline_geomdetale), allocatable :: geom_detale(:) !���������� � ��������� �������� ������� - ����� �����������, ���������������
      integer(4) ngeom_detale0 !������ ���������� ������������������� �������� � geom_detale0
      integer(4) ngeom_detale !���������� ��������� � geom_detale
		!�������������� ��������������
	    integer(4) npanel    !����� �������
        real(8), allocatable :: x(:), y(:)   !(npanel+1) ���������� ������ �������
        real(8), allocatable :: xc(:), yc(:) !(npanel)   ���������� ������� �������
        real(8), allocatable :: s(:)         !(npanel+1) ������� �������� ������
        real(8), allocatable :: sc(:)        !(npanel)   ������� �������� �������
        real(8), allocatable :: l(:)         !(npanel)   ����� �������
        complex(8), allocatable :: z(:)      !(npanel+1) ����������� ���������� ������ �������
        complex(8), allocatable :: zc(:)     !(npanel) ����������� ���������� ������� �������
        complex(8), allocatable :: ett(:)    !(npanel) cdexp(ii*tt)
	  
	    !�����������
        real(8), allocatable :: psiom(:,:)     !(npanel,umax) - �������� ������� ������� � ����������� ������
        integer(4), allocatable :: psiind(:,:) !(npanel,umax) - ������� �����������
                             !1-� ������ - ����� ����������� �����
                             !2-� ������ - 1- psi, 2- dpsi/dn, 3- om, 4- dom/dn
		real(8), allocatable :: bb_psioml(:,:) !(npanel+1,umax)
		real(8), allocatable :: cc_psioml(:,:,:) !(4,npanel+1,umax)

		!***boundLineType=2 - �������� ������� � �������������� ������� �������
		type(TBoundline2), allocatable :: line2(:)  !������� �������
    
    !��� ����������
    type(TCash) self_cash !����������� ���
    type(TCash), pointer :: cash !������ �� ��� ������ ������� ��� ����������� ���
    logical skip_getfun !���������� ������� ��� ���������� �������
  endtype
    
    type areapart
	    !��������� �������
      integer(4) i !�����
      integer(4) n_begin !������ ������ � ������� zm �������
      integer(4) ntr_begin !������ ������ � ������� trm �������
      integer(4) n_end !��������� ������ � ������� zm �������
      integer(4) ntr_end !��������� ������ � ������� trm �������
      integer(4) npe !����� ����� � �������� (3,4) - ��� ��������� ������
      integer(4) n !����� �����
      integer(4) ntr !���������� ������������� � �����
      complex(8), allocatable :: zm(:) !(n) !���������� ����� ����� � �������
      integer(4), allocatable :: trm(:,:) !(npe,ntr) !������� ������ ����������� �����
      !������ �� ������ �������
      type(areapart), pointer :: apref
      complex(8) a_ref,b_ref,c_ref
    endtype
    
    type areatype
	  !���������
      integer(4) npart !���������� ��������
      type(areapart), allocatable :: part(:) !(npart) ������� �������
      logical geom_inited
      integer(4) npe !������������ ����� ����� � �������� (3,4) - ��� ��������� ������
      integer(4), allocatable :: npe_ar(:) !(ntr) ����� ����� � �������� (3,4)
      integer(4) n !����� �����
      complex(8), allocatable :: zm(:) !(n) !���������� ����� ����� � �������
      integer(4), allocatable :: trm(:,:) !(npe,ntr) !������� ������ ����������� �����
      integer(4) ntr !���������� ������������� � �����
      complex(8), allocatable :: zmc(:) !(ntr) !���������� ������� �������������
	  !�������� � �������������
      type(TAreaValue), allocatable :: areaval(:) !(k) !�������� � ������������� ��� ������ �������� � �.�. (k=1,2,3,4 � ����������� �� ������)
                                           !1 - �������� �������� ������� � ������ ����� ������������� ���������
    									   !2 - ����������� (����� ���� ����������) ��� ������� ������� � ��������� �����������
    									   !3 - �������� ������� v � ������ ��������� ������� ���� \Delta\psi=v(x,y)\omega
										   !4 - �������� ������� vx � ��������� �������� ��� ���� ��� df/dy � ������� ������ ���������������� ��������� (eq=18,19)
										                                                !��� ���� ��� 1/y*df/dy � ������� ������ ���������������� ��������� (eq=20) 
										   !5 - �������� ������� vy � ��������� ��������
                       
                       !��� type_eq=22,24,25 (1 - ������������ ��, 2 - ���� ��� psi, 3 - ���� ��� eta, 4 - dpsi/dx, 5 - dpsi/dy, 6 - deta/dx, 7 - deta/dy)
                       !��� type_eq=26 (1 - ������������ ��, 2 - ���� ��� psi, 3 - dpsi/dx, 4 - dpsi/dy)
      type(TAreaValue), allocatable :: arval_funp(:) !(nu) !�������� ������� � ������ ������������� ��� ���������� ������� ��������� �������� ��� �����������
                                              !������ ������ - ����� ������� ������� � ������� (��� ��������� 4 ������� ��� ������)
      !����������� (���������������� ������ ��� ����� � ������������ � �������)
      integer(4) areaEq                    !��� ��������� (��� �������)
      integer(4), allocatable :: areaind(:,:) !(ntr,umaxtr) !������ ����������� (����� �������) � ������� � ������� ���������      
      real(8), allocatable :: psiarea(:,:) !(ntr,umaxtr) !�������� ������� ������� � �������
	                                       !������ ������=1,2
										   !1 - ������� ������� � ���������� ���� �����������
										   !1,2 - df/dx,df/dy � ��������� ��������
      integer(4), allocatable :: eq_ind(:) !(umaxtr) ������ �������, ��� ������� �������� ������� ��������� (TArea_Collocate_points%iu)
                       
      !��� ����������
      type(TCash) self_cash !����������� ���
      type(TCash), pointer :: cash !������ �� ��� ������ ������� ��� ����������� ���
      
      !����� ���������� 
		  type(TArea_Collocate_points) cppa  !����� ���������� - ����� ����������� ��������� � �������
    endtype
    
    type fict_var !������ ��� ��������� ��������� �����������
      integer mode !1 - ��� gu_gen, 2 - ��� gu_genGlob
      type(TBoungline_GU_gen), pointer :: gu_gen !����� ��������� ������� ��� ��������� �����������, ��� ��� �������� ����� �������� �����������
      type(TBoungline_GU_genGlobal), pointer :: gu_genGlob !����� ��������� ������� ��� ��������� �����������, ��� ��� �������� ����� �������� �����������
      integer(4), allocatable :: psiind(:) !(gu_gen%n) ������ �����������
      real(8), allocatable :: psiom(:)     !(gu_gen%n) - �������� �������, ���� ��� ������
      integer(4) i !������ ��� ������� ���������� ������������� (��������� � ������� ������ �� �������)
    endtype
  
  type col_info
    integer(4) ncol  !���������� ����������� � ������� ��� �����
    integer(4) icol_area !������ ���������� �������� � ������� col, ������� ����� ����� ������
    integer(4), allocatable :: col(:) !(ncol)
    integer(4) ncol_fict
    integer(4), allocatable :: col_fict(:) !(ncol_fict)
  end type
    
  type matrix_mb
	  integer(4) ix !����� ������� ������� ����������� � ���������� �������
    integer(4) ix2 !����� ���������� ������� ����������� � ���������� �������
    integer(4) iu !����� ������� ��������� � ���������� �������
    integer(4) iu2 !����� ���������� ��������� � ���������� �������
    integer(4) nx !����� ����������� (��������)
    integer(4) nu !����� ��������� (�����) ������������ ��� ���������� �� ��������� � pg_get_matrix
    integer(4) nu_all !����� ���������, ������� ������ ����
	  integer(4) nub !���������� ��������� �� ������� 
    integer(4) area_sm_count  !���������� ����������� � ������� (��� matrix_main%sparse%m)
    type(col_info), pointer :: ci
    type(col_info) self_ci
    integer(4) nnt !������� ����������� ��������� � ���� �������
    integer(4) nu_add !����� �������������� ��������� ��� ���������������� �������
  endtype
  
  type TOMP_buffer
    real(8), allocatable :: eq(:) !gs%m%nx_all
    real(8) b
    integer(4) i  !����� �������� ���������
    integer(4), allocatable :: buffc(:) !(nx) ������ ��� ��������
  endtype
  
  type TOMP_buffers
    integer(4) max_threads
    type(TOMP_buffer), allocatable :: b(:) !(max_threads)
  endtype
    
	type matrix_main
      integer(4) matrix_store_type  !0 - ������� ���������� �������
                                    !1 - ����������� ���������� �������
                                    !2 - ������-������������ ������-����������� �������
      real(8), allocatable :: m(:,:)
      type(sparse_matrix) sparse
      real(8), allocatable :: b(:)
      integer(4) nx !����� ����������� (��������)
      integer(4) nu !����� ��������� (�����)
      integer(4) ix_fict !����� �������, � �������� ����� ���������� ��������� ����������
      integer(4) nx_fict !����� ��������� ����������
      integer(4) nx_all  !����� ����� ����������� c ���������
      type(fict_var), allocatable :: fvar(:) !(nx_fict) ��������� ����������
      real(8) norm     !����� ������� ||Ax-b|| ��� ���������������� �������
      logical find_norm !����� ����� �������
      real(8) sparse_percent !������� ��������� (��� ������������ ���������) �� ���� ���������
      real(8) memory_mb !������ � ���������� ��������, ����������� ��� �������� �������
      integer(4) allocated_matrix_store_type !-1 - ���� �� ���������������� ��� =matrix_store_type, ���� ����������������
      integer(4) count_closing  !���������� ���������, ���������� ��� ���������� ���������, ���� 0, �� �������������� ��� ���������� (=nx)
      integer(4) nc_global          !���������� ��������� ������������� (��� ��������������� ��������� ������ ��� sparse)
      type(main_block_matrix) bm   !������-������������ ������� ��� matrix_store_type=2
      logical have_initial_approxsolv
      integer(4) max_threads
      integer(4) max_ncol
      type(TOMP_buffers) omp_buff
      real(8), allocatable :: res(:)  !������� ����
    endtype
    
	type TConstArea
    real(8) k_helm     !����������� � ������ ����������� � ���������� �������������
		real(8) k_oss      !����������� � ��������������� ������ (=1 ��� ��������� �������, =-1 ��� ��������� ��� ������� ����)
    real(8) Re         !����� ���������� ��� ��������� �����-������
    real(8) k          !�������������
    real(8) mub        !�������� ���������
	endtype

    type TArea !�������
	    integer(4) i !�����
      integer(4) mode !������� �������, ���������� �������������
	    integer(4), allocatable :: type_eq(:)  !(nu) ���� ��������� � �������
          !1 - ������������� �������
		      !2 - ������������� ������� (������ ��������� � ��������������� ���������)
          !3 - ��������������� �������
          !4 - ��������� ��������
          !5 - ���������� ��������� ����������� (����� ����������� � �������)
          !6 - ������������ ��������� ����������� (����� ����������� � �������)
          !7 - ��������� �������� (������ ��������� � ������������ ��������������� ���������)
          !8 - ������� ���� \Delta\psi=v(x,y)\omega  � ��������� �� \omega. ����������� ������ �� ������� � ��� \omega � ������� 
          !9 - ���������� ��������� ����������� ����� ������� ������� \Delta\psi+k_helm^2\psi=0
          !10 - ���������� ��������� ����������� ����� ������� ������� \Delta\psi-k_helm^2\psi=0
          !11 - ������� ���� \Delta\psi=v(x,y)\omega  � ��������� �� \omega. ����������� ������ �� �������
                !��� ��������� �� \omega ������������ ��������� type_eq(2)
          !12 - ������ ���������� ��������� ����������� ����� ������� ������� \Delta\psi+k_helm^2\psi=0 (� �������, ����� � ������� �������� ��������� ���������)
          !13 - ������ ���������� ��������� ����������� ����� ������� ������� \Delta\psi-k_helm^2\psi=0 (� �������, ����� � ������� �������� ��������� ���������)
          !14 - ������ ���������� ��������� ����������� (����� ����������� � �������)
		      !15 - ��������� ��������� \Delta^2\psi-k_helm^2\Delta\psi=0 ����� ������ ������� �����, ���������� ������� �������
		      !16 - ���������� ��������� �������� \Delta C = vx*dC/dx + vy*dC/dy
		      !17 - ������������ ��������� �������� \Delta C = vx*dC/dx + vy*dC/dy + f
		      !18 - ������������ ��������� ��� ���������������� ������ \Delta\psi = vy*d\psi/dy + f
		      !19 - ��������� ������� ��� ���������������� ������ \Delta\psi = -(k_oss/y)*d\psi/dy  (k_oss=1 - ��-�� �������, k_oss=-1 - ��-�� ��� ������� ����)
		      !20 - ��������� ������� ��� ���������������� ������ \Delta\psi = -(k_oss/y)*d\psi/dy  (k_oss=1 - ��-�� �������, k_oss=-1 - ��-�� ��� ������� ����) - ������ �������������� G_1/y
          !21 - ������������ ��������������� ��������� (����������� ������ �� �������) ��� ��������� �����-������
          !22 - ������������ ��������������� ��������� � ������������ � ������ �����
                !\Delta^2\psi=f1+f2\psi+f3\eta+f4 d\psi/dx+f5 d\psi/dy+f6 d\eta/dx+f7 d\eta/dy
          !23 - ������ ��������� ��� (22) ������������� ���������������� ��������� � ������������ � ������ �����
          !24 - ��������� �����-������ - ������������ ��������������� � ������������ d\psi/dx,d\psi/dy
          !25 - ��������� �����-������ - ������������ ��������������� � ������������ d\om/dx,d\om/dy
          !26 - ��������� �������� � ������������ � ������ �����
                !\Delta\psi=f1+f2\psi+f3 d\psi/dx+f4 d\psi/dy
      integer(4) n_eq_var !���������� ����������� � ������ ����� � ���������� ���� 22,24,25,26
      logical, allocatable :: eq_var(:) !(n_eq_var) true - ����������� ������������
      integer(4), allocatable :: var_ind(:) !(n_eq_var) ������ ������ ��� ����������� � areaind
      integer(4) nb  !����� ������ (��������� �������)
      type(TBound), allocatable :: bnd(:)
      integer(4) nu ! =1,2 ������� ��������� ��� ������� (����� ��������� ������� �������)
      integer(4) umax !=nu*2 - ���������� ����������� � ������ ����� �������
		  integer(4) umaxtr !���������� ����������� � ������ ������������
      type(matrix_mb) m  !�������
      type(areatype) a   !������� � ��������������
		  logical haveAreaEq
		  type(TConstArea) const   !���������

      !����� ���������� 
		  type(TBound_Collocate_points) cpp  !����� ���������� - ����� ����������� ��������� �� �������

		  !real(8), pointer :: func_inf
		  !pointer (pfunc_inf,func_inf)
      real(8) cash_size_mb
      logical need_cash_bb(0:max_int),need_cash_ab(0:max_int),need_cash_aa(0:max_int),need_cash_ba(0:max_int)
      !���������� � �������, ���������� ������������ ������ �������
      real(8) dx,dy !��������
      real(8) sina,cosa !�������
      type(TArea), pointer :: a_ref !������ �� �������, �� ������� ���� ����������� ��������� (�� ���� �� ���������� ��������, �� �������!!!)
      integer(4) type_rotate !��� ����������� �������
                       !-1 - ��� ������
                       !0 - ������ �������
                       !1 - ������� �� 180 ����
                       !2 - ������� �� 90
                       !3 - ������� �� 270
                       !4 - ������������ �������
      logical oss_resolve !������������� ��� ������� ��������
      integer(4) iorder !����� ������� ��� ����������� ����� ��� �������������� ������ �������� � ������-������������ �������
    endtype
    
	type TConstGlobal
		real(8) ds0numint  !��� ���������� ��������������
    logical use_numerical_int   !������������ ��������� ��������������
		logical test_izlom !��������� �� ����� ����� ������� ��� ���������� �������� ������� �� ������
		real(8) angle_izlom !������������ ���� ������ � ��������
		integer(4) di !������� ����� ����� �� ����� ��� �����������, �.�. ��� ���������� ������� ������������
    real(8) epsMinDistBound !����������� ���������� �� ����� �� �������
    real(8) epsMinDistBound2 !����������� ���������� �� ����� �� ������� (���� ��������������� ������� � ���� eps, �� ���������� ����� ������������)
    real(8) epsMinDistBound3 !����������� ������������� ���������� �� ����� �� ������� (��� ����������� in_bound=2)
  endtype
  
  type Tsndns_type
    integer(4) n !����� ����� ����������
    real(8) shift !����� ������� �������� ������ �����
    logical use_shift2 !true - ����� ���������� ��� ������� ��������� �������� ������������ �������
    real(8), allocatable :: s_ndns(:,:) !������� �������� ����� ���������� ����������� �� ��� boundLineType=2
                                        !������ ������ - ���������
                                        !������ ������ - ����� ����������
  endtype
    
    type TGreenSolve   !��������� ����� - ���� ������
	    integer(4) i !�����
        integer(4) na     !���������� ��������
        type(TArea), allocatable :: a(:)  !�������
        type(matrix_main) m  !����� �������
        type(TConstGlobal) const   !���������
        real(8) cash_size_mb
        integer(4) constvaln 
        integer(4) constvalk 
		    integer(4), allocatable :: constvalinds(:)
		    real(8), allocatable :: constvals(:)
        integer(4), allocatable :: constvala(:) !������� ��������, � �������� ������ ���� ��������� ���������� (��� matrix_store_type=2)
                                               !=0 ���� �� ���������������� 
        integer(4), allocatable :: constvalinds_fict(:) !��������� (�������������) ������� ��� �����������, ������� ����������� � ����� �������, � ������ ���� �������� � ������ 
        integer(4) constvalinds_fict_i !������� ��� ��������� ��������
    endtype
    
    type TGreenSolves   !��������� �����
	    integer(4) ngs  !����� �����
      type(TGreenSolve), allocatable :: ggs(:)  !������
      integer(4) nsn !����� ��������� sn
      type(Tsndns_type), allocatable :: sn(:) !������� �������� ����� ���������� ����������� �� ��� boundLineType=2
    endtype
	
    type(TGreenSolves), target :: gsMain
	type(TGreenSolve), pointer :: gs  !������� ������
	type(TArea), pointer :: gsarea    !������� ������� ������� ������
  type(areapart), pointer :: gsareapart    !������� ������� �������
	type(TBound), pointer :: gsbnd    !������� �������    
	type(TConstArea), pointer :: gsareaconst    !��������� ������� �������
	type(TBoundline), pointer :: gsbndl    !������� ������� �������
	type(TBoundline2), pointer :: gsbndl2    !������� ������� �������

	!��� ������ �� �����
  logical gs_use_mkl_progress
	real(8) gs_time0, gs_time_interval, gs_time1
    integer gs_use_parallel_matrix_build
    integer gs_use_parallel_get_fun
    integer gs_use_parallel_num_int
    logical gs_stop_on_int_dxdy_inbound
    logical gs_use_cash
    logical gs_cash_presolve              !��������������� ���������� ���� ����� ������������ �������
    logical gs_use_cash_lock_get
    logical gs_write_matrix
    logical gs_rebuild_matrix_for_iterate
    logical first_matrix_build
    logical use_global_sparsem   !������������ ����� ���������� ����������� �������
    logical max_equation_coef_eq_1 !��������� ������������ ������� � ������������� �� ������ �������� =1
    real(8) equation_coef_eq_0_eps !�������� �� �������� ���� ��� ���������� ������������� ��������� � ������� (=d0 - �� ���������)
    real(8) gs_block_matrix_solver_eps  !eps ��� �������� ������-������������ ������� (calc_slau_block_diagonal)
    integer(4) gs_block_matrix_solver_maxIter  !������������ ����� �������� ��� �������� ������-������������ ������� (calc_slau_block_diagonal)
    integer(4) gs_block_matrix_solver_maxIter_convergence !����������� ���������� ��������, �� ������� �� ����������� ���������� ����������� (calc_slau_block_diagonal)
    real(8) gs_block_matrix_solver_lambda  !�������� ���������� x=lam*x[i]+(1-lam)*x[i-1] (calc_slau_block_diagonal)
    integer(4) gs_pardiso_eps !(��. slau.f90 pardiso_eps) 
    integer(4) gs_system_count,gs_system_i !���������� � ������� �������� ���� (��� mkl_progreess)
    logical gs_write_log_file
    logical :: gs_used_pg_start=.false.
    logical :: gs_cash_is_solve !���������� ����������, ������������ ��� ���������� ����
    logical gs_profiling        !��� �������� ���������� ���� ������������ (����� ����������� - ������� �����)
    logical gs_test_areaval_eq_0 !��������� areaval ����� �� �������� ������ �����������
    logical gs_use_dual_integral_solve !��������� ��� ���������, ���������� Re � Im ������ ������������ ���������, �� ���� ��� ��� ��������������� ���������� ����
    logical gs_test_point_collenear_for_BesselG  !������������� ���������� ���������� � ����������� ��������� (��� ����� ������� �� ����� ������ � �.�.)
    logical gs_test_point_near_bound !���������� ������� ������ ������� �� ���������� ������ ����� ������
    
    
    interface
   
    subroutine init_bound_cash(cash,nb)
    import :: TCash
    integer(4) nb
    type(TCash), target :: cash
    end
    
    subroutine init_col_info_closing(a,ci,b1,bl1,k1,b2,bl2,k2)
    import :: TArea,col_info,TBound,TBoundline
    Type(TArea) a
    type(col_info), target :: ci
    Type(TBound) b1,b2
    type(TBoundline) bl1,bl2
    integer(4) k1(:),k2(:)
    end
    
    subroutine init_boundline_geomlineds3_array_mode1(ds1,ds2,kk,mode,x1,y1,x2,y2,min_np,vx,vy,npanel,npanel_k)
    import :: TAreaValue
    integer(4) npanel !���������� ������� � ������� �������
    integer(4) npanel_k  !���������� ������� �� ����������� �������
    real(8) x1,y1,x2,y2  !���������� ������ �������
    real(8) ds1          !��� � ������ �������
    real(8) ds2          !��� � ����� (��� � ��������) �������
    real(8) kk           !����������� ��������� ����� (>1!!!)
    integer(4) mode      !��� ��������
    					 !1 - -- --- ---- ---- ----
    					 !3 - -- --- --- --- --- -- - (���������� ��������������� ������ ����� �������)
    integer(4) min_np    !����������� ����� �������
    					 !���� ���������� ����� �������,������ ������������, 
    					 !�� �������� ����������� ��������� � ����������� ������ �������
    type(TAreaValue), target :: vx,vy
    end
    
    recursive subroutine init_areapart_geom_ref(aref,aa,bb,cc)
    !���������������� ������� ������� �� ������ ������� � ��������������� ���������
    !z1=A*z+B*conj(z)+C
    import :: areapart
    complex(8) aa,bb,cc
    type(areapart), target :: aref
    end
    
    subroutine init_boundline_geomlineds2_array(ds1,ds2,x1,y1,x2,y2,vx,vy,npanel,min_np)
    !���������������� ��������� ��� ������� ������� � ���� ������ ����� � ���������� ds
    !���������� �� ds1 �� ds2
    import :: TAreaValue
    integer(4) npanel !���������� ������� � ������� �������
    real(8) x1,y1,x2,y2  !���������� ������ �������
    real(8) ds1          !��� � ������ �������
    real(8) ds2          !��� � ����� �������
    integer(4) min_np    !����������� ����� ������� (���� 0 - �� �����������)
    type(TAreaValue), target :: vx,vy
    end
    
    function get_bb_cash_integral_2(a,knd,j,knd0,i,nf) result(res)
    !��� �������������� �� ������ ����������� ����� �� ������
    !���������� ������� ��� ��������
    Import :: TArea
    !j - ����� ��� ������, �� ������� ���� ��������������
    !i - ����� ����������� ����� (������ � ������� ���������)
    !nf - ����� ���������
    !knd0 - ������ ������� � ������� � ����������� ������
    !a   - ������� ��� ��������� ������ ��������������
    !knd - ����� �������, ��� ��������� ������ ��������������
    integer(4) knd,j,knd0,i,nf
    real(8) res
    !DEC$ IF DEFINED (DEBUG)
    type(TArea), target :: a
    !DEC$ ELSE
    type(TArea) a
    !DEC$ ENDIF
    end
    
    function get_ba_cash_integral_2(a,knd,j,i,nf) result (res)
    !��� �������������� �� ������ ����������� ����� � ������������
    !���������� ������� ��� ��������
    Import :: TArea
    !j - ����� ��� ������, �� ������� ���� ��������������
    !i - ����� ����������� ����� (������������)
    !nf - ����� ���������
    !a   - ������� ��� ��������� ������ ��������������
    !knd - ����� �������, ��� ��������� ������ ��������������
    integer(4) knd,j,i,nf
    real(8) res
    !DEC$ IF DEFINED (DEBUG)
    type(TArea), target :: a
    !DEC$ ELSE
    type(TArea) a
    !DEC$ ENDIF
    end
    
    function get_ab_cash_integral_2(a,j,knd0,i,nf) result(res)
    !��� �������������� �� ������������ ����������� ����� �� ������
    !���������� ������� ��� ��������
    Import :: TArea
    !j - ����� ��� ������������, �� ������� ���� ��������������
    !i - ����� ����������� ����� (������ � ������� ���������)
    !nf - ����� ���������
    !knd0 - ����� ������� � ����������� ������
    !a - �������
    integer(4) j,knd0,i,nf
    real(8) res
    !DEC$ IF DEFINED (DEBUG)
    type(TArea), target :: a
    !DEC$ ELSE
    type(TArea) a
    !DEC$ ENDIF
    end
    
    function get_aa_cash_integral_2(a,j,i,nf) result (res)
    !��� �������������� �� ������������ ����������� ����� �� ������������
    !���������� ������� ��� ��������
    Import :: TArea
    !j - ����� ��� ������������, �� ������� ���� ��������������
    !i - ����� ����������� ����� (������������)
    !nf - ����� ���������
    !a - �������
    integer(4) j,i,nf
    real(8) res
    !DEC$ IF DEFINED (DEBUG)
    type(TArea), target :: a
    !DEC$ ELSE
    type(TArea) a
    !DEC$ ENDIF
    end
    
    subroutine gaintf_prepare(nn,x1,y1,x2,y2,ds,xx,yy,mode,ifs)
    Import :: intf_struct
    integer(4) nn,j
    real(8) x1,y1,x2,y2,ds
    real(8),target :: xx(nn),yy(nn)
    integer(4) mode !1 - ������ ����� x1,y1,x2,y2,ds
                    !2 - �������������� ����� x,y1,n
                    !3 - ����� �������� ��������� x,y,n
    type(intf_struct), target :: ifs
    end
    
    subroutine init_boundline_gu_gen_var(gu,i,var)
    Import :: TBoundline_GU
    type(TBoundline_GU), target :: gu
    integer(4) i    !����� �������, =0 ��� c0
    real(8) var(gu%bndl%npanel)   !������ �������������
    end
    
    subroutine init_boundline_gu_genGlobal_var(gu,i,ind,var)
    Import :: TBoundline_GU
    type(TBoundline_GU), target :: gu
    integer(4) i    !����� �������, =0 ��� c0
    integer(4) ind  !������ ���������� ���������� � gs%constvalinds
    real(8) var(gu%bndl%npanel)   !������ �������������
    end
    
    subroutine init_boundline_gu_genGlobal_const(gu,i,ind,varconst)
    Import :: TBoundline_GU
    type(TBoundline_GU), target :: gu
    integer(4) i    !����� �������, =0 ��� c0
    integer(4) ind  !������ ���������� ���������� � gs%constvalinds
    real(8) varconst     !��������� ��� ������ �������������
    end
    
    end interface
        
    end

subroutine pg_start
!dec$ attributes dllexport:: pg_start
!������������� ���������� ��������
use pgmod
gs_use_mkl_progress=.false.
gs_time_interval=d1
gs_use_parallel_matrix_build=1 !1
gs_use_parallel_get_fun=1      !1
gs_use_parallel_num_int=0
gs_use_cash=.false.
gs_cash_presolve=.false.
gs_use_cash_lock_get=.true.
gs_write_matrix=.false.
gs_rebuild_matrix_for_iterate=.true.
first_matrix_build=.true.
gs_stop_on_int_dxdy_inbound=.true.
use_global_sparsem=.false.
max_equation_coef_eq_1=.false.
equation_coef_eq_0_eps=1.0d-16  !=1.0d-16 �������� double
gs_pardiso_eps=6   
gs_block_matrix_solver_eps=1.0d-12
gs_block_matrix_solver_maxIter=200
gs_block_matrix_solver_lambda=0.5d0 
gs_block_matrix_solver_maxIter_convergence=20
gs_write_log_file=.false.
gs_used_pg_start=.true.
gs_profiling=.false.
gs_test_areaval_eq_0=.false.
gs_use_dual_integral_solve=.true.
gs_test_point_collenear_for_BesselG=.true.
gs_test_point_near_bound=.true.
call sb_init
end

subroutine pg_finish
!dec$ attributes dllexport:: pg_finish
!����������
call pg_deallocate_mem 
call pg_stop_write_log_file
!DEC$ IF DEFINED (DEBUG)
call system("pause")
!DEC$ ENDIF
end

subroutine pg_start_write_log_file
!dec$ attributes dllexport:: pg_start_write_log_file
use pgmod
gs_write_log_file=.true.
OPEN (gs_log_file_i,FILE='pg.log')
end

subroutine pg_stop_write_log_file
!dec$ attributes dllexport:: pg_stop_write_log_file
use pgmod
if (gs_write_log_file) then
  close (gs_log_file_i)
  gs_write_log_file=.false.
endif
end

subroutine pg_allocate_problems(n)
!dec$ attributes dllexport:: pg_allocate_problems
!���������������� ���������� �����
use pgmod
integer(4) n !���������� �����
integer(4) i
if (.not.gs_used_pg_start) call pg_start
gsMain%ngs=n
gsMain%nsn=0
allocate(gsMain%ggs(n))
do i=1,n
  call pg_bind_problem(i)
  call initconst_GreenSolve
enddo
end

function pg_get_next_constind
!dec$ attributes dllexport:: pg_get_next_constind
!�������� ��������� ������ � ������� ���������� ��������
use pgmod
integer(4) pg_get_next_constind
gs%constvalk=gs%constvalk+1
pg_get_next_constind=gs%constvalk
end

subroutine pg_set_constvala(ia,i1,i2)
!dec$ attributes dllexport:: pg_set_constvala
!���������������� ������� �������� �� ���������
use pgmod
integer(4) ia !����� �������
integer(4) i1,i2 !�������� � ������� gs%constvala
gs%constvala(i1:i2)=ia
end

subroutine pg_allocate_constvalind(n)
!dec$ attributes dllexport:: pg_allocate_constvalind
!���������������� ����������, ���������� �� ���������� �������
use pgmod
integer(4) n !���������� �����������
gs%constvaln=n
gs%constvalk=0
allocate(gs%constvalinds(n))
allocate(gs%constvals(n))
allocate(gs%constvala(n))
allocate(gs%constvalinds_fict(n))
gs%constvalinds=0
gs%constvals=d0
gs%constvala=0
gs%constvalinds_fict=0
gs%constvalinds_fict_i=0
end

subroutine pg_bind_problem(k)
!dec$ attributes dllexport:: pg_bind_problem
!������� ������� ������
use pgmod
integer(4) k  !����� ������� ������
gs=>gsMain%ggs(k)
gs%i=k
end

subroutine pg_bind_domain(k)
!dec$ attributes dllexport:: pg_bind_domain
!������� ������� �������
use pgmod
integer(4) k  !����� ������� � ������� ������
gsarea=>gs%a(k)
gsarea%i=k
end

subroutine pg_bind_areapart(k)
!dec$ attributes dllexport:: pg_bind_areapart
!������� ������� ������� �������
use pgmod
integer(4) k  !����� ������� � ������� �������
gsareapart=>gsarea%a%part(k)
gsareapart%i=k
end

subroutine pg_bind_bound(k) 
!dec$ attributes dllexport:: pg_bind_bound
!������� ������� �������
use pgmod
integer(4) k  !����� ������� � ������� �������
gsbnd=>gsarea%bnd(k)
gsbnd%i=k
end

subroutine pg_bind_boundline(k) 
!dec$ attributes dllexport:: pg_bind_boundline
!������� ������� ������� �������
use pgmod
integer(4) k  !����� ������� � ������� �������
if (gsbnd%boundLineType==1) then
  gsbndl=>gsbnd%line(k)
  gsbndl%i=k
elseif (gsbnd%boundLineType==2) then
  gsbndl2=>gsbnd%line2(k)
  gsbndl2%i=k
endif
end

subroutine bind_AreaConst(k)
!������� ������� �������
use pgmod
integer(4) k  !����� ������� � ������� ������ ������
gsareaconst=>gs%a(k)%const
end

subroutine initconst_GreenSolve
!���������������� ���������
use pgmod
gs%const%ds0numint=ds0numint
gs%const%use_numerical_int=.false.
gs%const%test_izlom=.true.
gs%const%angle_izlom=45
gs%const%epsMinDistBound=1.0d-2
gs%const%epsMinDistBound2=1.0d-6
gs%const%epsMinDistBound3=1.0d-4
gs%const%di=2 
gs%na=0
gs%m%norm=d0
gs%m%find_norm=.false.
gs%m%count_closing=0
gs%cash_size_mb=d0
gs%m%omp_buff%max_threads=0
end

subroutine pg_set_matrix_n_add(n)
!dec$ attributes dllexport:: pg_set_matrix_n_add
!������� ���������� �������������� ��������� ��� ���������������� ����
use pgmod
integer(4) n
gsarea%m%nu_add=n
end

subroutine pg_set_matrix_sparse_count_closing(n)
!dec$ attributes dllexport:: pg_set_matrix_sparse_count_closing
!������� ���������� ���������, ���������� � ����������� ���� ��� ���������� ���������
use pgmod
integer(4) n
gs%m%count_closing=n+1  !+1 ��� ����� ������������� ��������
end

subroutine pg_set_areaconst(j,val)
!dec$ attributes dllexport:: pg_set_areaconst
!������� �������� �������
use pgmod
integer(4) j !��� ���������
                    !1-k_helm
					!2-k_oss
          !3-Re
real(8) val !�������� ���������
selectcase (j) 
case (1)
  gsarea%const%k_helm=val
case (2)
  gsarea%const%k_oss=val
case (3)
  gsarea%const%Re=val
case (4)
  gsarea%const%k=val
case (5)
  gsarea%const%mub=val
endselect
end

subroutine pg_allocate_domains(n)
!dec$ attributes dllexport:: pg_allocate_domains
!���������������� ���������� �������� � ������ � ������� ������ ������
use pgmod
integer(4) n !���������� ��������
integer(4) i
gs%na=n
allocate(gs%a(n))
do i=1,n
    call nullarea(i)
enddo
gs%m%matrix_store_type=0
gs%m%allocated_matrix_store_type=-1
gs%constvaln=0
end

subroutine pg_deallocate_mem
!dec$ attributes dllexport:: pg_deallocate_mem
!���������� ������ �� ���� �������
use pgmod
integer(4) i
if (allocated(gsMain%ggs)) then
  do i=1,gsMain%ngs
    call pg_bind_problem(i)
    call deallocate_domains
  enddo
  deallocate(gsMain%ggs)
endif
if (allocated(gsMain%sn)) then
  do i=1,gsMain%nsn
    if (allocated(gsMain%sn(i)%s_ndns)) deallocate(gsMain%sn(i)%s_ndns)
  enddo
  deallocate(gsMain%sn)
endif
end

subroutine deallocate_domains
!���������� ������ � ������
use pgmod
integer(4) i
call deallocate_mm(.true.) !������ ���� ����� ���������, ����� ���������� ������� ��������
if (allocated(gs%a)) then
    do i=1,gs%na
        call deallocate_bounds(i)
    enddo
    deallocate(gs%a)
endif
call deallocate_domains_constval
end

subroutine deallocate_domains_constval
!���������� ������ gs%constval � ������
use pgmod
if (allocated(gs%constvalinds)) deallocate(gs%constvalinds)
if (allocated(gs%constvals)) deallocate(gs%constvals)
if (allocated(gs%constvala)) deallocate(gs%constvala)
if (allocated(gs%constvalinds_fict)) deallocate(gs%constvalinds_fict)
end

subroutine deallocate_col_info(ci)
use pgmod
type(col_info) ci
if (allocated(ci%col)) deallocate(ci%col)
if (allocated(ci%col_fict)) deallocate(ci%col_fict)
end

subroutine init_omp_buffer
use pgmod
integer(4) omp_get_max_threads,i
type(TOMP_buffer), pointer :: b
gs%m%omp_buff%max_threads=omp_get_max_threads()
allocate(gs%m%omp_buff%b(0:gs%m%omp_buff%max_threads-1))
do i=0,gs%m%omp_buff%max_threads-1
  b=>gs%m%omp_buff%b(i)
  allocate(b%eq(gs%m%nx_all))
  b%eq=d0
  if (gs%m%matrix_store_type==1.or.gs%m%matrix_store_type==2) allocate(b%buffc(gs%m%nx))
enddo
end

subroutine free_omp_buffer
use pgmod
integer(4) i
type(TOMP_buffer), pointer :: b
if (allocated(gs%m%omp_buff%b)) then
  do i=0,gs%m%omp_buff%max_threads-1
    b=>gs%m%omp_buff%b(i)
    if (allocated(b%eq))deallocate(b%eq)
    if (allocated(b%buffc))deallocate(b%buffc)
  enddo
  deallocate(gs%m%omp_buff%b)
  gs%m%omp_buff%max_threads=0
endif
end

subroutine deallocate_mm(need_mm)
!������������ ������ ���������� �������
use pgmod
type(matrix_main), pointer :: m
logical need_mm
integer(4) i
m=>gs%m
if (need_mm) then
  if (allocated(m%m)) deallocate(m%m)
  if (allocated(m%b)) deallocate(m%b)
  call bm_deallocate_sparse_mrows(m%sparse)
endif
call free_omp_buffer
call bm_deallocate(need_mm,m%bm)
call bm_deallocate_sparse_m(m%sparse)
if (allocated(m%fvar)) then
  do i=1,UBOUND(m%fvar, 1) !m%nx_fict
    if (allocated(m%fvar(i)%psiind)) deallocate(m%fvar(i)%psiind)
    if (allocated(m%fvar(i)%psiom)) deallocate(m%fvar(i)%psiom)
  enddo
  deallocate(m%fvar)
endif
do i=1,gs%na
  call deallocate_col_info(gs%a(i)%m%self_ci)
enddo
end

subroutine nullarea(ia)
!�������� ������, ��������� � ��������
use pgmod
integer(4) ia !����� �������
type(TArea), pointer :: a
a=>gs%a(ia)
a%type_eq=0
a%n_eq_var=0
a%nb=0
a%nu=0
a%umax=0
a%umaxtr=0
a%a%ntr=0
a%a%n=0
a%a%npart=0
a%a%geom_inited=.false.
a%haveAreaEq=.false.
a%a%self_cash%bnd_inited=.false.
a%a%self_cash%area_inited=.false.
a%a%cash=>a%a%self_cash
a%cpp%inited=.false.
a%a%cppa%inited=.false.
a%cash_size_mb=d0
a%need_cash_bb=.false.
a%need_cash_ab=.false.
a%need_cash_aa=.false.
a%need_cash_ba=.false.
a%m%nu_add=0
a%dx=d0
a%dy=d0
a%sina=d0
a%cosa=d0
a%a_ref=>null()
a%type_rotate=-1
a%oss_resolve=.false.
a%iorder=0
a%mode=0
end

subroutine pg_set_domain_equation(type_eq)
!dec$ attributes dllexport:: pg_set_domain_equation
!������ ��� ��������� ��� �������
use pgmod
integer(4) type_eq  !��� ������ � ������� �������
integer(4) type_eq2  !��� ������ ������ � ������� �������
if (type_eq==3) then
  type_eq2=2
elseif (type_eq==15) then
  type_eq2=13
elseif (type_eq==21) then
  type_eq2=7
elseif (type_eq==22.or.type_eq==24.or.type_eq==25) then
  type_eq2=23
else
  type_eq2=0
endif
call set_domain_equation(type_eq,type_eq2)
end

subroutine pg_set_domain_equation_syst(type_eq,type_eq2)
!dec$ attributes dllexport:: pg_set_domain_equation_syst
!������ ��� ��������� ��� ������� (c������ ���������)
use pgmod
integer(4) type_eq  !��� ������ � ������� �������
integer(4) type_eq2  !��� ������ ������ � ������� �������
if (type_eq.ne.8.and.type_eq.ne.11) call gs_print_stop('!!! type_eq=8,11')
if (type_eq==8) then
  if (type_eq2.ne.2.and.type_eq2.ne.7.and.type_eq2.ne.12.and.type_eq2.ne.13.and.type_eq2.ne.14) call gs_print_stop('!!! type_eq=2,7,12,13,14')
elseif (type_eq==11) then 
  if (type_eq2.ne.2.and.type_eq2.ne.12.and.type_eq2.ne.13) call gs_print_stop('!!! type_eq=2,12,13')
endif
call set_domain_equation(type_eq,type_eq2)
end 

subroutine set_domain_equation(type_eq,type_eq2)
!������ ��� ��������� ��� �������
use pgmod
integer(4) type_eq  !��� ������ � ������� �������
integer(4) type_eq2  !��� ������ ������ � ������� �������
logical eq_var(7)
gsarea%nu=1
if (type_eq==3.or.type_eq==8.or.type_eq==11.or.type_eq==15.or.type_eq==21.or.type_eq==22.or.type_eq==24.or.type_eq==25) gsarea%nu=2
gsarea%umax=gsarea%nu*2
allocate(gsarea%type_eq(gsarea%nu))
if (type_eq==22.or.type_eq==24.or.type_eq==25) gsarea%n_eq_var=7
if (type_eq==26) gsarea%n_eq_var=4
if (gsarea%n_eq_var>0) then
  allocate(gsarea%eq_var(gsarea%n_eq_var))
  allocate(gsarea%var_ind(gsarea%n_eq_var))
  gsarea%eq_var=.false.
  gsarea%var_ind=0
endif
gsarea%type_eq(1)=type_eq
if (gsarea%nu==2) gsarea%type_eq(2)=type_eq2
call init_mesh_gu
if (type_eq==24.or.type_eq==25) then
  eq_var=.false.
  if (type_eq==24) then
    eq_var(4:5)=.true.
  else
    eq_var(6:7)=.true.
  endif
  call pg_set_eq_var(eq_var,7)
endif
end

subroutine pg_set_eq_var(eq_var,n)
!dec$ attributes dllexport:: pg_set_eq_var
!������ ������ gsarea%eq_var
use pgmod
integer(4) n !���������� ��������� � eq_var
logical eq_var(n) 
integer(4) i,k
if (n<gsarea%n_eq_var) call gs_print_stop("Error pg_set_eq_var! (n<gsarea%n_eq_var)")
gsarea%eq_var(1:gsarea%n_eq_var)=eq_var(1:gsarea%n_eq_var)
k=0
gsarea%var_ind=0
do i=2,gsarea%n_eq_var
  if (gsarea%eq_var(i)) then
    k=k+1
    gsarea%var_ind(i)=k
  endif
enddo
end

subroutine pg_allocate_bounds(inb)
!dec$ attributes dllexport:: pg_allocate_bounds
!�������� ������ ��� �������
use pgmod
integer(4) inb  !���������� ������ � ������� �������
integer(4) i
gsarea%nb=inb
allocate(gsarea%bnd(inb))
do i=1,inb
  call nullbound(i)
enddo
end

subroutine deallocate_bounds(ia)
!���������� ������ ��� �������
use pgmod
integer(4) ia !����� �������
integer(4) i
type(TArea), pointer :: a
a=>gs%a(ia)
if (allocated(a%bnd)) then
    do i=1,a%nb
        call deallocate_boundlines(ia,i)
    enddo
    deallocate(a%bnd)    
endif
if (allocated(a%type_eq)) deallocate(a%type_eq)
if (allocated(a%eq_var)) deallocate(a%eq_var)
if (allocated(a%var_ind)) deallocate(a%var_ind)
call deallocate_area(ia)
call deallocate_bound_collocate_points(ia)
call deallocate_cash(a%a%self_cash)
end

subroutine nullbound(ibnd)
!�������� ������ ��� �������
use pgmod
integer(4) ibnd !����� �������
type(TBound), pointer :: b
b=>gsarea%bnd(ibnd)
b%nline=0
b%npanel=0
b%self_cash%bnd_inited=.false.
b%self_cash%area_inited=.false.
b%cash=>b%self_cash
b%skip_getfun=.false.
end

subroutine pg_allocate_boundlines(nline)
!dec$ attributes dllexport:: pg_allocate_boundlines
!�������� ������ ��� �������
!��� boundLineType=1
use pgmod
integer(4) nline !����� �������� �������
integer(4) i
gsbnd%boundLineType=1
call deallocate_boundlines(gsarea%i,gsbnd%i)
gsbnd%nline=nline
allocate(gsbnd%line(nline))
allocate(gsbnd%geom_detale0(0))
gsbnd%ngeom_detale0=0
do i=1,nline
  call nullboundline(i)
enddo
end

subroutine allocate_geom_detale_for_next 
use pgmod
type(TBoundline_geomdetale), allocatable :: geom_detale1(:)
integer(4) n,n1,i
n=size(gsbnd%geom_detale0)
if (n==gsbnd%ngeom_detale0) then
  allocate(geom_detale1(n))
  geom_detale1=gsbnd%geom_detale0
  deallocate(gsbnd%geom_detale0)
  n1=n+gsbnd%nline
  allocate(gsbnd%geom_detale0(n1))
  gsbnd%geom_detale0(1:n)=geom_detale1
  deallocate(geom_detale1)
  do i=n+1,n1
    call null_geom_detale(gsbnd%i,i)
  enddo
endif
gsbnd%ngeom_detale0=gsbnd%ngeom_detale0+1
end

subroutine null_geom_detale(ib,igd)
use pgmod
integer(4) igd,ib
type(TBoundline_geomdetale), pointer :: gd
gd=>gsarea%bnd(ib)%geom_detale0(igd)
gd%mode=0
gd%s1=d0
gd%s2=d0
gd%x=d0
gd%y=d0
gd%x2=d0
gd%y2=d0
gd%r=d0
gd%gam1=d0
gd%gam2=d0
gd%i_begin=0
gd%i_end=0
gd%ibndl=0
gd%panelLmax=d0
end

subroutine pg_allocate_boundlines2(nline)
!dec$ attributes dllexport:: pg_allocate_boundlines2
!�������� ������ ��� �������
!��� boundLineType=2
use pgmod
integer(4) nline !����� �������� �������
integer(4) i
gsbnd%boundLineType=2
call deallocate_boundlines(gsarea%i,gsbnd%i)
gsbnd%nline=nline
allocate(gsbnd%line2(nline))
do i=1,nline
  call nullboundline2(i)
enddo
end
 
subroutine deallocate_boundlines(ia,ibnd)
!���������� ������ ��� �������
use pgmod
integer(4) ia   !����� �������
integer(4) ibnd !����� �������
integer(4) i
type(TBound), pointer :: b
type(TArea), pointer :: a
a=>gs%a(ia)
b=>a%bnd(ibnd)
if (allocated(b%line)) then
    do i=1,b%nline
        call deallocate_boundline(ia,ibnd,i)
    enddo
    deallocate(b%line)    
endif
if (allocated(b%geom_detale0)) deallocate(b%geom_detale0)
if (allocated(b%geom_detale)) deallocate(b%geom_detale)
if (allocated(b%line2)) then
    do i=1,b%nline
        call deallocate_boundline2(ia,ibnd,i)
    enddo
    deallocate(b%line2)    
endif
if (allocated(b%x)) deallocate(b%x)
if (allocated(b%y)) deallocate(b%y)
if (allocated(b%xc)) deallocate(b%xc)
if (allocated(b%yc)) deallocate(b%yc)
if (allocated(b%s)) deallocate(b%s)
if (allocated(b%sc)) deallocate(b%sc)
if (allocated(b%l)) deallocate(b%l)
if (allocated(b%z)) deallocate(b%z)
if (allocated(b%zc)) deallocate(b%zc)
if (allocated(b%ett)) deallocate(b%ett)
if (allocated(b%psiom)) deallocate(b%psiom)
if (allocated(b%psiind)) deallocate(b%psiind)
if (allocated(b%bb_psioml)) deallocate(b%bb_psioml)
if (allocated(b%cc_psioml)) deallocate(b%cc_psioml)
call deallocate_cash(b%self_cash)
end

subroutine allocate_bound_geomgu(b,inpanel)
!�������� ������ ��� ������� (���������)
!��� BoundLineType=1
use pgmod
type(TBound) b
integer(4) inpanel !���������� ������� � ������� �������
integer(4) npanel1
npanel1=inpanel+1
b%npanel=inpanel
allocate(b%x(npanel1))
allocate(b%y(npanel1))
allocate(b%xc(inpanel))
allocate(b%yc(inpanel))
allocate(b%s(npanel1))
allocate(b%sc(inpanel))
allocate(b%l(inpanel))
allocate(b%z(npanel1))
allocate(b%zc(inpanel))
allocate(b%ett(inpanel))
b%x=d0
b%y=d0
b%xc=d0
b%yc=d0
b%s=d0
b%sc=d0
b%l=d0
b%z=c0
b%zc=c0
b%ett=c0
allocate(b%psiom(b%npanel,gsarea%umax)) 
allocate(b%psiind(b%npanel,gsarea%umax)) 
b%psiom=d0
b%psiind=d0
end

subroutine allocate_bound_collocate_points(ia,incp)
!�������� ������ ��� ������� (����� ����������)
use pgmod
integer(4) incp !����� ����� ����������
integer(4) ia
type(TBound_Collocate_points), pointer :: cp
cp=>gs%a(ia)%cpp
cp%ncp=incp
allocate(cp%i(incp)) 
allocate(cp%ibnd(incp)) 
allocate(cp%iu(incp)) 
allocate(cp%ibndl(incp)) 
allocate(cp%j(incp)) 
cp%i=0
cp%ibnd=0
cp%iu=0
cp%ibndl=0
cp%j=0
end

subroutine allocate_area_collocate_points(ia,incp)
!�������� ������ ��� ������� (����� ����������)
use pgmod
integer(4) incp !����� ����� ����������
integer(4) ia
type(TArea_Collocate_points), pointer :: cp
cp=>gs%a(ia)%a%cppa
cp%ncp=incp
allocate(cp%itr(incp)) 
allocate(cp%iu(incp)) 
allocate(cp%j(incp)) 
cp%itr=0
cp%iu=0
cp%j=0
end

subroutine nullboundline(ibndl)
!�������� ������ ��� ������� �������
use pgmod
integer(4) ibndl !����� ������� �������
type(TBoundline), pointer :: b
b=>gsbnd%line(ibndl)
b%npanel=0
b%gu_mode=0
b%igd_begin=0
b%igd_end=0
b%gu_inited=.false.
b%ci%is_internal=.false.
b%ci%ia=0
b%ci%ibnd=0
b%ci%ibndl=0
end

subroutine pg_set_boundline_gu_mode(ibndl,imode)
!dec$ attributes dllexport:: pg_set_boundline_gu_mode
!������ ���������������� ��� ���������� �������
use pgmod
integer(4) ibndl !����� ������� �������
integer(4) imode !���������������� ��� ���������� �������
type(TBoundline), pointer :: b
b=>gsbnd%line(ibndl)
b%gu_mode=imode
end

subroutine nullboundline2(ibndl)
!�������� ������ ��� ������� �������
use pgmod
integer(4) ibndl !����� ������� �������
type(TBoundline2), pointer :: b
b=>gsbnd%line2(ibndl)
!!!���� ������
end

subroutine allocate_boundline_geom(ibndl,inpanel)
!�������� ������ ��� ������� ������� (���������)
!��� BoundLineType=1
use pgmod
integer(4) ibndl !����� ������� �������
integer(4) inpanel !���������� ������� � ������� �������
integer(4) npanel1
type(TBoundline), pointer :: b
b=>gsbnd%line(ibndl)
npanel1=inpanel+1
b%npanel=inpanel
allocate(b%x(npanel1))
allocate(b%y(npanel1))
b%x=d0
b%y=d0
end

subroutine allocate_boundline_geom2(ibndl)
!�������� ������ ��� ������� ������� (���������)
!��� BoundLineType=2
use pgmod
integer(4) ibndl !����� ������� �������
type(TBoundline2), pointer :: b
b=>gsbnd%line2(ibndl)
allocate(b%gc(gsarea%nu))
end

subroutine pg_allocate_bound_gu
!dec$ attributes dllexport:: pg_allocate_bound_gu
!�������� ������ ��� ������� (��������� �������)
!��� BoundLineType=1
use pgmod
integer(4) i,j
type(TBoundline), pointer :: bl
do j=1,gsbnd%nline
  bl=>gsbnd%line(j)
  allocate(bl%gu(gsarea%nu))
  allocate(bl%cp_line(gsarea%nu))
  do i=1,gsarea%nu
    call nullgu(j,i,.false.)
  enddo
  bl%cp_line=.true.
enddo
end

subroutine pg_allocate_boundline_gu_add
!dec$ attributes dllexport:: pg_allocate_boundline_gu_add
!�������� ������ ��� ������� ������� (�������������� ��������� �������)
!��� BoundLineType=1
use pgmod
call pg_allocate_boundline_gu_add2(.false.)
end

subroutine pg_allocate_boundline_gu_add2(cp_line)
!dec$ attributes dllexport:: pg_allocate_boundline_gu_add2
!�������� ������ ��� ������� ������� (�������������� ��������� �������)
!��� BoundLineType=1
use pgmod
logical cp_line !true - ���� �� ����� ��������� ����� ���������� ��� ����������� �������
integer(4) i
allocate(gsbndl%gu_add(gsarea%nu))
gsbndl%cp_line=cp_line
do i=1,gsarea%nu
  call nullgu(0,i,.true.)
enddo
end

subroutine pg_allocate_bound_gu2
!dec$ attributes dllexport:: pg_allocate_bound_gu2
!�������� ������ ��� ������� (��������� �������)
!��� BoundLineType=2
use pgmod
integer(4) j
type(TBoundline2), pointer :: b
do j=1,gsbnd%nline
  b=>gsbnd%line2(j)
  allocate(b%ga(gsarea%umax)) 
enddo
end

subroutine deallocate_boundline(ia,ibnd,ibndl)
!���������� ������ ��� �������
use pgmod
integer(4) ia,ibnd,ibndl
type(TBoundline), pointer :: b
type(TArea), pointer :: a
a=>gs%a(ia)
b=>a%bnd(ibnd)%line(ibndl)
if (allocated(b%x)) deallocate(b%x)
if (allocated(b%y)) deallocate(b%y)
call deallocate_boundline_gu(ia,ibnd,ibndl)
end

subroutine deallocate_domain_gu(ia)
use pgmod
integer(4) ia,i
type(TArea), pointer :: a
a=>gs%a(ia)
do i=1,a%nb
  call deallocate_bound_gu(ia,i)
enddo
end

subroutine deallocate_bound_gu(ia,ibnd)
use pgmod
integer(4) ia,ibnd,i
type(TBound), pointer :: b
type(TArea), pointer :: a
a=>gs%a(ia)
b=>a%bnd(ibnd)
do i=1,b%nline
  select case (b%boundLineType)
  case (1)
    call deallocate_boundline_gu(ia,ibnd,i)
  case (2)
    call deallocate_boundline2_ga(ia,ibnd,i)
  endselect
enddo
end

subroutine deallocate_boundline_gu(ia,ibnd,ibndl)
use pgmod
integer(4) ia,ibnd,ibndl,i
!type(TBound), pointer :: b0
type(TBoundline), pointer :: b
type(TArea), pointer :: a
a=>gs%a(ia)
!b0=>a%bnd(ibnd)
!if (b0%boundLineType/=1) return
!b=>b0%line(ibndl)
b=>a%bnd(ibnd)%line(ibndl)
if (allocated(b%gu)) then
  do i=1,a%nu
    call deallocate_gu(ia,ibnd,ibndl,i,.false.)
  enddo
  deallocate(b%gu)
  if (allocated(b%use_gu)) deallocate(b%use_gu)
endif
if (allocated(b%gu_add)) then
  do i=1,a%nu
    call deallocate_gu(ia,ibnd,ibndl,i,.true.)
  enddo
  deallocate(b%gu_add)
endif
if (allocated(b%cp_line)) deallocate(b%cp_line)
end

subroutine deallocate_boundline2(ia,ibnd,ibndl)
!���������� ������ ��� �������
use pgmod
integer(4) ia,ibnd,ibndl,i
type(TBoundline2), pointer :: b
type(TArea), pointer :: a
a=>gs%a(ia)
b=>a%bnd(ibnd)%line2(ibndl)
call deallocate_boundline2_ga(ia,ibnd,ibndl)
if (allocated(b%gc)) then
  do i=1,a%nu
    call deallocate_gc(ia,ibnd,ibndl,i)
  enddo
  deallocate(b%gc)
endif
end

subroutine deallocate_boundline2_ga(ia,ibnd,ibndl)
!���������� ������ ��� �������
use pgmod
integer(4) ia,ibnd,ibndl,i
!type(TBound), pointer :: b0
type(TBoundline2), pointer :: b
type(TArea), pointer :: a
a=>gs%a(ia)
!b0=>a%bnd(ibnd)
!if (b0%boundLineType/=2) return
!b=>b0%line2(ibndl)
b=>a%bnd(ibnd)%line2(ibndl)
if (allocated(b%ga)) then
  do i=1,a%umax
    call deallocate_ga(ia,ibnd,ibndl,i)
  enddo
  deallocate(b%ga)
endif
end

subroutine deallocate_bound_collocate_points(ia)
!���������� ������ ��� ����� ����������
use pgmod
integer(4) ia
type(TBound_Collocate_points), pointer :: cp
cp=>gs%a(ia)%cpp
if (allocated(cp%i)) deallocate(cp%i)
if (allocated(cp%ibndl)) deallocate(cp%ibndl)
if (allocated(cp%iu)) deallocate(cp%iu)
if (allocated(cp%ibnd)) deallocate(cp%ibnd)
if (allocated(cp%j)) deallocate(cp%j)
cp%inited=.false.
end

subroutine deallocate_area_collocate_points(ia)
!���������� ������ ��� ����� ����������
use pgmod
integer(4) ia
type(TArea_Collocate_points), pointer :: cp
cp=>gs%a(ia)%a%cppa
if (allocated(cp%itr)) deallocate(cp%itr)
if (allocated(cp%iu)) deallocate(cp%iu)
if (allocated(cp%j)) deallocate(cp%j)
cp%inited=.false.
end

subroutine nullgu(ibndl,igu,is_add)
!�������� ������, ��������� � ��������
use pgmod
integer(4) igu !����� ���������� �������
integer(4) ibndl
type(TBoundline_GU), pointer :: gu
logical is_add
if (is_add) then
  gu=>gsbndl%gu_add(igu)
  gu%bndl=>gsbndl
else
  gu=>gsbnd%line(ibndl)%gu(igu)
  gu%bndl=>gsbnd%line(ibndl)
endif
gu%bndu=0
gu%constvalind=0
gu%constval=d0
end

subroutine allocate_gu(igu,ival)
!�������� ������ ��� ���������� �������
use pgmod
integer(4) igu !����� ���������� �������
integer(4) ival !����������� ������� �������� (���� = 1)
type(TBoundline_GU), pointer :: gu
if (igu>gsarea%nu) then
  gu=>gsbndl%gu_add(igu-gsarea%nu)
else
  gu=>gsbndl%gu(igu)
endif
if (.not.allocated(gu%bndval)) allocate(gu%bndval(gsbndl%npanel,ival))
gu%bndval=d0
end

subroutine allocate_gu_gen(gu,n)
!�������� ������ ��� ���������� �������
use pgmod
integer(4) n !����������� ��������
type(TBoundline_GU) gu !��������� �������
integer(4) i
call deallocate_gu_gen(gu%bndg)
gu%bndg%is_direct=.true.
gu%bndg%ia=0
gu%bndg%ibnd=0
gu%bndg%ibndl=0
gu%bndg%n=n
allocate(gu%bndg%carr(0:n))
do i=0,n
  gu%bndg%carr(i)%c=d0
enddo
allocate(gu%bndg%second_bnd(n))
allocate(gu%bndg%bndf(n))
gu%bndg%second_bnd=.false.
gu%bndg%bndf=0
end

subroutine allocate_gu_genGlobal(gu,n)
!�������� ������ ��� ���������� �������
use pgmod
integer(4) n !����������� ��������
type(TBoundline_GU) gu !��������� �������
integer(4) i
call deallocate_gu_genGlobal(gu%guglob)
gu%guglob%n=n
allocate(gu%guglob%carr(0:n))
do i=0,n
  gu%guglob%carr(i)%c=d0
enddo
allocate(gu%guglob%gu_indsi(n))
gu%guglob%gu_indsi=0
end

subroutine deallocate_gu_gen(gug)
!���������� ������ ��� ���������� �������
use pgmod
type(TBoungline_GU_gen) gug
integer(4) i
if (allocated(gug%second_bnd)) deallocate(gug%second_bnd)
if (allocated(gug%bndf)) deallocate(gug%bndf)
if (allocated(gug%carr)) then
  do i=0,gug%n
    if (allocated(gug%carr(i)%v)) deallocate(gug%carr(i)%v)
  enddo
  deallocate(gug%carr)
endif
end

subroutine deallocate_gu_genGlobal(gug)
!���������� ������ ��� ���������� �������
use pgmod
type(TBoungline_GU_genGlobal) gug
integer(4) i
if (allocated(gug%carr)) then
  do i=0,gug%n
    if (allocated(gug%carr(i)%v)) deallocate(gug%carr(i)%v)
  enddo
  deallocate(gug%carr)
endif
if (allocated(gug%gu_indsi)) deallocate(gug%gu_indsi)
end

subroutine allocate_gu_use
!�������� ������ ��� use_gu �� ������� ������� �������
use pgmod
allocate(gsbndl%use_gu(gsbndl%npanel,gsarea%nu))
gsbndl%use_gu=.true.
end

subroutine pg_allocate_and_set_ga2(iga,n)
!dec$ attributes dllexport:: pg_allocate_and_set_ga2
!�������� ������ ��� ������������� ����������������� �������
!���������� ������� �� ��������� (-1 - ����� ������������� � get_ind_area)
!��� BoundLineType=2
use pgmod
integer(4) iga !����� ���������������� �������
integer(4) n
type(TBounLineFuncApprox), pointer :: ga
ga=>gsbndl2%ga(iga)
call pg_allocate_and_set_ga2_base(ga,n)
end

subroutine pg_allocate_and_set_ga2_base(ga,n)
!dec$ attributes dllexport:: pg_allocate_and_set_ga2_base
!�������� ������ ��� ������������� ����������������� �������
!���������� ������� �� ��������� (-1 - ����� ������������� � get_ind_area)
!��� BoundLineType=2
use pgmod
integer(4) n
type(TBounLineFuncApprox) ga
allocate(ga%val(n))
allocate(ga%ind(n))
allocate(ga%kk(n))
allocate(ga%valkk(n))
ga%n=n
ga%kk=d1
ga%ga_ref=0
call set_ga2_base(ga)
end

subroutine set_ga2_base(ga)
!���������� ������� �� ��������� (-1 - ����� ������������� � get_ind_area)
!��� BoundLineType=2
use pgmod
type(TBounLineFuncApprox) ga
ga%val=d0
ga%valkk=d0
ga%ind=0
if (ga%typea==1) then
  if (ga%bnda==2) then
    ga%ind(1)=-1
  elseif (ga%bnda==3) then
    ga%ind(2:3)=-1
  elseif (ga%bnda==4) then
    ga%ind(1:3)=-1
  endif
endif
end

subroutine pg_allocate_gc2(ibndl,igc,n)
!dec$ attributes dllexport:: pg_allocate_gc2
!�������� ������ ��� ����� ����������
!��� BoundLineType=2
use pgmod
integer(4) ibndl
integer(4) igc !����� ���������������� �������
integer(4) n   !���������� ����� ����������
type(TBoundline_Collocate), pointer :: gc
gc=>gsbnd%line2(ibndl)%gc(igc)
allocate(gc%g(n))
allocate(gc%z(n))
allocate(gc%ztc(n))
gc%g=d0
gc%z=c0
gc%ztc=c0
gc%n=n
end

subroutine deallocate_gu(ia,ibnd,ibndl,igu,is_add)
!���������� ������ ��� ���������� �������
use pgmod
integer(4) ia !����� �������
integer(4) ibnd  !����� �������
integer(4) ibndl !����� ������� �������
integer(4) igu !����� ���������� �������
type(TBoundline_GU), pointer :: gu
logical is_add
if (is_add) then
  gu=>gs%a(ia)%bnd(ibnd)%line(ibndl)%gu_add(igu)
else
  gu=>gs%a(ia)%bnd(ibnd)%line(ibndl)%gu(igu)
endif
if (allocated(gu%bndval)) deallocate(gu%bndval)
select case (gu%bndu)
case (3)
  call deallocate_gu_gen(gu%bndg)
case (5)
  call deallocate_gu_genGlobal(gu%guglob)
endselect
end

subroutine deallocate_gc(ia,ibnd,ibndl,igu)
!���������� ������ ��� ���������� �������
use pgmod
integer(4) ia !����� �������
integer(4) ibnd  !����� �������
integer(4) ibndl !����� ������� �������
integer(4) igu !����� ���������� �������
type(TBoundline_Collocate), pointer :: gc
gc=>gs%a(ia)%bnd(ibnd)%line2(ibndl)%gc(igu)
if (allocated(gc%g)) deallocate(gc%g)
if (allocated(gc%z)) deallocate(gc%z)
if (allocated(gc%ztc)) deallocate(gc%ztc)
end

subroutine deallocate_ga(ia,ibnd,ibndl,iga)
!���������� ������ ��� ���������� �������
use pgmod
integer(4) ia !����� �������
integer(4) ibnd  !����� �������
integer(4) ibndl !����� ������� �������
integer(4) iga !����� ���������� �������
type(TBounLineFuncApprox), pointer :: ga
ga=>gs%a(ia)%bnd(ibnd)%line2(ibndl)%ga(iga)
call pg_deallocate_ga_base(ga)
end

subroutine pg_deallocate_ga_base(ga)
!dec$ attributes dllexport:: pg_deallocate_ga_base
!���������� ������ ��� ���������� �������
use pgmod
type(TBounLineFuncApprox) ga
if (allocated(ga%val)) deallocate(ga%val)
if (allocated(ga%ind)) deallocate(ga%ind)
if (allocated(ga%kk)) deallocate(ga%kk)
if (allocated(ga%valkk)) deallocate(ga%valkk)
end

subroutine pg_allocate_area(inpart)
!dec$ attributes dllexport:: pg_allocate_area
!���������������� ���������� �������� � ������� ���������
use pgmod
integer(4) inpart   !���������� ��������
integer(4) i
type(areatype), pointer :: a
a=>gsarea%a
a%npart=inpart
allocate(a%part(inpart))
do i=1,inpart
  call null_areapart(i)
enddo
end

subroutine allocate_area_geom(in,intr,inpe)
!�������� ������ ��� ��������� ����� � �������
use pgmod
integer(4) in   !���������� �����
integer(4) intr !���������� �������������
integer(4) inpe !���������� ������ � ����� �������� (3 - �����������, 4 - ���������������)
type(areatype), pointer :: a
a=>gsarea%a
a%n=in
a%ntr=intr
a%npe=inpe
allocate(a%zm(in))
allocate(a%trm(inpe,intr))
allocate(a%zmc(intr))
allocate(a%npe_ar(intr))
a%zm=c0
a%trm=0
a%zmc=c0
a%npe_ar=0
end

subroutine null_areapart(ia)
!�������� ������, ��������� � ��������
use pgmod
integer(4) ia !����� ������� �������
type(areapart), pointer :: a
a=>gsarea%a%part(ia)
a%n=0
a%ntr=0
a%npe=0
end

subroutine pg_allocate_areapart_geom(in,intr,inpe)
!dec$ attributes dllexport:: pg_allocate_areapart_geom
!�������� ������ ��� ��������� ����� � ������� �������
use pgmod
integer(4) in   !���������� �����
integer(4) intr !���������� �������������
integer(4) inpe !���������� ������ � ����� �������� (3 - �����������, 4 - ���������������)
call pg_allocate_areapart_geom_n(in)
call pg_allocate_areapart_geom_tr(intr,inpe)
end

subroutine pg_allocate_areapart_geom_n(in)
!dec$ attributes dllexport:: pg_allocate_areapart_geom_n
!�������� ������ ��� ��������� ����� � ������� �������
!������ ����
use pgmod
integer(4) in   !���������� �����
type(areapart), pointer :: a
a=>gsareapart
a%n=in
allocate(a%zm(in))
a%zm=c0
end

subroutine pg_allocate_areapart_geom_tr(intr,inpe)
!dec$ attributes dllexport:: pg_allocate_areapart_geom_tr
!�������� ������ ��� ��������� ����� � ������� �������
!������ ��������
use pgmod
integer(4) intr !���������� �������������
integer(4) inpe !���������� ������ � ����� �������� (3 - �����������, 4 - ���������������)
type(areapart), pointer :: a
a=>gsareapart
a%ntr=intr
a%npe=inpe
allocate(a%trm(inpe,intr))
a%trm=0
end

subroutine pg_init_areapart_geom_xy(x,y)
!dec$ attributes dllexport:: pg_init_areapart_geom_xy
!���������������� ��������� ����� � ������� ������� (���������� �����)
use pgmod
real(8) x(gsareapart%n),y(gsareapart%n) !����������
gsareapart%zm=dcmplx(x,y)
end

subroutine pg_init_areapart_geom_tr(tr)
!dec$ attributes dllexport:: pg_init_areapart_geom_tr
!���������������� ��������� ����� � ������� ������� (������������ ��� ����������������)
use pgmod
integer(4) tr(gsareapart%npe*gsareapart%ntr)
gsareapart%trm=reshape(tr, (/gsareapart%npe,gsareapart%ntr/))
end

subroutine allocate_areaval(i,mode)
use pgmod
integer(4) i,mode
type(TAreaValue), pointer :: av
if (mode==0) then
  av=>gsarea%a%areaval(i)
else
  av=>gsarea%a%arval_funp(i)
endif
allocate(av%v(gsarea%a%ntr))
av%v=d0
end

subroutine pg_allocate_area_gu
!dec$ attributes dllexport:: pg_allocate_area_gu
!�������� ������ ��� ������� � ������� (������ �� �� �������) � �������
use pgmod
integer(4) ku   !��� ��������� ��������� � �������
integer(4) ku2   !��� ������� ��������� � �������
integer(4) n,k,i
type(areatype), pointer :: a
logical allocated1
ku=gsarea%type_eq(1)
ku2=0
if (gsarea%nu==2) ku2=gsarea%type_eq(2)
a=>gsarea%a
n=2
if((ku==8).or.(ku==11)) n=3
if((ku==16).or.(ku==17)) n=5
if((ku==18).or.(ku==19).or.(ku==20)) n=4
if (ku==22.or.ku==24.or.ku==25) n=7
if (ku==26) n=4
allocate(a%areaval(n))
!������������ ���������
allocated1=.false.
select case (ku)
case (4,6,17,18,21)
  call allocate_areaval(1,0)
  allocated1=.true.
end select
if (.not.allocated1) then
  select case (ku2)
  case (7)
    call allocate_areaval(1,0)
  end select
endif
!��������� � ������������ � �������
select case (ku)
case (5,6)
  call allocate_areaval(2,0)
case (8)
  call allocate_areaval(3,0)
  if (ku2==14) call allocate_areaval(2,0)
case (11)
  call allocate_areaval(3,0)
case (16,17)
  call allocate_areaval(4,0)
  call allocate_areaval(5,0)
case (18:20)
  call allocate_areaval(4,0)
case (22,24:26)
  do i=1,gsarea%n_eq_var
    if (gsarea%eq_var(i)) call allocate_areaval(i,0)
  enddo
endselect
allocate(a%arval_funp(gsarea%nu))
select case (ku)
case (4:6,8,11,16:22,24:26)
  call allocate_areaval(1,1)
end select
select case (ku2)
case (7,14,23)
  call allocate_areaval(2,1)
end select
call init_mesh_gu2
if (gsarea%haveAreaEq) then
  allocate(a%areaind(a%ntr,gsarea%umaxtr))
  allocate(a%psiarea(a%ntr,gsarea%umaxtr))
  a%areaind=d0
  a%psiarea=d0
endif
select case (ku)
case (22,24:26)
  allocate(a%eq_ind(gsarea%umaxtr))
  k=0
  do i=2,gsarea%n_eq_var
    if (gsarea%eq_var(i)) then
      k=k+1
      a%eq_ind(k)=i-1
    endif
  enddo
endselect
end

subroutine deallocate_area(ia)
!����������� ������ � �������
use pgmod
integer(4) ia   !����� �������
type(areatype), pointer :: a
integer(4) i
a=>gs%a(ia)%a
if (allocated(a%part)) then
  do i=1,a%npart
    call deallocate_areapart(ia,i)
  enddo
  deallocate(a%part)
endif
if (allocated(a%zm)) deallocate(a%zm)
if (allocated(a%trm)) deallocate(a%trm)
if (allocated(a%zmc)) deallocate(a%zmc)
if (allocated(a%npe_ar)) deallocate(a%npe_ar)
if (allocated(a%areaval)) then
  do i=1,ubound(a%areaval,1)
    if (allocated(a%areaval(i)%v)) deallocate(a%areaval(i)%v)
  enddo
  deallocate(a%areaval)
endif
if (allocated(a%arval_funp)) then
  do i=1,ubound(a%arval_funp,1)
    if (allocated(a%arval_funp(i)%v)) deallocate(a%arval_funp(i)%v)
  enddo
  deallocate(a%arval_funp)
endif
if (allocated(a%areaind)) deallocate(a%areaind)
if (allocated(a%psiarea)) deallocate(a%psiarea)
if (allocated(a%eq_ind)) deallocate(a%eq_ind)
call deallocate_area_collocate_points(ia)
end

subroutine deallocate_areapart(ia,ipart)
!����������� ������ � �������
use pgmod
integer(4) ia   !����� �������
integer(4) ipart   !����� ������� �������
type(areapart), pointer :: a
a=>gs%a(ia)%a%part(ipart)
if (allocated(a%zm)) deallocate(a%zm)
if (allocated(a%trm)) deallocate(a%trm)
end

subroutine init_mesh_gu
!������ ��� ��������� � ������� �������
use pgmod
integer(4) type_eq
type_eq=gsarea%type_eq(1)
selectcase (type_eq)
case (8)
  gsArea%a%areaEq=gsarea%type_eq(2)
  gsarea%haveAreaEq=.true.
case (5,6,16:20,22,24:26)
  gsArea%a%areaEq=type_eq
  gsarea%haveAreaEq=.true.
endselect
end

subroutine init_mesh_gu2
!������ ��� ��������� � ������� �������
use pgmod
integer(4) type_eq,i
type_eq=gsarea%type_eq(1)
if (gsarea%haveAreaEq) then
  gsarea%umaxtr=1
  selectcase (type_eq)
  case (16,17)
    gsarea%umaxtr=2
  case (22,24:26)
    gsarea%umaxtr=0
    do i=2,gsarea%n_eq_var
      if (gsarea%eq_var(i)) gsarea%umaxtr=gsarea%umaxtr+1
    enddo
  endselect
endif
end

subroutine pg_init_area_gu_const(val,k)
!dec$ attributes dllexport:: pg_init_area_gu_const
!�������� � ������������� ��� ������ ��������, ����������� �.�.�
!���������� �������� �� ���� �������������
use pgmod
real(8) val !�������� � ������������� 
integer(4) k !��� ��������, ��. �������� areatype%areaval
real(8) g(gsarea%a%ntr) !������ ��������
g=val
call pg_init_area_gu(g,k)
end

subroutine pg_init_areapart_gu_const(val,k)
!dec$ attributes dllexport:: pg_init_areapart_gu_const
!�������� � ������������� ��� ������ ��������, ����������� �.�.�
!���������� �������� �� ���� �������������
!��� �������
use pgmod
real(8) val !�������� � ������������� 
integer(4) k !��� ��������, ��. �������� areatype%areaval
real(8) g(gsareapart%ntr) !������ ��������
g=val
call pg_init_areapart_gu(g,k)
end

subroutine pg_init_area_gu_oss
!dec$ attributes dllexport:: pg_init_area_gu_oss
!�������� � ������������� ��� ��������������� ����� 19,20
use pgmod
integer(4) i
integer(4) k !��� ��������, ��. �������� areatype%areaval
real(8) g(gsarea%a%ntr) !������ ��������
k=4
do i=1,gsarea%a%ntr
  g(i)=-gsarea%const%k_oss
  if (gsarea%type_eq(1)==19) g(i)=g(i)/dimag(gsarea%a%zmc(i))
enddo
!g=dimag(gsarea%a%zmc)
call pg_init_area_gu(g,k)
end

subroutine pg_init_area_gu(g,k)
!dec$ attributes dllexport:: pg_init_area_gu
!�������� � ������������� ��� ������ ��������, ����������� �.�.�
use pgmod
integer(4) k !��� ��������, ��. �������� areatype%areaval
real(8) g(gsarea%a%ntr) !������ ��������
integer(4) i
type(areapart), pointer :: a
do i=1,gsarea%a%npart
  a=>gsarea%a%part(i)
  call init_areapart_gu(i,g(a%ntr_begin:a%ntr_end),k)
enddo
end

subroutine pg_init_areapart_gu(g,k)
!dec$ attributes dllexport:: pg_init_areapart_gu
!�������� � ������������� ��� ������ ��������, ����������� �.�.�
!��� �������
use pgmod
integer(4) k !��� ��������, ��. �������� areatype%areaval
real(8) g(gsareapart%ntr) !������ ��������
call init_areapart_gu(gsareapart%i,g,k)
end

subroutine pg_init_area_gu_test(g,k,err,ierr,erre)
!dec$ attributes dllexport:: pg_init_area_gu_test
!�������� � ������������� ��� ������ ��������, ����������� �.�.�
use pgmod
integer(4) k !��� ��������, ��. �������� areatype%areaval
real(8) g(gsarea%a%ntr) !������ ��������
real(8) err !���������� �����������
real(8) erre !������������� �����������
integer(4) ierr !������ ������, �� ������� ���������� ���������� �����������
integer(4) i
type(areapart), pointer :: a
err=d0
erre=d0
ierr=-1
do i=1,gsarea%a%npart
  a=>gsarea%a%part(i)
  call init_areapart_gu_test(i,g(a%ntr_begin:a%ntr_end),k,.true.,err,ierr,erre)
enddo
erre=erre/gsarea%a%ntr
end

subroutine pg_init_areapart_gu_test(g,k,err,ierr,erre)
!dec$ attributes dllexport:: pg_init_areapart_gu_test
!�������� � ������������� ��� ������ ��������, ����������� �.�.�
!��� �������
use pgmod
integer(4) k !��� ��������, ��. �������� areatype%areaval
real(8) g(gsareapart%ntr) !������ ��������
real(8) err !���������� �����������
real(8) erre !������������� �����������
integer(4) ierr !������ ������ (� 1 ��� �������), �� ������� ���������� ���������� �����������
err=d0
erre=d0
ierr=-1
call init_areapart_gu_test(gsareapart%i,g,k,.true.,err,ierr,erre)
erre=erre/gsareapart%ntr
ierr=ierr-gsareapart%ntr_begin+1
end

function test_area_gu(k)
!�������� � ������������� ��� ������ ��������, ����������� �.�.�
use pgmod
logical test_area_gu
integer(4) k !��� ��������, ��. �������� areatype%areaval
type(areatype), pointer :: a
type(TAreaValue), pointer :: av
a=>gsarea%a
test_area_gu=.false.
if (.not.allocated(a%areaval)) call gs_print_stop("Erorr test_area_gu")
if (k>ubound(a%areaval,1)) call gs_print_stop("Erorr test_area_gu")
av=>a%areaval(k)
if (.not.allocated(av%v)) call gs_print_stop("Erorr test_area_gu")
test_area_gu=.true.
end

subroutine init_areapart_gu(ipart,g,k)
use pgmod
integer(4) ipart !����� �������
integer(4) k !��� ��������, ��. �������� areatype%areaval
real(8) g(gsarea%a%part(ipart)%ntr) !������ ��������
real(8) err,erre
integer(4) ierr
call init_areapart_gu_test(ipart,g,k,.false.,err,ierr,erre)
end

subroutine init_areapart_gu_test(ipart,g,k,need_test,err,ierr,erre)
use pgmod
integer(4) ipart !����� �������
integer(4) k !��� ��������, ��. �������� areatype%areaval
real(8) g(gsarea%a%part(ipart)%ntr) !������ ��������
logical need_test
real(8) err,terr,erre
integer(4) ierr,i
logical test_area_gu
type(areatype), pointer:: aa
type(areapart), pointer:: a
type(TAreaValue), pointer :: av
if (.not.test_area_gu(k)) return
aa=>gsarea%a
a=>aa%part(ipart)
av=>aa%areaval(k)
if (need_test) then
  do i=a%ntr_begin,a%ntr_end
    terr=dabs(av%v(i)-g(i))
    if (ierr<0.or.terr>err) then
      err=terr
      ierr=i
    endif
    erre=erre+terr
  enddo
endif
av%v(a%ntr_begin:a%ntr_end)=g
end

subroutine init_gsMain_s_ndns_one(k)
use pgmod
integer(4) i,j,k
real(8) shift,shift2
type(Tsndns_type), pointer :: sn
sn=>gsMain%sn(k)
shift2=d5/sn%n
if (.not.allocated(sn%s_ndns)) allocate(sn%s_ndns(2,sn%n))
do j=1,2
  shift=d0
  if (j==2.and.sn%use_shift2) shift=shift2
  do i=1,sn%n
    sn%s_ndns(j,i)=(i-d1)/3.0d0+shift+sn%shift
  enddo
enddo
end

subroutine init_gsMain_s_ndns(n)
use pgmod
integer(4) n !����� ����� ����������
type(Tsndns_type), pointer :: sn
if (.not.allocated(gsMain%sn)) then
  allocate(gsMain%sn(1))
  gsMain%nsn=1
  sn=>gsMain%sn(1)
  sn%n=n
  sn%shift=d0
  sn%use_shift2=.false.
  call init_gsMain_s_ndns_one(1)
endif
end

subroutine init_gsMain_s_ndns2(n)
use pgmod
integer(4) n !����� ����� ����������
integer(4) i
type(Tsndns_type), pointer :: sn
if (.not.allocated(gsMain%sn)) then
  allocate(gsMain%sn(2))
  gsMain%nsn=2
  do i=1,2
    sn=>gsMain%sn(i)
    sn%n=n
    sn%shift=d0
    if (i==2) sn%shift=d5/n
    sn%use_shift2=.false.
    call init_gsMain_s_ndns_one(i)
  enddo
endif
end

subroutine initCashIntegral(c)
use pgmod
type(TCashIntegral) c
allocate(c%inited(0:max_int))
allocate(c%solved(0:max_int))
allocate(c%intvals(0:max_int))
c%inited=.false.
c%solved=.false.
end

subroutine init_bb_cash(ia,knd)
use pgmod
integer(4) ia,knd
type(TBound), pointer :: b
type(TArea), pointer :: a
a=>gs%a(ia)
b=>a%bnd(knd)
b%self_cash%ia=ia
call init_bound_cash(b%cash,a%nb)
end

subroutine init_ab_cash(ia)
use pgmod
integer(4) ia
type(TArea), pointer :: a
a=>gs%a(ia)
a%a%self_cash%ia=ia
call init_bound_cash(a%a%cash,a%nb)
end

subroutine init_ba_cash(ia,knd)
use pgmod
integer(4) ia,knd
type(TBound), pointer :: b
b=>gs%a(ia)%bnd(knd)
b%self_cash%ia=ia
call init_area_cash(b%cash)
end

subroutine init_aa_cash(ia)
use pgmod
integer(4) ia
type(TArea), pointer :: a
a=>gs%a(ia)
a%a%self_cash%ia=ia
call init_area_cash(a%a%cash)
end

subroutine init_bound_cash(cash,nb)
use pgmod
integer(4) i,nb
type(TCashIntegral), pointer :: c
type(TCash), target :: cash
if (.not.cash%bnd_inited) then
  allocate(cash%bnd(nb))
  cash%bnd_inited=.true.
  do i=1,nb
    c=>cash%bnd(i)
    call initCashIntegral(c)
  enddo
endif
end

subroutine init_area_cash(cash)
use pgmod
type(TCash) cash
if (.not.cash%area_inited) then
  cash%area_inited=.true.    
  call initCashIntegral(cash%area)
endif
end

subroutine initCashIntegral_nf(c,nf,n1,n2,ia)
use pgmod
integer(4) nf,n1,n2,ia
type(TCashIntegral) c
real(8) cash_size_mb
allocate(c%intvals(nf)%i(n1,n2))
c%inited(nf)=.true.
c%intvals(nf)%i=real8_inf
cash_size_mb=8.0d0*n1*n2/1024**2
gs%a(ia)%cash_size_mb=gs%a(ia)%cash_size_mb+cash_size_mb
gs%cash_size_mb=gs%cash_size_mb+cash_size_mb
end

function get_bb_cash_integral(ia,knd,j,knd0,i,nf)
!��� �������������� �� ������ ����������� ����� �� ������
use pgmod
!j - ����� ��� ������, �� ������� ���� ��������������
!i - ����� ����������� ����� (������ � ������� ���������)
!nf - ����� ���������
!knd0 - ������ ������� � ������� � ����������� ������
!ia, knd - ������ ������� � �������, ��� ��������� ������ ��������������
integer(4) ia,knd,j,knd0,i,nf
real(8) get_bb_cash_integral,get_bb_cash_integral_
if (gs_use_cash_lock_get) then
  !$omp critical (lock_cash)
  get_bb_cash_integral=get_bb_cash_integral_(ia,knd,j,knd0,i,nf)
  !$omp end critical (lock_cash)
else
  get_bb_cash_integral=get_bb_cash_integral_(ia,knd,j,knd0,i,nf)
endif
end

function get_bb_cash_integral_(ia,knd,j,knd0,i,nf)
!��� �������������� �� ������ ����������� ����� �� ������
use pgmod
!j - ����� ��� ������, �� ������� ���� ��������������
!i - ����� ����������� ����� (������ � ������� ���������)
!nf - ����� ���������
!knd0 - ������ ������� � ������� � ����������� ������
!ia, knd - ������ ������� � �������, ��� ��������� ������ ��������������
integer(4) ia,knd,j,knd0,i,nf
real(8) get_bb_cash_integral_
type(TCashIntegral), pointer :: c
type(TBound), pointer :: b
type(TArea), pointer :: a
type(TCashValue), pointer :: cv
!call init_bb_cash(ia,knd)
a=>gs%a(ia)
b=>a%bnd(knd)
c=>b%cash%bnd(knd0)
if (.not.c%inited(nf)) then
  get_bb_cash_integral_=real8_inf
else
  cv=>c%intvals(nf)
  get_bb_cash_integral_=cv%i(j,i)
endif
end

function get_bb_cash_integral_2(a,knd,j,knd0,i,nf) result(res)
!��� �������������� �� ������ ����������� ����� �� ������
!���������� ������� ��� ��������
use pgmod
!j - ����� ��� ������, �� ������� ���� ��������������
!i - ����� ����������� ����� (������ � ������� ���������)
!nf - ����� ���������
!knd0 - ������ ������� � ������� � ����������� ������
!a   - ������� ��� ��������� ������ ��������������
!knd - ����� �������, ��� ��������� ������ ��������������
integer(4) knd,j,knd0,i,nf
real(8) res
!DEC$ IF DEFINED (DEBUG)
type(TCashIntegral), pointer :: c
type(TBound), pointer :: b
type(TArea), target :: a
type(TCashValue), pointer :: cv
b=>a%bnd(knd)
c=>b%cash%bnd(knd0)
cv=>c%intvals(nf)
res=cv%i(j,i)
!DEC$ ELSE
type(TArea) a
res=a%bnd(knd)%cash%bnd(knd0)%intvals(nf)%i(j,i)
!DEC$ ENDIF
end

subroutine set_bb_cash_integral(ia,knd,j,knd0,i,nf,val)
!��� �������������� �� ������ ����������� ����� �� ������
use pgmod
!j - ����� ��� ������, �� ������� ���� ��������������
!i - ����� ����������� ����� (������ � ������� ���������)
!nf - ����� ���������
!knd0 - ������ ������� � ������� � ����������� ������
!ia, knd - ������ ������� � �������, ��� ��������� ������ ��������������
integer(4) ia,knd,j,knd0,i,nf
real(8) val
!$omp critical (lock_cash)
call set_bb_cash_integral_(ia,knd,j,knd0,i,nf,val)
!$omp end critical (lock_cash)
end

subroutine set_bb_cash_integral_(ia,knd,j,knd0,i,nf,val)
!��� �������������� �� ������ ����������� ����� �� ������
use pgmod
!j - ����� ��� ������, �� ������� ���� ��������������
!i - ����� ����������� ����� (������ � ������� ���������)
!nf - ����� ���������
!knd0 - ������ ������� � ������� � ����������� ������
!ia, knd - ������ ������� � �������, ��� ��������� ������ ��������������
integer(4) ia,knd,j,knd0,i,nf
real(8) val
type(TCashIntegral), pointer :: c
type(TBound), pointer :: b,b0
type(TArea), pointer :: a
type(TCashValue), pointer :: cv
a=>gs%a(ia)
b=>a%bnd(knd)
b0=>a%bnd(knd0)
c=>b%cash%bnd(knd0)
if (.not.c%inited(nf)) call initCashIntegral_nf(c,nf,b%npanel,b0%npanel,b%cash%ia)
cv=>c%intvals(nf)
cv%i(j,i)=val
end

subroutine set_bb_cash_integral_2(ia,knd,j,knd0,i,nf,val)
!��� �������������� �� ������ ����������� ����� �� ������
!���������� ������� ��� ��������
use pgmod
!j - ����� ��� ������, �� ������� ���� ��������������
!i - ����� ����������� ����� (������ � ������� ���������)
!nf - ����� ���������
!knd0 - ������ ������� � ������� � ����������� ������
!ia, knd - ������ ������� � �������, ��� ��������� ������ ��������������
integer(4) ia,knd,j,knd0,i,nf
real(8) val
!DEC$ IF DEFINED (DEBUG)
type(TCashIntegral), pointer :: c
type(TBound), pointer :: b,b0
type(TArea), pointer :: a
type(TCashValue), pointer :: cv
a=>gs%a(ia)
b=>a%bnd(knd)
b0=>a%bnd(knd0)
c=>b%cash%bnd(knd0)
cv=>c%intvals(nf)
cv%i(j,i)=val
!DEC$ ELSE
gs%a(ia)%bnd(knd)%cash%bnd(knd0)%intvals(nf)%i(j,i)=val
!DEC$ ENDIF
end

function get_ba_cash_integral(ia,knd,j,i,nf)
!��� �������������� �� ������ ����������� ����� � ������������
use pgmod
!j - ����� ��� ������, �� ������� ���� ��������������
!i - ����� ����������� ����� (������������)
!nf - ����� ���������
!ia, knd - ������ ������� � �������, ��� ��������� ������ ��������������
integer(4) ia,knd,j,i,nf
real(8) get_ba_cash_integral,get_ba_cash_integral_
if (gs_use_cash_lock_get) then
  !$omp critical (lock_cash)
  get_ba_cash_integral=get_ba_cash_integral_(ia,knd,j,i,nf)
  !$omp end critical (lock_cash)
else
  get_ba_cash_integral=get_ba_cash_integral_(ia,knd,j,i,nf)
endif
end

function get_ba_cash_integral_(ia,knd,j,i,nf)
!��� �������������� �� ������ ����������� ����� � ������������
use pgmod
!j - ����� ��� ������, �� ������� ���� ��������������
!i - ����� ����������� ����� (������������)
!nf - ����� ���������
!ia, knd - ������ ������� � �������, ��� ��������� ������ ��������������
integer(4) ia,knd,j,i,nf
real(8) get_ba_cash_integral_
type(TCashIntegral), pointer :: c
type(TBound), pointer :: b
type(TArea), pointer :: a
type(TCashValue), pointer :: cv
!call init_ba_cash(ia,knd)
a=>gs%a(ia)
b=>a%bnd(knd)
c=>b%cash%area
if (.not.c%inited(nf)) then
  get_ba_cash_integral_=real8_inf
else
  cv=>c%intvals(nf)
  get_ba_cash_integral_=cv%i(j,i)
endif
end

function get_ba_cash_integral_2(a,knd,j,i,nf) result (res)
!��� �������������� �� ������ ����������� ����� � ������������
!���������� ������� ��� ��������
use pgmod
!j - ����� ��� ������, �� ������� ���� ��������������
!i - ����� ����������� ����� (������������)
!nf - ����� ���������
!a   - ������� ��� ��������� ������ ��������������
!knd - ����� �������, ��� ��������� ������ ��������������
integer(4) knd,j,i,nf
real(8) res
!DEC$ IF DEFINED (DEBUG)
type(TCashIntegral), pointer :: c
type(TBound), pointer :: b
type(TArea), target :: a
type(TCashValue), pointer :: cv
b=>a%bnd(knd)
c=>b%cash%area
cv=>c%intvals(nf)
res=cv%i(j,i)
!DEC$ ELSE
type(TArea) a
res=a%bnd(knd)%cash%area%intvals(nf)%i(j,i)
!DEC$ ENDIF
end

subroutine set_ba_cash_integral(ia,knd,j,i,nf,val)
!��� �������������� �� ������ ����������� ����� � ������������
use pgmod
!j - ����� ��� ������, �� ������� ���� ��������������
!i - ����� ����������� ����� (������������)
!nf - ����� ���������
!ia, knd - ������ ������� � �������, ��� ��������� ������ ��������������
integer(4) ia,knd,j,i,nf
real(8) val
!$omp critical (lock_cash)
call set_ba_cash_integral_(ia,knd,j,i,nf,val)
!$omp end critical (lock_cash)
end

subroutine set_ba_cash_integral_(ia,knd,j,i,nf,val)
!��� �������������� �� ������ ����������� ����� � ������������
use pgmod
!j - ����� ��� ������, �� ������� ���� ��������������
!i - ����� ����������� ����� (������������)
!nf - ����� ���������
!ia, knd - ������ ������� � �������, ��� ��������� ������ ��������������
integer(4) ia,knd,j,i,nf
real(8) val
type(TCashIntegral), pointer :: c
type(TBound), pointer :: b
type(TArea), pointer :: a
type(TCashValue), pointer :: cv
a=>gs%a(ia)
b=>a%bnd(knd)
c=>b%cash%area
if (.not.c%inited(nf)) call initCashIntegral_nf(c,nf,b%npanel,a%a%ntr,b%cash%ia)
cv=>c%intvals(nf)
cv%i(j,i)=val
end

subroutine set_ba_cash_integral_2(ia,knd,j,i,nf,val)
!��� �������������� �� ������ ����������� ����� � ������������
!���������� ������� ��� ��������
use pgmod
!j - ����� ��� ������, �� ������� ���� ��������������
!i - ����� ����������� ����� (������������)
!nf - ����� ���������
!ia, knd - ������ ������� � �������, ��� ��������� ������ ��������������
integer(4) ia,knd,j,i,nf
real(8) val
!DEC$ IF DEFINED (DEBUG)
type(TCashIntegral), pointer :: c
type(TBound), pointer :: b
type(TArea), pointer :: a
type(TCashValue), pointer :: cv
a=>gs%a(ia)
b=>a%bnd(knd)
c=>b%cash%area
cv=>c%intvals(nf)
cv%i(j,i)=val
!DEC$ ELSE
gs%a(ia)%bnd(knd)%cash%area%intvals(nf)%i(j,i)=val
!DEC$ ENDIF
end

function get_ab_cash_integral(ia,j,knd0,i,nf)
!��� �������������� �� ������������ ����������� ����� �� ������
use pgmod
!j - ����� ��� ������������, �� ������� ���� ��������������
!i - ����� ����������� ����� (������ � ������� ���������)
!nf - ����� ���������
!knd0 - ������ ������� � ������� � ����������� ������
!ia - ������ �������
integer(4) ia,j,knd0,i,nf
real(8) get_ab_cash_integral,get_ab_cash_integral_
if (gs_use_cash_lock_get) then
  !$omp critical (lock_cash)
  get_ab_cash_integral=get_ab_cash_integral_(ia,j,knd0,i,nf)
  !$omp end critical (lock_cash)
else
  get_ab_cash_integral=get_ab_cash_integral_(ia,j,knd0,i,nf)
endif
end

function get_ab_cash_integral_(ia,j,knd0,i,nf)
!��� �������������� �� ������������ ����������� ����� �� ������
use pgmod
!j - ����� ��� ������������, �� ������� ���� ��������������
!i - ����� ����������� ����� (������ � ������� ���������)
!nf - ����� ���������
!knd0 - ������ ������� � ������� � ����������� ������
!ia - ������ �������
integer(4) ia,j,knd0,i,nf
real(8) get_ab_cash_integral_
type(TCashIntegral), pointer :: c
type(TArea), pointer :: a
type(TCashValue), pointer :: cv
!call init_ab_cash(ia)
a=>gs%a(ia)
c=>a%a%cash%bnd(knd0)
if (.not.c%inited(nf)) then
  get_ab_cash_integral_=real8_inf
else
  cv=>c%intvals(nf)
  get_ab_cash_integral_=cv%i(j,i)
endif
end

function get_ab_cash_integral_2(a,j,knd0,i,nf) result(res)
!��� �������������� �� ������������ ����������� ����� �� ������
!���������� ������� ��� ��������
use pgmod
!j - ����� ��� ������������, �� ������� ���� ��������������
!i - ����� ����������� ����� (������ � ������� ���������)
!nf - ����� ���������
!knd0 - ����� ������� � ����������� ������
!a - �������
integer(4) j,knd0,i,nf
real(8) res
!DEC$ IF DEFINED (DEBUG)
type(TCashIntegral), pointer :: c
type(TArea), target :: a
type(TCashValue), pointer :: cv
c=>a%a%cash%bnd(knd0)
cv=>c%intvals(nf)
res=cv%i(j,i)
!DEC$ ELSE
type(TArea) a
res=a%a%cash%bnd(knd0)%intvals(nf)%i(j,i)
!DEC$ ENDIF
end

subroutine set_ab_cash_integral(ia,j,knd0,i,nf,val)
!��� �������������� �� ������������ ����������� ����� �� ������
use pgmod
!j - ����� ��� ������������, �� ������� ���� ��������������
!i - ����� ����������� ����� (������ � ������� ���������)
!nf - ����� ���������
!knd0 - ������ ������� � ������� � ����������� ������
!ia - ������ �������
integer(4) ia,j,knd0,i,nf
real(8) val
!$omp critical (lock_cash)
call set_ab_cash_integral_(ia,j,knd0,i,nf,val)
!$omp end critical (lock_cash)
end

subroutine set_ab_cash_integral_(ia,j,knd0,i,nf,val)
!��� �������������� �� ������������ ����������� ����� �� ������
use pgmod
!j - ����� ��� ������������, �� ������� ���� ��������������
!i - ����� ����������� ����� (������ � ������� ���������)
!nf - ����� ���������
!knd0 - ������ ������� � ������� � ����������� ������
!ia - ������ �������
integer(4) ia,j,knd0,i,nf
real(8) val
type(TCashIntegral), pointer :: c
type(TBound), pointer :: b0
type(TArea), pointer :: a
type(TCashValue), pointer :: cv
a=>gs%a(ia)
c=>a%a%cash%bnd(knd0)
b0=>a%bnd(knd0)
if (.not.c%inited(nf)) call initCashIntegral_nf(c,nf,a%a%ntr,b0%npanel,a%a%cash%ia)
cv=>c%intvals(nf)
cv%i(j,i)=val
end

subroutine set_ab_cash_integral_2(ia,j,knd0,i,nf,val)
!��� �������������� �� ������������ ����������� ����� �� ������
!���������� ������� ��� ��������
use pgmod
!j - ����� ��� ������������, �� ������� ���� ��������������
!i - ����� ����������� ����� (������ � ������� ���������)
!nf - ����� ���������
!knd0 - ������ ������� � ������� � ����������� ������
!ia - ������ �������
integer(4) ia,j,knd0,i,nf
real(8) val
!DEC$ IF DEFINED (DEBUG)
type(TCashIntegral), pointer :: c
type(TBound), pointer :: b0
type(TArea), pointer :: a
type(TCashValue), pointer :: cv
a=>gs%a(ia)
c=>a%a%cash%bnd(knd0)
b0=>a%bnd(knd0)
cv=>c%intvals(nf)
cv%i(j,i)=val
!DEC$ ELSE
gs%a(ia)%a%cash%bnd(knd0)%intvals(nf)%i(j,i)=val
!DEC$ ENDIF
end

function get_aa_cash_integral(ia,j,i,nf)
!��� �������������� �� ������������ ����������� ����� �� ������������
use pgmod
!j - ����� ��� ������������, �� ������� ���� ��������������
!i - ����� ����������� ����� (������������)
!nf - ����� ���������
!ia - ������ �������
integer(4) ia,j,i,nf
real(8) get_aa_cash_integral,get_aa_cash_integral_
if (gs_use_cash_lock_get) then
  !$omp critical (lock_cash)
  get_aa_cash_integral=get_aa_cash_integral_(ia,j,i,nf)
  !$omp end critical (lock_cash)
else
  get_aa_cash_integral=get_aa_cash_integral_(ia,j,i,nf)
endif
end

function get_aa_cash_integral_(ia,j,i,nf)
!��� �������������� �� ������������ ����������� ����� �� ������������
use pgmod
!j - ����� ��� ������������, �� ������� ���� ��������������
!i - ����� ����������� ����� (������������)
!nf - ����� ���������
!ia - ������ �������
integer(4) ia,j,i,nf
real(8) get_aa_cash_integral_
type(TCashIntegral), pointer :: c
type(TArea), pointer :: a
type(TCashValue), pointer :: cv
!call init_aa_cash(ia)
a=>gs%a(ia)
c=>a%a%cash%area
if (.not.c%inited(nf)) then
  get_aa_cash_integral_=real8_inf
else
  cv=>c%intvals(nf)
  get_aa_cash_integral_=cv%i(j,i)
endif
end

function get_aa_cash_integral_2(a,j,i,nf) result (res)
!��� �������������� �� ������������ ����������� ����� �� ������������
!���������� ������� ��� ��������
use pgmod
!j - ����� ��� ������������, �� ������� ���� ��������������
!i - ����� ����������� ����� (������������)
!nf - ����� ���������
!a - �������
integer(4) j,i,nf
real(8) res
!DEC$ IF DEFINED (DEBUG)
type(TArea), target :: a
type(TCashValue), pointer :: cv
cv=>a%a%cash%area%intvals(nf)
res=cv%i(j,i)
!DEC$ ELSE
type(TArea) a
res=a%a%cash%area%intvals(nf)%i(j,i)
!DEC$ ENDIF
end

subroutine set_aa_cash_integral(ia,j,i,nf,val)
!��� �������������� �� ������������ ����������� ����� �� ������������
use pgmod
!j - ����� ��� ������������, �� ������� ���� ��������������
!i - ����� ����������� ����� (������������)
!nf - ����� ���������
!ia - ������ �������
integer(4) ia,j,i,nf
real(8) val
!$omp critical (lock_cash)
call set_aa_cash_integral_(ia,j,i,nf,val)
!$omp end critical (lock_cash)
end

subroutine set_aa_cash_integral_(ia,j,i,nf,val)
!��� �������������� �� ������������ ����������� ����� �� ������������
use pgmod
!j - ����� ��� ������������, �� ������� ���� ��������������
!i - ����� ����������� ����� (������������)
!nf - ����� ���������
!ia - ������ �������
integer(4) ia,j,i,nf
real(8) val
type(TCashIntegral), pointer :: c
type(TArea), pointer :: a
type(TCashValue), pointer :: cv
a=>gs%a(ia)
c=>a%a%cash%area
if (.not.c%inited(nf)) call initCashIntegral_nf(c,nf,a%a%ntr,a%a%ntr,a%a%cash%ia)
cv=>c%intvals(nf)
cv%i(j,i)=val
end

subroutine set_aa_cash_integral_2(ia,j,i,nf,val)
!��� �������������� �� ������������ ����������� ����� �� ������������
!���������� ������� ��� ��������
use pgmod
!j - ����� ��� ������������, �� ������� ���� ��������������
!i - ����� ����������� ����� (������������)
!nf - ����� ���������
!ia - ������ �������
integer(4) ia,j,i,nf
real(8) val
!DEC$ IF DEFINED (DEBUG)
type(TArea), pointer :: a
type(TCashValue), pointer :: cv
a=>gs%a(ia)
cv=>a%a%cash%area%intvals(nf)
cv%i(j,i)=val
!DEC$ ELSE
gs%a(ia)%a%cash%area%intvals(nf)%i(j,i)=val
!DEC$ ENDIF
end

subroutine deallocate_cash(cash)
use pgmod
integer(4) i
type(TCash) cash
if (cash%bnd_inited) then
  do i=1,ubound(cash%bnd,dim=1)
    call deallocate_cashIntegrals(cash%bnd(i))
  enddo
  deallocate(cash%bnd)
endif
if (cash%area_inited) call deallocate_cashIntegrals(cash%area)
end

subroutine deallocate_cashIntegrals(c)
use pgmod
type(TCashIntegral) c
integer(4) i
do i=0,max_int
  if (c%inited(i)) deallocate(c%intvals(i)%i)
enddo
deallocate(c%inited)
deallocate(c%solved)
deallocate(c%intvals)
end

subroutine pg_set_cash_ref(ip,ia)
!dec$ attributes dllexport:: pg_set_cash_ref
!���������� ������ �� cash ������ �������, �� ������� ������� ������� �������� ������������ ������� � ����� � ������� �� ��������� � ���������
use pgmod
integer(4) ip !������ ������ � ��������, �� ������� ����� ������
integer(4) ia !������ �������, �� ������� ����� ������
integer(4) i
type(TArea), pointer :: a2
type(TBound), pointer :: b,b2
if (gsarea%type_rotate<0) then
  call gs_print("Error pg_set_cash_ref. Set cash ref for an uncopied domain (ia="//itoa(ia)//")!")
  call gs_print_stop("Use pg_copy_subdomain_dxdy_rotate");
else if (gsarea%type_rotate==0.and.(.not.gsarea%oss_resolve)) then
  !������� (����� ��������������� ����� � ��������� �� y)
  a2=>gsMain%ggs(ip)%a(ia)
  do i=1,gsarea%nb
    b=>gsarea%bnd(i)
    b2=>a2%bnd(i)
    b%cash=>b2%self_cash
  enddo
  gsarea%a%cash=>a2%a%self_cash
endif
end

subroutine gs_print(s)
use pgmod
character(*) s
print*, trim(s)
if (gs_write_log_file) write(gs_log_file_i,*) trim(s)
end

subroutine gs_print_stop(s)
use pgmod
character(*) s
call gs_print(s)
!DEC$ IF DEFINED (DEBUG)
call system("pause")
!DEC$ ENDIF
close(gs_log_file_i)
stop
end

subroutine pg_clear_skip_allbounds
!dec$ attributes dllexport:: pg_set_cash_ref
!�������� �������� skip_getfun � ���� ������ ������� �������
use pgmod
integer(4) i
do i=1,gsarea%nb
  gsarea%bnd(i)%skip_getfun=.false.
enddo
end

function get_gug_c(av,k) result(cc0)
!�������� i-� ����������� (0 - �0) ��� k-� ������
use pgmod
type(TAreaValue_c) av
integer(4) k
real(8) cc0
if (allocated(av%v)) then
  cc0=av%v(k)
else
  cc0=av%c
endif
end
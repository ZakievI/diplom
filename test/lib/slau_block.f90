!блочно диагональная матрица
!блоки - плотная матрица, остальная часть - разреженная
module slau_block
use gen_mod
!use slau_int   

  type sparse_matrix_row
      integer(4) nc                 !количество ненулевых коэффициентов
      real(8), allocatable :: m(:)  !коэффициенты
      integer(4), allocatable :: c(:)  !индексы столбцов
      integer(4) sm_ibegin !начальный индекс в sparse%m при выделении сразу глобальной матрицы
      integer(4) sm_count  !количество неизвестных в области - для предварительного выделения памяти
    endtype
    
  type sparse_matrix
    integer(4) nc                 !количество ненулевых коэффициентов
    real(8), allocatable :: m(:)  !коэффициенты
    integer(4), allocatable :: c(:)  !индексы столбцов
    integer(4), allocatable :: r(:) !индексы элементов, начинающих новые строки
    integer(4), allocatable :: r2(:) !индексы элементов, начинающих новые строки после блока
    type(sparse_matrix_row), allocatable :: mrows(:)
  endtype
  
  type block_matrix
    integer(4) nu !количество уравнений в блоке
    integer(4) nx !количество неизвестных в блоке
    integer(4) iu !индекс первого уравнения в глобальной матрице
    integer(4) iu2 !индекс последнего уравнения в глобальной матрице
    integer(4) ix !индекс первой неизвестной блока в глобальной матрице
    integer(4) ix2 !номер последнего столбца неизвестных в глобальной матрице
    real(8), allocatable :: m(:,:) !диагональный блок (плотная матрица)
    type(sparse_matrix) sparse     !коэффициенты вне диагонального блока (разреженная матрица)
    integer(4), allocatable :: p(:) !столбец перестановок
  endtype
  
  type main_block_matrix
    integer(4) nblock  !количество блоков
    integer(4) nu_all  !общее число уравнений
    integer(4) nx_all  !общее число неизвестных
    type(block_matrix), allocatable :: blockm(:) !(nblock) блоки
    integer(4), allocatable :: ord (:) !(nblock)порядок блоков в итерационном процессе
    integer(4), allocatable :: klevel (:) !(nlevel+1) индексы блоков, с которых начинаются уровни
    integer(4) nlevel !количество уровней
    integer(4), allocatable :: eq_block(:) !(nu_all) номер блока для уравнения
    logical r2_inited !инициализированы массивы blockm(i)%sparse%r2
    logical p_inited !инициализированы массивы blockm(i)%p
    type(main_block_matrix), pointer :: bm_ref
    real(8), allocatable :: x_limit(:,:) !(2,nx_all)
    logical have_limit
  endtype
  
  integer(4) sb_write_debug_res
  procedure(), pointer :: sb_iter_res=>null()
  logical sb_solver_write_all_iteration !показывать все итерации в методе решения блочных СЛАУ
  real(8) sb_solver_erra_max
  
  real(8), parameter :: sp_empty_real = -0.99925d10
  
  interface

  subroutine bm_main_block_sparse(bm)
  Import :: main_block_matrix
  type(main_block_matrix), target :: bm
  end
    
  subroutine bm_calc_slau_block_diagonal2(bm,rsh,x,lam,Imax,eps,maxIter_convergence,system_i,write_iter)
  Import :: main_block_matrix
  type(main_block_matrix), target :: bm
  real(8) rsh(bm%nu_all)
  real(8), target :: x(bm%nx_all)
  real(8) lam,eps
  integer(4) system_i,Imax,maxIter_convergence
  logical write_iter
  external write_iter
  end
  
  subroutine bm_init_block_p(bm)
  Import :: main_block_matrix
  type(main_block_matrix), target :: bm
  end
  
  subroutine bm_A_mult_x2(alf,bm,x,bet,bb,res,buff)
  !res=alf*bm*x+bet*bb
  Import :: main_block_matrix
  type(main_block_matrix), target :: bm
  real(8) alf,bet
  real(8) x(:),res(:),buff(:),bb(:)
  end
  
  subroutine bm_ILU_sum(bmILU,b,bILU,k1,j,imax,sum)
  !подсчитать сумму -\sum_1^imax(l_ki*u_ij)
  Import :: main_block_matrix,block_matrix
  type(main_block_matrix), target :: bmILU
  type(block_matrix), pointer :: b,bILU
  integer(4) j,k1,imax
  real(8) sum
  end
  
  subroutine bm_swap_vector_P(bm,x)
  !переставить элементы вектора в соответствии с с массивами b%p
  Import :: main_block_matrix
  type(main_block_matrix), target :: bm !исходная матрица, содержащая портрет (массивы sparse%r и sparse%c в блоках)
  real(8) x(:) !(bm%nu_all)
  end
  
  subroutine bm_init_block_r2(bm)
  Import :: main_block_matrix
  type(main_block_matrix), target :: bm
  end
  
  subroutine bm_init_block_matrix(bm)
  Import :: main_block_matrix
  type(main_block_matrix), target :: bm
  end
  
  subroutine bm_ILUP_sum(bmILU,b,bILU,k1,j,imax,sum)
  !подсчитать сумму -\sum_1^imax(l_ki*u_ij)
  Import :: main_block_matrix,block_matrix
  type(main_block_matrix), target :: bmILU
  type(block_matrix), pointer :: b,bILU
  integer(4) j,k1,imax
  real(8) sum
  end
  
  subroutine bm_init_ILU(bm,bmILU,init_p)
  Import :: main_block_matrix
  type(main_block_matrix), target :: bm !исходная матрица
  type(main_block_matrix), target :: bmILU !факторизованная матрица
  logical init_p
  end
  
  subroutine bm_init_block_p_ipiv(bm,ipiv)
  Import :: main_block_matrix
  type(main_block_matrix), target :: bm
  integer(4) ipiv(:) !(bm%nu_all)
  end
  
  subroutine bm_ILUP(bmILU)
  !построение ILU-факторизации для блочно-диагональной матрицы
  Import :: main_block_matrix
  type(main_block_matrix), target :: bmILU !факторизованная матрица
  end
  
  subroutine bm_ILUP_solve(bmILU,x)
  !вычисление x=M^-1*b, где M - LU матрица (полученная, например, при ILU факторизации)
  !!!для случая, когда диагональные элементы принадлежат блокам!
  !Ly=b, Ux=y
  Import :: main_block_matrix
  type(main_block_matrix), target :: bmILU !факторизованная матрица
  real(8) x(:) !(bm%nu_all) на входе правая часть, на выходе - решение
  end
  
  subroutine bm_A_mult_x(alf,bm,x,res,buff)
  !res=alf*bm*x
  Import :: main_block_matrix
  type(main_block_matrix) bm
  real(8) alf
  real(8) x(:),res(:),buff(:)
  end
  
  subroutine bm_m_to_block_diagonal(m,bm)
  Import :: main_block_matrix
  real(8) m(:,:)
  type(main_block_matrix), target :: bm
  end
  
  subroutine bm_BiCGStab(bm,rsh,x,Imax,eps,maxIter_convergence,write_iter)
  !https://ru.wikipedia.org/wiki/Стабилизированный_метод_бисопряжённых_градиентов
  Import :: main_block_matrix
  type(main_block_matrix) bm
  real(8) rsh(bm%nu_all)
  real(8), target :: x(bm%nx_all)
  integer(4) Imax,maxIter_convergence
  real(8) eps
  logical write_iter
  external write_iter
  end
  
  subroutine bm_ILU(bmILU)
  !построение ILU-факторицации для блочно-диагональной матрицы
  Import :: main_block_matrix
  type(main_block_matrix), target :: bmILU !факторизованная матрица
  end
  
  subroutine bm_ILU_solve(bm,bmILU,x)
  !вычисление x=M^-1*b, где M - LU матрица (полученная, например, при ILU факторизации)
  !!!для случая, когда диагональные элементы принадлежат блокам!
  !Ly=b, Ux=y
  Import :: main_block_matrix
  type(main_block_matrix), target :: bm !исходная матрица, содержащая портрет (массивы sparse%r и sparse%c в блоках)
  type(main_block_matrix), target :: bmILU !факторизованная матрица
  real(8) x(:) !(bm%nu_all) на входе правая часть, на выходе - решение
  end

  
  end interface
  
end module

subroutine sb_init
use slau_block
sb_write_debug_res=0
sb_solver_write_all_iteration=.false.
sb_solver_erra_max=1.0d20
end

subroutine bm_allocate_blocks(bm,nb)
use slau_block
type(main_block_matrix) bm
integer(4) nb,i
if (bm%nblock==0) then
  allocate(bm%blockm(nb))
  allocate(bm%ord(nb))
  bm%nlevel=1
  allocate(bm%klevel(bm%nlevel+1))
  do i=1,nb
    call null_block(bm%blockm(i))
    bm%ord(i)=i
  enddo
  bm%klevel(1)=1
  bm%klevel(2)=nb+1
  bm%nblock=nb
endif
bm%nu_all=0
bm%nx_all=0
bm%have_limit=.false.
end

subroutine null_block(b)
use slau_block
type(block_matrix) b
b%nu=0
b%nx=0
b%iu=0
b%ix=0
end

subroutine bm_deallocate(need_mm,bm)
use slau_block
logical need_mm
type(main_block_matrix) bm
integer(4) i
do i=1,bm%nblock
  call bm_deallocate_sparse_m(bm%blockm(i)%sparse)
enddo
bm%r2_inited=.false.
bm%p_inited=.false.
if (need_mm) then
  do i=1,bm%nblock
    call bm_deallocate_block(bm%blockm(i))
  enddo
  if (allocated(bm%blockm)) deallocate(bm%blockm)
  if (allocated(bm%ord)) deallocate(bm%ord)
  if (allocated(bm%klevel)) deallocate(bm%klevel)
  bm%nblock=0
  if (allocated(bm%eq_block)) deallocate(bm%eq_block)
  if (allocated(bm%x_limit)) deallocate(bm%x_limit)
endif
end

subroutine bm_deallocate_block(b)
use slau_block
type(block_matrix) b
if (allocated(b%m)) deallocate(b%m)
if (allocated(b%p)) deallocate(b%p)
call bm_deallocate_sparse_mrows(b%sparse)
b%nu=0
b%nx=0
end

subroutine bm_reallocate_block(b,nu,nx)
use slau_block
integer(4) nu,nx,j
type(block_matrix) b
if (b%nx/=nx.or.b%nu/=nu) then
  call bm_deallocate_block(b)
  allocate(b%m(nu,nx))
  allocate(b%sparse%mrows(nu))
  b%nu=nu
  b%nx=nx
endif
!b%m=d0
do j=1,b%nu
  call bm_null_mrow(b%sparse%mrows(j))
enddo
end

subroutine bm_deallocate_sparse_m(sparse)
use slau_block
type(sparse_matrix) sparse
if (allocated(sparse%m)) deallocate(sparse%m)
if (allocated(sparse%c)) deallocate(sparse%c)
if (allocated(sparse%r)) deallocate(sparse%r)
if (allocated(sparse%r2)) deallocate(sparse%r2)
end

subroutine bm_null_mrow(mrow)
use slau_block
type(sparse_matrix_row) mrow
if (allocated(mrow%m)) deallocate(mrow%m)
if (allocated(mrow%c)) deallocate(mrow%c)
mrow%nc=0
mrow%sm_count=0
end

subroutine bm_deallocate_sparse_mrows(sparse)
use slau_block
type(sparse_matrix) sparse
integer(4) i
if (allocated(sparse%mrows)) then
  do i=1,ubound(sparse%mrows,1)
    call bm_null_mrow(sparse%mrows(i))
  enddo
  deallocate(sparse%mrows)
endif
end

subroutine bm_init_block_matrix(bm)
use slau_block
integer(4) i,nu,ku,kx,nx
type(main_block_matrix), target :: bm
type(block_matrix), pointer :: b
nu=0
nx=0
do i=1,bm%nblock
  b=>bm%blockm(i)
  nu=nu+b%nu
  nx=nx+b%nx
enddo
bm%nu_all=nu
bm%nx_all=nx
if (allocated(bm%eq_block)) deallocate(bm%eq_block)
allocate(bm%eq_block(bm%nu_all))
ku=1
kx=1
do i=1,bm%nblock
  b=>bm%blockm(i)
  b%iu=ku
  b%ix=kx
  ku=ku+b%nu
  kx=kx+b%nx
  b%ix2=kx-1
  b%iu2=ku-1
  bm%eq_block(b%iu:b%iu2)=i
enddo
bm%r2_inited=.false.
bm%p_inited=.false.
end

subroutine bm_main_block_sparse(bm)
use slau_block
type(main_block_matrix), target :: bm
integer(4) ia,i,j,j1
type(block_matrix), pointer :: b
type(sparse_matrix_row), pointer :: mrow
do ia=1,bm%nblock
  b=>bm%blockm(ia)
  b%sparse%nc=0
  do i=1,b%nu
    b%sparse%nc=b%sparse%nc+b%sparse%mrows(i)%nc
  enddo
  allocate(b%sparse%m(b%sparse%nc))
  allocate(b%sparse%c(b%sparse%nc))
  allocate(b%sparse%r(b%nu+1))
  j=1
  do i=1,b%nu
    mrow=>b%sparse%mrows(i)
    j1=j+mrow%nc-1
    b%sparse%m(j:j1)=mrow%m(1:mrow%nc)
    b%sparse%c(j:j1)=mrow%c(1:mrow%nc)
    b%sparse%r(i)=j
    j=j1+1
    call bm_null_mrow(mrow)
  enddo
  b%sparse%r(b%nu+1)=j
enddo
end

subroutine bm_init_block_r2(bm)
use slau_block
type(main_block_matrix), target :: bm
type(block_matrix), pointer :: b
integer(4) ia,i,j,k,j1,j2
if (bm%r2_inited) return
do ia=1,bm%nblock
  b=>bm%blockm(ia)
  allocate(b%sparse%r2(b%nu))
  do i=1,b%nu
    j1=b%sparse%r(i)
    j2=b%sparse%r(i+1)-1
    b%sparse%r2(i)=b%sparse%r(i+1)
    do j=j1,j2
      k=b%sparse%c(j)
      if (k>b%ix) then
        b%sparse%r2(i)=j
        exit
      endif
    enddo
  enddo
enddo
bm%r2_inited=.true.
end

subroutine bm_init_block_p(bm)
use slau_block
type(main_block_matrix), target :: bm
type(block_matrix), pointer :: b
integer(4) ia,i
if (bm%p_inited) return
do ia=1,bm%nblock
  b=>bm%blockm(ia)
  allocate(b%p(b%nu))
  forall (i=1:b%nu) b%p(i)=i
enddo
bm%p_inited=.true.
end 

subroutine bm_init_block_p_ipiv(bm,ipiv)
use slau_block
type(main_block_matrix), target :: bm
type(block_matrix), pointer :: b
integer(4) ipiv(:) !(bm%nu_all)
integer(4) ia,i,i1
if (bm%p_inited) return
call bm_init_block_p(bm)
do ia=1,bm%nblock
  b=>bm%blockm(ia)
  i1=0
  do i=b%iu,b%iu2
    i1=i1+1
    if (ipiv(i)/=i1) call swap_int(b%p(i1),b%p(ipiv(i)))
  enddo
enddo
bm%p_inited=.true.
end

subroutine bm_calc_slau_block_diagonal2(bm,rsh,x,lam,Imax,eps,maxIter_convergence,system_i,write_iter)
use slau_block
type(main_block_matrix), target :: bm
integer(4) ia,ia1,info,i,i1,i2,j,j1,j2,system_i,Imax,maxIter_convergence
real(8) norm,dnrm2,t,norma
logical is_x_ne_0,test_array_ne_0
type(block_matrix), pointer :: b
integer(4), allocatable :: ipiv(:)
real(8), allocatable :: b2(:),b1(:),normx(:)
real(8), allocatable, target :: xold(:)
real(8) rsh(bm%nu_all)
real(8), target :: x(bm%nx_all)
character(80) str
real(8) lam,lam1,best_err,eps !,lam0
integer ibest_err
logical last_iter_writed,write_iter,is_maxIter_convergence
external write_iter
type(main_block_matrix) bmILU
logical use_precond,need_limit
!real(8) prev_norm
!integer(4) norm_less_count
!norm_less_count=0

use_precond=.false.
if (use_precond) call bm_init_ILU(bm,bmILU,.false.)

!lam0=lam
!lam=0.05d0
lam1=d1-lam
!делаем LU разложение для блоков
allocate(ipiv(bm%nu_all))
do ia=1,bm%nblock
  system_i=ia
  b=>bm%blockm(ia)
  if (b%nx==b%nu) then
    call dgetrf(b%nu, b%nx, b%m, b%nu, ipiv(b%iu), info)
  else
  endif
enddo
if (use_precond) then
  call bm_init_block_p_ipiv(bm,ipiv)
  call bm_ILUP(bmILU)
  x=rsh
  call bm_ILUP_solve(bmILU,x)
  call bm_deallocate(.true.,bmILU)
endif
!итерационный процесс
need_limit=.true.
is_maxIter_convergence=.false.
allocate(b1(bm%nu_all),xold(bm%nx_all),b2(bm%nu_all),normx(bm%nx_all))
ibest_err=0
write(str,"('Iteration block matrix. MaxIter=',i0,', iterConv=',i0,', eps=',E9.2)") Imax,maxIter_convergence,eps
last_iter_writed=write_iter(trim(str),.false.)
do i=1,Imax
  xold=x
  do ia1=1,bm%nblock
    ia=bm%ord(ia1)
    b=>bm%blockm(ia)
    i1=b%iu
    i2=b%iu2
    b1(i1:i2)=rsh(i1:i2)
    is_x_ne_0=.false.
    do j=1,2
      if (j==1) then !перед блоком
        j1=1
        j2=b%ix-1
      else           !после блока
        j1=b%ix+b%nx
        j2=bm%nx_all
      endif
      if (j2<j1) cycle
      !проверяем, что в текущем решении есть ненулевые значения
      is_x_ne_0=test_array_ne_0(x(j1),j2-j1+1)
      if (is_x_ne_0) exit
    enddo
    !переносим переменные вне блока в правую часть с учетом текущего приближения для x
    if (is_x_ne_0) then
      call mkl_dcsrgemv('N',b%nu,b%sparse%m,b%sparse%r,b%sparse%c,x,b2)
      forall (j=i1:i2) b1(j)=b1(j)-b2(j-i1+1)
    endif
    !проверяем, что правая часть не нулевая (иначе будет нулевое решение)
    is_x_ne_0=test_array_ne_0(b1(i1),b%nu)
    !решаем СЛАУ для блока с подправленной правой частью
    if (is_x_ne_0) then
      x(i1:i2)=b1(i1:i2)
      call dgetrs( 'N', b%nu, 1, b%m, b%nu, ipiv(i1), x(i1), b%nu, info )
    else
      x(i1:i2)=d0
    endif
    if (bm%have_limit.and.need_limit) then
      do j=i1,i2
        if (bm%x_limit(1,j)/=sp_empty_real) x(j)=max(bm%x_limit(1,j),x(j))
        if (bm%x_limit(2,j)/=sp_empty_real) x(j)=min(bm%x_limit(2,j),x(j))
      enddo
    endif
    if (bm%nblock>1) forall (j=i1:i2) x(j)=lam*x(j)+lam1*xold(j)
    
    if (sb_write_debug_res==2) call sb_iter_res
    
  enddo
  if (bm%nblock>1) then
    !вычисляем невязку
    normx=x-xold
    norma=dnrm2(bm%nx_all, normx, 1)
    do j=1,bm%nx_all
      t=dabs(x(j))
      if (t>d1) normx(j)=normx(j)/t
    enddo
    norm=dnrm2(bm%nx_all, normx, 1)
  else
    norm=d0
    norma=d0
  endif
  write(str,"('i=',i0,' err=',E11.3,' erra=',E11.3)") i,norm,norma
  !last_iter_writed=write_iter(trim(str),.true.)
  last_iter_writed=write_iter(trim(str),.not.sb_solver_write_all_iteration)
  
  if (norm<eps) exit
  !need_limit=norm>1.0d-1
  
  !if (norm<1.0d-2.and.lam<lam0) then
  !  write(*,*)"lam!"
  !  lam=lam0
  !  lam1=d1-lam
  !endif
  
  if (sb_write_debug_res==1) call sb_iter_res
  
  if (norma>sb_solver_erra_max) exit
  
  if ((norm<d1.and.ibest_err==0).or.norm<best_err) then
    best_err=norm
    ibest_err=i
  else
    is_maxIter_convergence=i-ibest_err>maxIter_convergence.and.ibest_err>0
    if (is_maxIter_convergence) exit
  endif
  
  !if (i>1) then
  !  if (norm<prev_norm) then
  !    norm_less_count=norm_less_count+1
  !    if (norm_less_count>10) then
  !      lam=lam+0.002d0
  !      lam1=d1-lam
  !      write(str,"('new lam=',E11.3)") lam
  !      last_iter_writed=write_iter(trim(str),.true.)
  !      norm_less_count=0
  !    endif
  !  else
  !    norm_less_count=0
  !  endif
  !endif
  !prev_norm=norm
  
enddo
if (.not.last_iter_writed) last_iter_writed=write_iter(str,.false.)
deallocate(ipiv,b1,xold,b2,normx)
if (i>Imax.or.norma>sb_solver_erra_max) then
  last_iter_writed=write_iter("Iterative process does not converge! (bm_calc_slau_block_diagonal)",.false.)
  if (i>Imax) then
    last_iter_writed=write_iter("MaxIter reached!!!",.false.)
  else
    last_iter_writed=write_iter("Max Abs Error reached!!!",.false.)
  endif
  !DEC$ IF DEFINED (DEBUG)
  call system("pause")
  !DEC$ ENDIF
  stop
endif
if (is_maxIter_convergence) last_iter_writed=write_iter("MaxIter_convergence reached!",.false.)
end

function write_iter_empty(str,test)
character(*) str
logical write_iter_empty,test
write(*,*) str
write_iter_empty=.true.
test=test
end

subroutine bm_calc_slau_block_diagonal(bm,rsh,x,lam,Imax,eps,maxIter_convergence)
use slau_block
type(main_block_matrix), target :: bm
real(8) x(bm%nx_all),rsh(bm%nu_all)
real(8) lam,eps
integer(4) Imax,maxIter_convergence
integer(4) system_i
logical write_iter_empty
external write_iter_empty
call bm_calc_slau_block_diagonal2(bm,rsh,x,lam,Imax,eps,maxIter_convergence,system_i,write_iter_empty)
end

subroutine bm_BiCGStab(bm,rsh,x,Imax,eps,maxIter_convergence,write_iter)
!https://ru.wikipedia.org/wiki/Стабилизированный_метод_бисопряжённых_градиентов
use slau_block
type(main_block_matrix) bm,bmILU
integer(4) Imax,maxIter_convergence
real(8) eps,ddot
real(8) rsh(bm%nu_all)
real(8), target :: x(bm%nx_all)
external write_iter
real(8), allocatable :: r0(:),r(:),v(:),p(:),s(:),t(:),buff(:),y(:),z(:)
real(8) ro,ro0,alf,bet,om,om_,normb,normr,dnrm2,norm
integer(4) i,j
character(80) str
real(8) best_err
integer ibest_err
logical last_iter_writed,is_maxIter_convergence,write_iter,use_procond
use_procond=.true.
is_maxIter_convergence=.false.
ibest_err=0
write(str,"('BiCGStab block matrix. MaxIter=',i0,', iterConv=',i0,', eps=',E9.2)") Imax,maxIter_convergence,eps
last_iter_writed=write_iter(trim(str),.false.)
allocate(r0(bm%nu_all),r(bm%nu_all),v(bm%nu_all),p(bm%nu_all),s(bm%nu_all),t(bm%nu_all),buff(bm%nu_all))
normb=dnrm2(bm%nu_all, rsh, 1)
if (use_procond) then
  call bm_init_ILU(bm,bmILU,.true.)
  call bm_ILUP(bmILU)
  allocate(y(bm%nu_all),z(bm%nu_all))
endif
!подготовка итерационного процесса
call bm_A_mult_x2(-d1,bm,x,d1,rsh,r0,buff) !0.2
r=r0  !0.3
ro=d1 !0.4
alf=d1
om=d1
v=d0 !0.5
p=d0
do i=1,Imax
  ro0=ro
  ro=ddot(bm%nu_all,r0,1,r,1) !1
  bet=ro/ro0*alf/om            !2
  forall (j=1:bm%nu_all) p(j)=r(j)+bet*(p(j)-om*v(j)) !3
  if (use_procond) then
    y=p
    call bm_ILUP_solve(bmILU,y)
    call bm_A_mult_x(d1,bm,y,v,buff) !4
  else
    call bm_A_mult_x(d1,bm,p,v,buff) !4
  endif
  alf=ddot(bm%nu_all,r0,1,v,1) !5
  alf=ro/alf
  forall (j=1:bm%nu_all) s(j)=r(j)-alf*v(j) !6
  if (use_procond) then
    z=s
    call bm_ILUP_solve(bmILU,z)
    call bm_A_mult_x(d1,bm,z,t,buff) !7
  else
    call bm_A_mult_x(d1,bm,s,t,buff) !7
  endif
  om=ddot(bm%nu_all,t,1,s,1) !8
  om_=ddot(bm%nu_all,t,1,t,1) 
  om=om/om_
  if (use_procond) then
    forall (j=1:bm%nu_all) x(j)=x(j)+om*z(j)+alf*y(j) !9
  else
    forall (j=1:bm%nu_all) x(j)=x(j)+om*s(j)+alf*p(j) !9
  endif
  forall (j=1:bm%nu_all) r(j)=s(j)-om*t(j) !10
  !критерий сходимости
  normr=dnrm2(bm%nu_all, r, 1)
  norm=normr/normb
  write(str,"('i=',i0,' err=',E11.3)") i,norm
  last_iter_writed=write_iter(trim(str),.true.)
  if (norm<eps) exit
  if ((norm<d1.and.ibest_err==0).or.norm<best_err) then
    best_err=norm
    ibest_err=i
  else
    is_maxIter_convergence=i-ibest_err>maxIter_convergence.and.ibest_err>0
    if (is_maxIter_convergence) exit
  endif
enddo
if (.not.last_iter_writed) last_iter_writed=write_iter(str,.false.)
deallocate(r0,r,v,p,s,t,buff)
if (use_procond) then
  call bm_deallocate(.true.,bmILU)
  deallocate(y,z)
endif
if (i>Imax) then
  last_iter_writed=write_iter("Iterative process does not converge! (bm_calc_slau_block_diagonal)",.false.)
  last_iter_writed=write_iter("MaxIter reached!!!",.false.)
  !DEC$ IF DEFINED (DEBUG)
  call system("pause")
  !DEC$ ENDIF
  stop
endif
if (is_maxIter_convergence) last_iter_writed=write_iter("MaxIter_convergence reached!",.false.)
end

subroutine bm_A_mult_x(alf,bm,x,res,buff)
!res=alf*bm*x
use slau_block
type(main_block_matrix) bm
real(8) alf
real(8) x(:),res(:),buff(:)
call bm_A_mult_x2(alf,bm,x,d0,buff,res,buff)
end

subroutine bm_A_mult_x2(alf,bm,x,bet,bb,res,buff)
!res=alf*bm*x+bet*bb
use slau_block
type(main_block_matrix), target :: bm
type(block_matrix), pointer :: b
real(8) alf,bet
real(8) x(:),res(:),buff(:),bb(:)
integer(4) i,j
if (bet/=d0) res=bb
do i=1,bm%nblock
  b=>bm%blockm(i)
  call dgemv('N', b%nu, b%nx, alf, b%m, b%nu, x(b%ix), 1, bet, res(b%iu), 1)
  call mkl_dcsrgemv('N',b%nu,b%sparse%m,b%sparse%r,b%sparse%c,x,buff(b%iu))
  if (alf/=d1) then
    forall(j=b%iu:b%iu2) res(j)=res(j)+alf*buff(j)
  else
    forall(j=b%iu:b%iu2) res(j)=res(j)+buff(j)
  endif
enddo
end

subroutine bm_m_to_block_diagonal(m,bm)
use slau_block
real(8) m(:,:)
type(main_block_matrix), target :: bm
type(block_matrix), pointer :: b
integer(4) i,j,k,n
do i=1,bm%nblock
  b=>bm%blockm(i)
  b%m=m(b%iu:b%iu2,b%ix:b%ix2)
  b%sparse%nc=0
  do j=b%iu,b%iu2
    do k=1,bm%nx_all
      if (k>=b%ix.and.k<=b%ix2) cycle
      if (m(j,k)/=d0) b%sparse%nc=b%sparse%nc+1
    enddo
  enddo
  allocate(b%sparse%m(b%sparse%nc))
  allocate(b%sparse%c(b%sparse%nc))
  allocate(b%sparse%r(b%nu+1))
  n=1
  do j=b%iu,b%iu2
    b%sparse%r(j-b%iu+1)=n
    do k=1,bm%nx_all
      if (k>=b%ix.and.k<=b%ix2) cycle
      if (m(j,k)/=d0) then
        b%sparse%m(n)=m(j,k)
        b%sparse%c(n)=k
        n=n+1
      endif
    enddo
  enddo
  b%sparse%r(b%nu+1)=n
enddo
end

subroutine bm_sparse_csr_csc(sp1,sp2,iu,iu2,m)
!преобразование sparse матрицы из csr в csc формат
use slau_block
type(sparse_matrix) sp1 !csr (input)
type(sparse_matrix) sp2 !csc (output)
integer(4) m !количество столбцов в общей матрице
integer(4) iu,iu2 !индексы первой и последней строки в sp1
integer(4), allocatable :: r(:)
integer(4) info,job(6)
allocate(r(m+1),sp2%c(m+1),sp2%r(sp1%nc))
sp2%nc=sp1%nc
if (iu>1) r(1:iu-1)=sp1%r(1)
r(iu:iu2)=sp1%r
r(iu2+1:m+1)=sp1%r(iu2-iu+2)
job(1)=0
job(2)=1
job(3)=1
job(6)=0
call mkl_dcsrcsc(job, m, sp1%m, sp1%c, r, null(), sp2%r, sp2%c, info)
deallocate(r)
end

function bm_get_u_ij(bmILU,i,j,u_ij) result(res)
!найти U_{ij}
!bm_get_u_ij=true - есть ненулевое значение U_{ij} 
use slau_block
logical res
type(main_block_matrix), target :: bmILU
type(main_block_matrix), pointer :: bm
type(block_matrix), pointer :: bU,bU_ILU
integer(4) i,j,i1,i2,i3,i4,j4
real(8) u_ij
res=.false.
bm=>bmILU%bm_ref
bU=>bm%blockm(bm%eq_block(i)) !содержит массивы c,r
bU_ILU=>bmILU%blockm(bm%eq_block(i)) !содержит массив m
i1=i-bU%iu+1 !локальный индекс строки для U
i2=-1
if (j<bU%ix) then
  i2=bU%sparse%r(i1)
  i3=bU%sparse%r2(i1)-1
elseif (j>bU%ix2) then
  i2=bU%sparse%r2(i1)
  i3=bU%sparse%r(i1+1)-1
endif
if (i2>0) then
  !U-sparse
  do i4=i2,i3
    if (bU%sparse%c(i4)==j) then
      u_ij=bU_ILU%sparse%m(i4)  
      res=.true.
      return
    elseif (bU%sparse%c(i4)>j) then
      return
    endif
  enddo
else
  !U-dens
  j4=j-bU%ix+1 !локальный индекс столбца для U
  u_ij=bU_ILU%m(i1,j4)  
  res=u_ij/=d0
endif
end

function bm_get_u_ij_P(bmILU,i,j,u_ij) result(res)
!найти U_{ij}
!bm_get_u_ij=true - есть ненулевое значение U_{ij} 
use slau_block
logical res
type(main_block_matrix), target :: bmILU
type(main_block_matrix), pointer :: bm
type(block_matrix), pointer :: bU,bU_ILU
integer(4) i,j,i1,i2,i3,i4,j4
real(8) u_ij
res=.false.
bm=>bmILU%bm_ref
bU=>bm%blockm(bm%eq_block(i)) !содержит массивы c,r
bU_ILU=>bmILU%blockm(bm%eq_block(i)) !содержит массив m
i1=i-bU%iu+1 !локальный индекс строки для U
i2=-1
if (j<bU%ix) then
  i2=bU%sparse%r(bU%p(i1))
  i3=bU%sparse%r2(bU%p(i1))-1
elseif (j>bU%ix2) then
  i2=bU%sparse%r2(bU%p(i1))
  i3=bU%sparse%r(bU%p(i1)+1)-1
endif
if (i2>0) then
  !U-sparse
  do i4=i2,i3
    if (bU%sparse%c(i4)==j) then
      u_ij=bU_ILU%sparse%m(i4)  
      res=.true.
      return
    elseif (bU%sparse%c(i4)>j) then
      return
    endif
  enddo
else
  !U-dens
  j4=j-bU%ix+1 !локальный индекс столбца для U
  u_ij=bU_ILU%m(bU%p(i1),j4)  
  res=u_ij/=d0
endif
end

subroutine bm_ILU_sum(bmILU,b,bILU,k1,j,imax,sum)
!подсчитать сумму -\sum_1^imax(l_ki*u_ij)
use slau_block
type(main_block_matrix), target :: bmILU
type(block_matrix), pointer :: b,bILU
integer(4) i,j,j1,j2,k1,i1,imax
real(8) sum,u_ij,l_ki
logical bm_get_u_ij
!k1 - локальный индекс уравнения в b,bILU
!L-sparse
j1=b%sparse%r(k1)
j2=b%sparse%r(k1+1)-1
!$omp parallel do reduction(- : sum) private (i1,i,l_ki,u_ij) shared (j1,j2,b,imax,bILU,j,bmILU)
do i1=j1,j2
  i=b%sparse%c(i1) !глобальный индекс столбца для L и строки для U
  !if (i>imax) exit
  if (i>imax) cycle
  l_ki=bILU%sparse%m(i1)
  if (bm_get_u_ij(bmILU,i,j,u_ij)) sum=sum-l_ki*u_ij
enddo
!omp end parallel do nowait
!L-dens
j1=min(imax,b%ix2)
!$omp parallel do reduction(- : sum) private (i,i1,l_ki,u_ij) shared (b,j1,bILU,k1,bmILU,j)
do i=b%ix,j1
  i1=i-b%ix+1 !локальный индекс столбца в bILU
  l_ki=bILU%m(k1,i1)
  if (l_ki/=d0) then
    if (bm_get_u_ij(bmILU,i,j,u_ij)) sum=sum-l_ki*u_ij
  endif
enddo
!omp end parallel do nowait
end

subroutine bm_ILUP_sum(bmILU,b,bILU,k1,j,imax,sum)
!подсчитать сумму -\sum_1^imax(l_ki*u_ij)
use slau_block
type(main_block_matrix), target :: bmILU
type(block_matrix), pointer :: b,bILU
integer(4) i,j,j1,j2,k1,i1,imax
real(8) sum,u_ij,l_ki
logical bm_get_u_ij_P
!k1 - локальный индекс уравнения в b,bILU
!L-sparse
j1=b%sparse%r(b%p(k1))
j2=b%sparse%r(b%p(k1)+1)-1
!$omp parallel do reduction(- : sum) private (i1,i,l_ki,u_ij) shared (j1,j2,b,imax,bILU,j,bmILU)
do i1=j1,j2
  i=b%sparse%c(i1) !глобальный индекс столбца для L и строки для U
  !if (i>imax) exit
  if (i>imax) cycle
  l_ki=bILU%sparse%m(i1)
  if (bm_get_u_ij_P(bmILU,i,j,u_ij)) sum=sum-l_ki*u_ij
enddo
!omp end parallel do nowait
!L-dens
j1=min(imax,b%ix2)
!$omp parallel do reduction(- : sum) private (i,i1,l_ki,u_ij) shared (b,j1,bILU,k1,bmILU,j)
do i=b%ix,j1
  i1=i-b%ix+1 !локальный индекс столбца в bILU
  l_ki=bILU%m(b%p(k1),i1)
  if (l_ki/=d0) then
    if (bm_get_u_ij_P(bmILU,i,j,u_ij)) sum=sum-l_ki*u_ij
  endif
enddo
!omp end parallel do nowait
end

subroutine bm_ILU(bmILU)
!построение ILU-факторицации для блочно-диагональной матрицы
use slau_block
type(main_block_matrix), pointer :: bm !исходная матрица
type(main_block_matrix), target :: bmILU !факторизованная матрица
type(block_matrix), pointer :: b,bILU
integer(4) j,k,j1,j2,k1,ia,j3,imax,k_1,i0
real(8) sum,u_jj
logical bm_get_u_ij
bm=>bmILU%bm_ref
if (bm%nu_all/=bm%nx_all) then
  write(*,*)"Error bm_ILU! (nu/=nx)"
  stop
endif
do ia=1,bm%nblock
  b=>bm%blockm(ia)
  bILU=>bmILU%blockm(ia)
  k1=0
  do k=b%iu,b%iu2
    k_1=k-1
    k1=k1+1 !локальный индекс уравнения
    do i0=1,2
      !LU-sparce
      if (i0==1) then
        j1=b%sparse%r(k1)
        j2=b%sparse%r2(k1)-1
      else
        j1=b%sparse%r2(k1)
        j2=b%sparse%r(k1+1)-1
      endif
      do j3=j1,j2 
        j=b%sparse%c(j3) !глобальный индекс столбца
        sum=bILU%sparse%m(j3)
        imax=min(j-1,k_1)
        if (imax>0) call bm_ILU_sum(bmILU,b,bILU,k1,j,imax,sum)
        if (j<=k_1) then !L делим на u_jj
          if (bm_get_u_ij(bmILU,j,j,u_jj)) then
            sum=sum/u_jj
          else
            write(*,*) "Error bm_ILU! u_jj=null (LU-sparce)"
            stop
          endif
        endif
        bILU%sparse%m(j3)=sum
      enddo
      if (i0==2) exit
      !LU-dens
      do j=b%ix,b%ix2
        j1=j-b%ix+1 !локальный индекс столбца
        sum=bILU%m(k1,j1)
        !if (full_LU_inBlock.or.sum/=d0) then !делать полную LU-факторизацию для блоков (учитывая и нулевые элементы)
          imax=min(j-1,k_1)
          if (imax>0) call bm_ILU_sum(bmILU,b,bILU,k1,j,imax,sum)
          if (j<=k_1) then !L делим на u_jj
            if (bm_get_u_ij(bmILU,j,j,u_jj)) then
              sum=sum/u_jj
            else
              write(*,*) "Error bm_ILU! u_jj=null (LU-dens)"
              stop
            endif
          endif
        !endif
        bILU%m(k1,j1)=sum
      enddo
    enddo
  enddo
enddo
end

subroutine bm_ILUP(bmILU)
!построение ILU-факторизации для блочно-диагональной матрицы
use slau_block
type(main_block_matrix), pointer :: bm !исходная матрица
type(main_block_matrix), target :: bmILU !факторизованная матрица
type(block_matrix), pointer :: b,bILU
integer(4) j,k,j1,j2,k1,ia,j3,imax,k_1,i0,jmax
real(8) sum,u_jj,amax,t
logical bm_get_u_ij_P
bm=>bmILU%bm_ref
if (bm%nu_all/=bm%nx_all) then
  write(*,*)"Error bm_ILU! (nu/=nx)"
  stop
endif
do ia=1,bm%nblock
  b=>bm%blockm(ia)
  bILU=>bmILU%blockm(ia)
  k1=0
  do k=b%iu,b%iu2
    k_1=k-1
    k1=k1+1 !локальный индекс уравнения
    !ищем главный элемент
    amax=dabs(bILU%m(b%p(k1),k1))
    jmax=k1
    do j1=k1+1,b%nu
      t=bILU%m(b%p(j1),k1)
      if (dabs(t)>amax) then
        jmax=j1
        amax=dabs(t)
      endif
    enddo
    if (jmax/=k1) call swap_int(b%p(k1),b%p(jmax))
    do i0=1,2
      !LU-sparce
      if (i0==1) then
        j1=b%sparse%r(b%p(k1))
        j2=b%sparse%r2(b%p(k1))-1
      else
        j1=b%sparse%r2(b%p(k1))
        j2=b%sparse%r(b%p(k1)+1)-1
      endif
      do j3=j1,j2 
        j=b%sparse%c(j3) !глобальный индекс столбца
        sum=bILU%sparse%m(j3)
        imax=min(j-1,k_1)
        if (imax>0) call bm_ILUP_sum(bmILU,b,bILU,k1,j,imax,sum)
        if (j<=k_1) then !L делим на u_jj
          if (bm_get_u_ij_P(bmILU,j,j,u_jj)) then
            sum=sum/u_jj
          else
            write(*,*) "Error bm_ILU! u_jj=null (LU-sparce)"
            stop
          endif
        endif
        bILU%sparse%m(j3)=sum
      enddo
      if (i0==2) exit
      !LU-dens
      do j=b%ix,b%ix2
        j1=j-b%ix+1 !локальный индекс столбца
        sum=bILU%m(b%p(k1),j1)
        !if (full_LU_inBlock.or.sum/=d0) then !делать полную LU-факторизацию для блоков (учитывая и нулевые элементы)
          imax=min(j-1,k_1)
          if (imax>0) call bm_ILUP_sum(bmILU,b,bILU,k1,j,imax,sum)
          if (j<=k_1) then !L делим на u_jj
            if (bm_get_u_ij_P(bmILU,j,j,u_jj)) then
              sum=sum/u_jj
            else
              write(*,*) "Error bm_ILU! u_jj=null (LU-dens)"
              stop
            endif
          endif
        !endif
        bILU%m(b%p(k1),j1)=sum
      enddo
    enddo
  enddo
enddo
end

subroutine bm_init_ILU(bm,bmILU,init_p)
use slau_block
type(main_block_matrix), target :: bm !исходная матрица
type(main_block_matrix), target :: bmILU !факторизованная матрица
type(block_matrix), pointer :: b,bILU
integer(4) i
logical init_p
call bm_init_block_r2(bm)
if (init_p) call bm_init_block_p(bm)
call bm_allocate_blocks(bmILU,bm%nblock)
do i=1,bm%nblock
  b=>bm%blockm(i)
  bILU=>bmILU%blockm(i)
  call bm_reallocate_block(bILU,b%nu,b%nx)
  allocate(bILU%sparse%m(b%sparse%nc))
  bILU%m=b%m
  bILU%sparse%m=b%sparse%m
enddo
call bm_init_block_matrix(bmILU)
bmILU%bm_ref=>bm
end

subroutine bm_ILU_solve(bm,bmILU,x)
!вычисление x=M^-1*b, где M - LU матрица (полученная, например, при ILU факторизации)
!!!для случая, когда диагональные элементы принадлежат блокам!
!Ly=b, Ux=y
use slau_block
type(main_block_matrix), target :: bm !исходная матрица, содержащая портрет (массивы sparse%r и sparse%c в блоках)
type(main_block_matrix), target :: bmILU !факторизованная матрица
real(8) x(:) !(bm%nu_all) на входе правая часть, на выходе - решение
real(8) v,ddot,ddoti,u_ii
type(block_matrix), pointer :: b,bILU
integer(4) ia,i,i1,j1,j2
!прямой проход по матрице L
do ia=1,bm%nblock
  b=>bm%blockm(ia)
  bILU=>bmILU%blockm(ia)
  do i=b%iu,b%iu2 !строка
    i1=i-b%iu+1   !локальный индекс строки
    !L-dens
    if (b%ix<=i-1) then !!!предполагаем, что диагональный элемент в блоке
      v=ddot(i-b%ix,bILU%m(i1,:),1,x(b%ix),1)
      x(i)=x(i)-v
    endif
    !L-sparse 
    j1=b%sparse%r(i1)
    j2=b%sparse%r2(i1) !!!предполагаем, что диагональный элемент в блоке
    if (j2>j1) then
      v=ddoti(j2-j1,bILU%sparse%m(j1),b%sparse%c(j1),x)
      x(i)=x(i)-v
    endif
  enddo
enddo
!обратный проход по матрице U
do ia=bm%nblock,1,-1
  b=>bm%blockm(ia)
  bILU=>bmILU%blockm(ia)
  do i=b%iu2,b%iu,-1 !строка
    i1=i-b%iu+1   !локальный индекс строки
    j1=i-b%ix+1   !локальный индекс столбца (диагонального элемента)
    u_ii=bILU%m(i1,j1)   !!!предполагаем, что диагональный элемент в блоке
    !U-dens
    if (b%ix2>=i+1) then !!!предполагаем, что диагональный элемент в блоке
      v=ddot(b%ix2-i,bILU%m(i1,j1+1:b%nx),1,x(i+1),1)
      x(i)=x(i)-v
    endif
    !U-sparse 
    j1=b%sparse%r2(i1)
    j2=b%sparse%r(i1+1) !!!предполагаем, что диагональный элемент в блоке
    if (j2>j1) then
      v=ddoti(j2-j1,bILU%sparse%m(j1),b%sparse%c(j1),x)
      x(i)=x(i)-v
    endif
    x(i)=x(i)/u_ii
  enddo
enddo
end

subroutine bm_swap_vector_P(bm,x)
use slau_block
!переставить элементы вектора в соответствии с с массивами b%p
type(main_block_matrix), target :: bm !исходная матрица, содержащая портрет (массивы sparse%r и sparse%c в блоках)
type(block_matrix), pointer :: b
real(8) x(:) !(bm%nu_all)
real(8), allocatable :: x1(:)
integer(4) ia,i,i1,j
allocate(x1(bm%nu_all))
x1=x
do ia=1,bm%nblock
  b=>bm%blockm(ia)
  do i=b%iu,b%iu2 !строка
    i1=i-b%iu+1   !локальный индекс строки
    j=b%p(i1)+b%iu-1
    x(i)=x1(j)
  enddo
enddo
deallocate (x1)
end

subroutine bm_ILUP_solve(bmILU,x)
!вычисление x=M^-1*b, где M - LU матрица (полученная, например, при ILU факторизации)
!!!для случая, когда диагональные элементы принадлежат блокам!
!Ly=b, Ux=y
use slau_block
type(main_block_matrix), pointer :: bm !исходная матрица, содержащая портрет (массивы sparse%r и sparse%c в блоках)
type(main_block_matrix), target :: bmILU !факторизованная матрица
real(8) x(:) !(bm%nu_all) на входе правая часть, на выходе - решение
real(8) v,ddot,ddoti,u_ii
type(block_matrix), pointer :: b,bILU
integer(4) ia,i,i1,j1,j2
bm=>bmILU%bm_ref
!переставляем x
call bm_swap_vector_P(bm,x)
!прямой проход по матрице L
do ia=1,bm%nblock
  b=>bm%blockm(ia)
  bILU=>bmILU%blockm(ia)
  do i=b%iu,b%iu2 !строка
    i1=i-b%iu+1   !локальный индекс строки
    !L-dens
    if (b%ix<=i-1) then !!!предполагаем, что диагональный элемент в блоке
      v=ddot(i-b%ix,bILU%m(b%p(i1),:),1,x(b%ix),1)
      x(i)=x(i)-v
    endif
    !L-sparse 
    j1=b%sparse%r(b%p(i1))
    j2=b%sparse%r2(b%p(i1)) !!!предполагаем, что диагональный элемент в блоке
    if (j2>j1) then
      v=ddoti(j2-j1,bILU%sparse%m(j1),b%sparse%c(j1),x)
      x(i)=x(i)-v
    endif
  enddo
enddo
!обратный проход по матрице U
do ia=bm%nblock,1,-1
  b=>bm%blockm(ia)
  bILU=>bmILU%blockm(ia)
  do i=b%iu2,b%iu,-1 !строка
    i1=i-b%iu+1   !локальный индекс строки
    j1=i-b%ix+1   !локальный индекс столбца (диагонального элемента)
    u_ii=bILU%m(b%p(i1),j1)   !!!предполагаем, что диагональный элемент в блоке
    !U-dens
    if (b%ix2>=i+1) then !!!предполагаем, что диагональный элемент в блоке
      v=ddot(b%ix2-i,bILU%m(b%p(i1),j1+1:b%nx),1,x(i+1),1)
      x(i)=x(i)-v
    endif
    !U-sparse 
    j1=b%sparse%r2(b%p(i1))
    j2=b%sparse%r(b%p(i1)+1) !!!предполагаем, что диагональный элемент в блоке
    if (j2>j1) then
      v=ddoti(j2-j1,bILU%sparse%m(j1),b%sparse%c(j1),x)
      x(i)=x(i)-v
    endif
    x(i)=x(i)/u_ii
  enddo
enddo
end

subroutine bm_BiCGStab_test
use slau_block
!x+y=3
!x-y=1
!x=2, y=1
type(main_block_matrix), target :: bm
type(block_matrix), pointer :: b
real(8) m(2,2),bb(2),x(2)
integer(4) i
logical write_iter_empty
external write_iter_empty
call bm_allocate_blocks(bm,2)
do i=1,2
  b=>bm%blockm(i)
  call bm_reallocate_block(b,1,1)
enddo
call bm_init_block_matrix(bm)
m(1,1)=d1
m(1,2)=d1
m(2,1)=d1
m(2,2)=-d1
bb=[3.0d0,d1]
x=d0
call bm_m_to_block_diagonal(m,bm)
call bm_BiCGStab(bm,bb,x,1000,1.0d-12,20,write_iter_empty)
end

subroutine bm_BiCGStab_test2
use slau_block
!x+z=4
!x+y=3
!y+z=5
!x=1, y=2, z=3
type(main_block_matrix), target :: bm
real(8) m(3,3),bb(3),x(3)
logical write_iter_empty
external write_iter_empty
call bm_allocate_blocks(bm,2)
call bm_reallocate_block(bm%blockm(1),2,1)
call bm_reallocate_block(bm%blockm(2),1,2)
call bm_init_block_matrix(bm)
m(1,:)=[d1,d0,d1]
m(2,:)=[d1,d1,d0]
m(3,:)=[d0,d1,d1]
bb=[4.0d0,3.0d0,5.0d0]
x=d0
call bm_m_to_block_diagonal(m,bm)
call bm_BiCGStab(bm,bb,x,1000,1.0d-12,20,write_iter_empty)
end

subroutine bm_ILU_test
use slau_block
real(8) m(7,7),x(7),bb(7),buff(7)
type(main_block_matrix) bm,bmILU
integer(4) precond
logical write_iter_empty
external write_iter_empty
m(2,:)=[9,0,0,3,1,0,1]
m(1,:)=[0,11,2,1,0,0,2]
m(4,:)=[0,1,10,2,0,0,0]
m(3,:)=[2,1,2,9,1,0,0]
m(7,:)=[1,0,0,1,12,0,1]
m(6,:)=[0,0,0,0,0,8,0]
m(5,:)=[2,2,0,0,3,0,8]
!-------
call bm_allocate_blocks(bm,3)
call bm_reallocate_block(bm%blockm(1),2,2)
call bm_reallocate_block(bm%blockm(2),2,2)
call bm_reallocate_block(bm%blockm(3),3,3)
!-------
!call bm_allocate_blocks(bm,1)
!call bm_reallocate_block(bm%blockm(1),7,7)
!--------
!call bm_allocate_blocks(bm,2)
!call bm_reallocate_block(bm%blockm(1),2,4)
!call bm_reallocate_block(bm%blockm(2),5,3)
!-------
call bm_init_block_matrix(bm)
call bm_m_to_block_diagonal(m,bm)
x=d1
x(1:3)=d2
call bm_A_mult_x(d1,bm,x,bb,buff) 
precond=0
select case (precond)
case (0)
  x=d0
case (1)
  !приближение ILU
  call bm_init_ILU(bm,bmILU,.true.)
  call bm_ILU(bmILU)
  x=bb
  call bm_ILU_solve(bm,bmILU,x)
case (2)
  !приближение ILUP
  call bm_init_ILU(bm,bmILU,.true.)
  call bm_ILUP(bmILU)
  x=bb
  call bm_ILUP_solve(bmILU,x)
end select
!тестируем решение СЛАУ
!call bm_calc_slau_block_diagonal(bm,bb,x,d1,1000,1.0d-12,20)
call bm_BiCGStab(bm,bb,x,1000,1.0d-12,20,write_iter_empty)
end

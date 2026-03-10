!!!дописать путь к dll
!!!в свойствах проекта Configuration Properties -> Debugging -> Environment
!!!PATH=D:\Users\Pehat\fort\common_fort;$(PATH);

!!!для подключения MKL: Configuration Properties -> Fortran -> Libraries -> Use Intel Math Kernel Library = Parallel
  
INCLUDE 'mkl_dss.f90'  
INCLUDE 'mkl_pardiso.f90'
!INCLUDE 'mkl_spblas.f90'    
!INCLUDE 'mkl_sparse_qr.f90'  

module slau_mod
use gen_mod
use kernel32
!use amgcl


interface
  subroutine VCL_LUFACT(a, b, n, nmax, x, toFloat)
  !dec$ attributes dllimport :: VCL_LUFACT
    integer(4) :: n
    integer(4) :: nmax
    real(8)    :: a(nmax,nmax), b(nmax), x(nmax)
    integer(4) :: toFloat
  end subroutine VCL_LUFACT
end interface

INTEGER (HANDLE) p, stat
pointer (pVCL_LUFact, VCL_LUFACT)
logical vcl_loaded

end module slau_mod

  
  
subroutine calc_slau(a, b, n, nmax, x, mode, find_norm, norm)
use slau_mod
integer(4) :: mode !0 - IMSL
                   !1 - MKL  
                   !2 - vcl OpenMP 
                   !3 - vcl OpenCL 
                   !4 - vcl OpenCL(toFloat)
integer(4) :: n
integer(4) :: nmax
integer(4) :: info,ipiv(n)
real(8)    b(nmax), x(nmax)
real(8) :: a(nmax,nmax)
real(8),allocatable    :: normv(:)
real(8),allocatable    :: a1(:,:)
logical find_norm   !найти норму ||Ax-b||
real(8) norm,dnrm2
if (find_norm) then
  allocate(a1(nmax,nmax))
  a1=a
endif
if (mode==0) then
  call DLSLRG(n,a,nmax,b,1,x) 
  !call DLSARG(n,a,nmax,b,1,x) !повышенная точность (на порядок по сравнению с DLSLRG), в 2-3 раза дольше работает
elseif (mode==1) then
  call dgetrf( n, n, a, nmax, ipiv, info )
  x=b
  call dgetrs( 'N', n, 1, a, nmax, ipiv, x, nmax, info )
elseif (mode==2) then
  call slau_vcl(a, b, n, nmax, x, 0, 0)
elseif (mode==3) then
  call slau_vcl(a, b, n, nmax, x, 0, 1)
elseif (mode==4) then
  call slau_vcl(a, b, n, nmax, x, 1, 1)
endif
!if (test_result) call test_slau_res(a,b,n,nmax,x)
if (find_norm) then
  allocate(normv(n))
  normv=b(1:n)
  call dgemv('N',n,n,d1,a1,nmax,x,1,-d1,normv,1)
  norm=dnrm2(n, normv, 1)
  deallocate(normv)
endif
if (allocated(a1)) deallocate(a1)
end
  
subroutine calc_slau_sparse(values, columns, rowIndex, rhs, nRows, nCols, nNonZeros , solution, mode, find_norm, norm, pardiso_eps)
!решение СЛАУ с разреженной матрицей
use gen_mod
integer(4) nRows   !число уравнений
integer(4) nCols   !число неизвестных
integer(4) nNonZeros   !число ненулевых коэффициентов
real(8) values(nNonZeros)  !коэффициенты
integer(4) columns(nNonZeros)  !индексы столбцов
integer(4) rowIndex(nRows+1) !индексы по строкам
real(8) rhs(nRows)   !столбец свободных членов
real(8) solution(nCols)   !решение
integer(4) mode !0 - MKL DSS
                !1 - MKL PARDISO
                !2 - amgcl
integer(4) pardiso_eps !погрешность для итераци процесса в  MKL PARDISO: =10:(-pardiso_eps)
integer(4) i
logical find_norm   !найти норму ||Ax-b||
real(8) norm,dnrm2
real(8), allocatable :: normv(:)
select case(mode)
case(0)
  call calc_slau_dss(values, columns, rowIndex, rhs, nRows, nCols, nNonZeros , solution)
case(1)
  call calc_slau_pardiso(values, columns, rowIndex, rhs, nRows, nCols, nNonZeros , solution, pardiso_eps)
case(2)
  !call amgcl_solve_slau(nRows, rowIndex, columns, values, rhs, solution, 1.0d-6, .true., norm)
  write(*,*) "AMGCL is comment!!!"
  call system("pause")
  stop
end select
if (find_norm.and.mode<2) then
  allocate(normv(nRows))
  call mkl_dcsrgemv('N',nRows,values,rowIndex,columns,solution,normv)
  forall (i=1:nRows) normv(i)=normv(i)-rhs(i)
  norm=dnrm2(nRows, normv, 1)
  deallocate(normv)
endif
end

subroutine calc_slau_dss(values, columns, rowIndex, rhs, nRows, nCols, nNonZeros, solution)
!решение СЛАУ с разреженной матрицей
!используя MKL DSS
!описание процедуры
!https://software.intel.com/en-us/mkl-developer-reference-c-direct-sparse-solver-dss-interface-routines#7C5F5C27-0E6D-4612-9460-E0F03A20DA83
!описание структуры данных
!https://software.intel.com/en-us/mkl-developer-reference-c-dss-nonsymmetric-matrix-storage#042CD407-36A0-4AE9-AC9E-34E8E0A004D3
!пример программы
!http://sep.stanford.edu/sep/claudio/Research/Prst_ExpRefl/ShtPSPI/intel/mkl/10.0.3.020/examples/solver/source/dss_sym_f90.f90
use gen_mod
use mkl_dss
integer(4) nRows   !число уравнений
integer(4) nCols   !число неизвестных
integer(4) nNonZeros   !число ненулевых коэффициентов
real(8) values(nNonZeros)  !коэффициенты
integer(4) columns(nNonZeros)  !индексы столбцов
integer(4) rowIndex(nRows+1) !индексы по строкам
real(8) rhs(nRows)   !столбец свободных членов
real(8) solution(nCols)   !решение
integer(4) err,perm(1),nRhs
TYPE(MKL_DSS_HANDLE) :: handle_dss ! Allocate storage for the solver handle.
nRhs=1
do while(.true.)
  ! Initialize the solver.
  err = DSS_CREATE( handle_dss, MKL_DSS_DEFAULTS )
  IF (err /= MKL_DSS_SUCCESS) exit
  ! Define the non-zero structure of the matrix.
  err = DSS_DEFINE_STRUCTURE( handle_dss, MKL_DSS_NON_SYMMETRIC, rowIndex,  nRows, nCols, columns, nNonZeros )
  IF (err /= MKL_DSS_SUCCESS) exit
  ! Reorder the matrix.
  err = DSS_REORDER( handle_dss, MKL_DSS_METIS_OPENMP_ORDER, perm )
  !err = DSS_REORDER( handle_dss, MKL_DSS_AUTO_ORDER, perm )
  IF (err /= MKL_DSS_SUCCESS) exit
  ! Factor the matrix.
  err = DSS_FACTOR_REAL( handle_dss, MKL_DSS_DEFAULTS, values )
  !err = DSS_FACTOR_REAL( handle_dss, MKL_DSS_INDEFINITE, values )
  IF (err /= MKL_DSS_SUCCESS) exit
  ! Solve the problem.
  err = DSS_SOLVE_REAL(handle_dss, MKL_DSS_DEFAULTS, rhs, nRhs, solution )
  IF (err /= MKL_DSS_SUCCESS) exit
  ! Deallocate solver storage.
  err = DSS_DELETE( handle_dss, MKL_DSS_DEFAULTS )
  exit
enddo
IF (err /= MKL_DSS_SUCCESS) then
  WRITE(*,*) "Error calc_slau_dss!"
  WRITE(*,*) "Solver returned error code ", err
  call system("pause")
  stop
endif
end

subroutine calc_slau_pardiso(values, columns, rowIndex, rhs, nRows, nCols, nNonZeros , solution, pardiso_eps)
!решение СЛАУ с разреженной матрицей
!используя MKL DSS
!описание процедуры
!https://software.intel.com/en-us/mkl-developer-reference-c-direct-sparse-solver-dss-interface-routines#7C5F5C27-0E6D-4612-9460-E0F03A20DA83
!описание структуры данных
!https://software.intel.com/en-us/mkl-developer-reference-c-dss-nonsymmetric-matrix-storage#042CD407-36A0-4AE9-AC9E-34E8E0A004D3
!пример программы
!http://sep.stanford.edu/sep/claudio/Research/Prst_ExpRefl/ShtPSPI/intel/mkl/10.0.3.020/examples/solver/source/dss_sym_f90.f90
use gen_mod
USE mkl_pardiso
integer(4) nRows   !число уравнений
integer(4) nCols   !число неизвестных
integer(4) nNonZeros   !число ненулевых коэффициентов
real(8) values(nNonZeros)  !коэффициенты
integer(4) columns(nNonZeros)  !индексы столбцов
integer(4) rowIndex(nRows+1) !индексы по строкам
real(8) rhs(nRows)   !столбец свободных членов
real(8) solution(nCols)   !решение
integer(4) pardiso_eps !погрешность для итераци процесса в  MKL PARDISO: =10:(-pardiso_eps)
integer(4) err,perm(1),iparm(64),i
TYPE(MKL_PARDISO_HANDLE), ALLOCATABLE  :: pt(:)
allocate(pt(64))
do i=1,64
   pt(i)%DUMMY=0 
   iparm(i)=0
enddo

iparm(1) = 1 ! no solver default
iparm(2) = 3 ! The parallel (OpenMP) version of the nested dissection algorithm
iparm(3) = 0 ! reserved, =0
iparm(4) = 0 ! no iterative-direct algorithm
iparm(5) = 0 ! no user fill-in reducing permutation
iparm(6) = 0 ! =0 solution on the first n compoments of x

iparm(8) = 10 ! numbers of iterative refinement steps
iparm(10) = pardiso_eps !13 ! perturbe the pivot elements with 1E-13

iparm(11) = 1 ! use nonsymmetric permutation and scaling MPS
iparm(12) = 0 ! Ax=b
iparm(13) = 0 ! maximum weighted matching algorithm is switched-off (default for symmetric). Try iparm(13) = 1 in case of inappropriate accuracy
iparm(14) = 0 ! Output: number of perturbed pivots
iparm(18) = -1 ! Output: number of nonzeros in the factor LU
iparm(19) = -1 ! Output: Mflops for LU factorization
iparm(20) = 0 ! Output: Numbers of CG Iterations


do while(.true.)
  CALL pardiso (pt, 1, 1, 11, 11, nRows, values, rowIndex, columns, perm, 1, iparm, 0, rhs, solution, err)
  IF (err /= 0) exit
  CALL pardiso (pt, 1, 1, 11, 23, nRows, values, rowIndex, columns, perm, 1, iparm, 0, rhs, solution, err)
  IF (err /= 0) exit
  !CALL pardiso (pt, 1, 1, 11, 33, nRows, values, rowIndex, columns, perm, 1, iparm, 0, rhs, solution, err)
  !IF (err /= 0) exit
  CALL pardiso (pt, 1, 1, 11, -1, nRows, values, rowIndex, columns, perm, 1, iparm, 0, rhs, solution, err)
  exit
enddo

IF (err /= 0) then
  WRITE(*,*) "Error calc_slau_pardiso!"
  WRITE(*,*) "Solver returned error code ", err
  call system("pause")
  stop
endif
deallocate(pt)
end

subroutine load_vcl(mode)
use slau_mod
integer(4) :: mode !0 - OpenMP, 1 - OpenCL 
character(12) libname
if (vcl_loaded) return
if (mode==0) then
  !DEC$ IF DEFINED(_M_IX86)
  libname="vclo32.dll"c
  !DEC$ ELSE
  libname="vclo64.dll"c
  !DEC$ ENDIF
  p = loadlibrary(libname)
  if (p==0) then
    libname="vcl.dll"c
    p = loadlibrary(libname)
  endif
else
  !DEC$ IF DEFINED(_M_IX86)
  libname="vcl32.dll"c
  !DEC$ ELSE
  libname="vcl64.dll"c
  !DEC$ ENDIF
  p = loadlibrary(libname)
  if (p==0) then
    libname="vcl.dll"c
    p = loadlibrary(libname)
  endif
endif
if (p/=0) then
  !DEC$ IF DEFINED(_M_IX86)
  pVCL_LUFact = getprocaddress(p, "_VCL_LUFACT@24"c)
  !DEC$ ELSE
  pVCL_LUFact = getprocaddress(p, "VCL_LUFACT"c)
  !DEC$ ENDIF
  if (pVCL_LUFact==0) then
    write(*,"('Error occurred finding function ''VCL_LUFact'' in ', A12)") libname
    call system("pause")
    stop
  end if
else
  write (*,"('Error occurred opening ', A12)") libname
  call system("pause")
  stop
end if
vcl_loaded=.true.
end

subroutine calc_slau_ls(a, b, m, n, nmax, x, mode,find_norm,norm)
!для переопределенной матрицы (m*n) m>n методом наименьших квадратов
use slau_mod
integer(4) :: mode !0 - IMSL
                   !1 - MKL  
integer(4) :: n,m,nmax,k,info,lwork
real(8)    :: a(nmax,n), b(nmax), x(n), tol
real(8), allocatable :: work(:)
real(8), allocatable :: a1(:,:),res(:)
real(8) norm,dnrm2
logical find_norm
if (find_norm.or.mode==0) allocate(res(m))
if (mode==0) then
  tol=1.0d-8
  call DLSQRR(m,n,a,nmax,b,tol,x,res,k)
  !call DLSBRR(m,n,a,nmax,b,tol,x,res,k) !с итерационным уточнением
elseif (mode==1) then
  lwork=n*64
  allocate(work(lwork))
  if (find_norm) then
    allocate(a1(m,n))
    a1=a(1:m,1:n)
    res=b(1:m)
  endif
  call dgels('N', m, n, 1, a, nmax, b, nmax, work,lwork, info)
  deallocate(work)
  x(1:n)=b(1:n)
  if (find_norm) call dgemv('N', m, n, d1, a1, m, x, 1, -d1, res, 1)
endif
if (find_norm) norm=dnrm2(m, res, 1)
if (allocated(res)) deallocate(res)
if (allocated(a1)) deallocate(a1)
end

!subroutine calc_slau_sparse_ls(values, columns, rowIndex, rhs, nRows, nCols, nNonZeros , solution, find_norm, norm)
!!решение СЛАУ с разреженной матрицей
!use slau_mod
!USE MKL_SPBLAS
!USE MKL_SPARSE_QR
!integer(4) nRows   !число уравнений
!integer(4) nCols   !число неизвестных
!integer(4) nNonZeros   !число ненулевых коэффициентов
!real(8) values(nNonZeros)  !коэффициенты
!integer(4) columns(nNonZeros)  !индексы столбцов
!integer(4) rowIndex(nRows+1) !индексы по строкам
!real(8) rhs(nRows)   !столбец свободных членов
!real(8) solution(nCols)   !решение
!logical find_norm   !найти норму ||Ax-b||
!real(8) norm,dnrm2
!TYPE(SPARSE_MATRIX_T) csrA
!TYPE(MATRIX_DESCR) descrA
!INTEGER info
!descrA % TYPE = SPARSE_MATRIX_TYPE_GENERAL
!!Create CSR matrix
!info = MKL_SPARSE_D_CREATE_CSR(csrA,SPARSE_INDEX_BASE_ONE,nRows,nCols,rowIndex(1),rowIndex(2),columns,values)
!!Solve Ax=b using Sparse QR decomposition
!info = MKL_SPARSE_D_QR(SPARSE_OPERATION_NON_TRANSPOSE, csrA, descrA, SPARSE_LAYOUT_COLUMN_MAJOR, 1, solution, ncols, rhs, nrows )
!!Release internal representation of CSR matrix
!info = MKL_SPARSE_DESTROY(csrA)
!end

subroutine free_vcl
use slau_mod
if (vcl_loaded) stat = freelibrary(p)
end

subroutine slau_vcl(a, b, n, nmax, x, toFloat, mode)
use slau_mod
integer(4) :: mode !0 - OpenMP, 1 - OpenCL 
integer(4) :: n
integer(4) :: nmax
real(8)    :: a(nmax,nmax), b(nmax), x(nmax)
integer(4) :: toFloat
call load_vcl(mode)
if (pVCL_LUFact/=0) then
  call VCL_LUFact(a, b, n, nmax, x, toFloat)
endif
end

subroutine test_slau_res(a,b,n,nmax,x)
use gen_mod
integer(4) n,nmax
real(8) a(nmax,nmax),b(nmax),x(nmax),s(n),smax
integer(4) i,j
smax=d0
do i=1,n
  s(i)=d0
  do j=1,n
    s(i)=s(i)+a(i,j)*x(j)
  enddo
  s(i)=s(i)-b(i)
  if (dabs(s(i))>smax) smax=dabs(s(i))
enddo
end

subroutine test_qr
!x+y=3  x=2
!2x+y=5  y=1
!2x+2y=6
use gen_mod
integer(4), parameter :: n=2, m=3
real(8) a(m,n),b(m),x(n),a1(m,n),b1(m),norm
data a/d1,d2,d2,d1,d1,d2/
data b/3.0d0,5.0d0,6.0d0/
b(3)=5.9d0
a1=a
b1=b
!call calc_slau(a, b, n, m, x, 1)
call calc_slau_ls(a, b, m, n, m, x, 0,.true.,norm)
end
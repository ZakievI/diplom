function ffunc_inf(xx,yy,ia)
use pgmod
real(8) ffunc_inf,xx,yy
integer(4) ia
ffunc_inf=d0
return
xx=xx
yy=yy
ia=ia
end

subroutine pg_get_matrix
!dec$ attributes dllexport:: pg_get_matrix
real(8) ffunc_inf
external ffunc_inf
call pg_get_matrix_finf(ffunc_inf)
end

subroutine pg_get_matrix_finf(func_inf)
!dec$ attributes dllexport:: pg_get_matrix2
!собрать матрицу
use pgmod
integer(4) i,j,k
type(matrix_mb), pointer :: m
real(8) func_inf 
external func_inf
!if (gs_profiling) gs_use_parallel_matrix_build=0
!gs_use_parallel_get_fun=0
if (first_matrix_build.or.gs_rebuild_matrix_for_iterate) then
  call geom_collocate_points
  call null_matrix
  j=1
  k=1
  do i=1,gs%na
    m=>gs%a(i)%m
    m%ix=j
    m%iu=k !чтобы замыкающие уравнения добавлялись к нужному блоку (используется в init_cpp_j)
    call get_ind_area(i)
    m%nnt=m%nu   !если к блоку будут добавляться closing уравнения
    m%nu_all=m%nx+m%nu_add
    j=j+m%nx
    k=k+m%nu_all
    m%ix2=j-1
    m%iu2=k-1
  enddo
  j=j-1
  k=k-1
  !учитываем неучтенные gs%constvalinds
  do i=1,gs%constvaln
    if (gs%constvalinds(i)==0) then
      !!!все неизвестные должны быть приписаны к какому-либо блоку
      call gs_print_stop("Error gs%constvalinds! (pg_get_matrix_finf)") 
      !j=j+1
      !gs%constvalinds(i)=j
    endif
  enddo
  if (gs%constvalinds_fict_i>0) call get_ind_area3
  call geom_collocate_points_area
  !выделяем память под общую матрицу
  call reallocate_mm(k,j)
  if (gs%m%matrix_store_type==2) then
    call bm_init_block_matrix(gs%m%bm)
    call init_bm_order
  endif
  !второй проход для фиктивных переменных
  if (gs%m%nx_fict>0) then
    do i=1,gs%na
      call get_ind_area2(i,j)
    enddo
  endif
  call init_cpp_j
  call prepare_global_sparsem
  !инициализируем кэш
  if (gs_use_cash) call initCash
  call init_omp_buffer
  first_matrix_build=.false.
endif
!составляем матрицу
call gs_print("Build matrix")
call init_gs_time
do i=1,gs%na
  call get_matrix_area(i,func_inf)
  !if (gs_profiling.and.i==10) exit
enddo
call print_total_time
gs%m%have_initial_approxsolv=.false.
end

subroutine init_bm_order
use pgmod
integer(4) i,j,k,imax,level(gs%na),levela(gs%na)
type(TArea), pointer :: a
imax=0
do i=1,gs%na
  a=>gs%a(i)
  if (a%iorder>imax) imax=a%iorder
enddo
if (imax==0) then
  forall (i=1:gs%na) gs%m%bm%ord(i)=i
  !gs%m%bm%nlevel=1
  !gs%m%bm%klevel=1
else
  !gs%m%bm%nlevel=imax
  j=0
  do k=1,imax
    do i=1,gs%na
      a=>gs%a(i)
      if (a%iorder==k) then
        j=j+1
        gs%m%bm%ord(j)=i
        !gs%m%bm%klevel(j)=k
        level(j)=k
        levela(i)=k
      endif
    enddo
  enddo
  if (j/=gs%na) call gs_print_stop("Error pg_init_subdmain_block_order! (ord)")
endif
end

subroutine reallocate_mm(nu,nx)
!выделение памяти глобальной матрицы
use pgmod
type(matrix_main), pointer :: m
type(matrix_mb), pointer :: ma
type(block_matrix), pointer :: b
integer(4) nx,nu,i
logical need_mm
m=>gs%m
need_mm=(nx/=m%nx).or.(nu/=m%nu).or.(m%allocated_matrix_store_type/=m%matrix_store_type)
call deallocate_mm(need_mm)
if (need_mm) then
  select case (m%matrix_store_type)
  case (0)
    gs%m%memory_mb=8.0d0*nx*nu/1024.0d0**2
    write(*,"('neq=',i0,' MatrixSize=',F9.1,'Mb')") nu,gs%m%memory_mb
    allocate(m%m(nu,nx))
  case (1) 
    allocate(m%sparse%mrows(nu))
  end select
	allocate(m%b(nu))	
  m%nx=nx
  m%nu=nu
endif
select case (m%matrix_store_type)
case (0)
  m%m=d0 !это нужно, т.к. потом строки присваиваются только частично
case (1)
  do i=1,nu
    call bm_null_mrow(m%sparse%mrows(i))
  enddo
case (2)
  call bm_allocate_blocks(gs%m%bm,gs%na)
  do i=1,gs%na
    ma=>gs%a(i)%m
    b=>gs%m%bm%blockm(i)
    call bm_reallocate_block(b,ma%nx,ma%nx)
  enddo
end select
m%b=d0
if (m%nx_fict>0) allocate(m%fvar(m%nx_fict))
m%ix_fict=nx+1
m%nx_all=nx+m%nx_fict
m%allocated_matrix_store_type=m%matrix_store_type
end

subroutine null_matrix
use pgmod
integer(4) i,j,k,i2
type(TArea), pointer :: a
type(TBound), pointer :: b
type(TBoundline2), pointer :: bl2
gs%m%nx_fict=0
do i=1,gs%na
  a=>gs%a(i)
  do j=1,a%nb
    b=>a%bnd(j)
    select case (b%boundLineType)
    case (1)
      b%psiom=d0
      b%psiind=0
    case (2)
      do k=1,b%nline
        bl2=>b%line2(k)
        do i2=1,a%umax
          call set_ga2_base(bl2%ga(i2))
        enddo
      enddo
    end select
  enddo
  if (a%haveAreaEq) then
    a%a%areaind=0
    a%a%psiarea=d0
    !a%a%arval_funp=d0
  endif
enddo
if (gs%constvaln>0) then
  gs%constvalinds=0
  gs%constvals=d0
  gs%constvalinds_fict=0
  gs%constvalinds_fict_i=0
endif
end

subroutine initCash
use pgmod
integer(4) i
gs_cash_is_solve=.false.
if (gs_cash_presolve) call init_gs_time
do i=1,gs%na
  call init_need_cash_for_matrix(i)
  call initCash_area(i)
enddo
if (gs_cash_presolve.and.gs_cash_is_solve) call print_total_time
end

subroutine gs_printInit_cash
use pgmod
if (gs_cash_is_solve) return
call gs_print("Init cash")
gs_cash_is_solve=.true.
end

subroutine init_need_cash_eq1(mode,der,need_cash_b,need_cash_a)
!инициализация кэша: уравнение Лапласа
use pgmod
integer(4) mode,der
logical need_cash_b(0:max_int) !интегралы по границе
logical need_cash_a(0:max_int) !интегралы по области
select case (mode)
case (1)
  need_cash_b(0:1)=.true.
case (2)
  select case (der)
  case (1)
    need_cash_b(4:6:2)=.true.
  case (2)
    need_cash_b(5:7:2)=.true.
  case default
    call gs_print_stop("Error init_need_cash_eq1")
  end select
case default
  call gs_print_stop("Error init_need_cash_eq1")
end select
return
need_cash_a(0)=need_cash_a(0)
end

subroutine init_need_cash_eq3(mode,der,need_cash_b,need_cash_a)
!инициализация кэша: бигармоническое уравнение (для главной функции)
use pgmod
integer(4) mode,der
logical need_cash_b(0:max_int) !интегралы по границе
logical need_cash_a(0:max_int) !интегралы по области
select case (mode)
case (1)
  need_cash_b(0:3)=.true.
case (2)
  select case (der)
  case (1)
    need_cash_b(4:10:2)=.true.
  case (2)
    need_cash_b(5:11:2)=.true.
  case default
    call gs_print_stop("Error init_need_cash_eq3")
  end select
case default
  call gs_print_stop("Error init_need_cash_eq3")
end select
return
need_cash_a(0)=need_cash_a(0)
end

subroutine init_need_cash_eq4(mode,der,need_cash_b,need_cash_a)
!инициализация кэша: уравнение Пуассона
                    !eq5 - уравнение Гельмгольца с переменными коэф, однородное и неоднородное (через неизвесные в области) eq=5,6
                    !eq8 - система вида \Delta\psi=v(x,y)\omega  и уравнение на \omega c неизвестными в области
                    !eq16 - однородное уравнение переноса \Delta C = vx*dC/dx + vy*dC/dy
                    !eq17 - неоднородное уравнение переноса \Delta C = vx*dC/dx + vy*dC/dy + f
                    !eq26 - уравнение Пуассона с неизвестными в правой части
use pgmod
integer(4) mode,der
logical need_cash_b(0:max_int) !интегралы по границе
logical need_cash_a(0:max_int) !интегралы по области
select case (mode)
case (1)
  need_cash_b(0:1)=.true.
  need_cash_a(0)=.true.
case (2)
  select case (der)
  case (1)
    need_cash_b(4:6:2)=.true.
    need_cash_a(4)=.true.
  case (2)
    need_cash_b(5:7:2)=.true.
    need_cash_a(5)=.true.
  case default
    call gs_print_stop("Error init_need_cash_eq4")
  end select
case default
  call gs_print_stop("Error init_need_cash_eq4")
end select
end

subroutine init_need_cash_eq9(mode,der,need_cash_b,need_cash_a)
!инициализация кэша: однородное уравнение Гельмгольца через функции Бесселя \Delta\psi+k_helm^2\psi=0
use pgmod
integer(4) mode,der
logical need_cash_b(0:max_int) !интегралы по границе
logical need_cash_a(0:max_int) !интегралы по области
select case (mode)
case (1)
  need_cash_b(12:13)=.true.
case (2)
  call gs_print_stop("Error init_need_cash_eq9. Derivatives are not implemented!!!")
case default
  call gs_print_stop("Error init_need_cash_eq9")
end select
return
need_cash_a(0)=der==der
end

subroutine init_need_cash_eq10(mode,der,need_cash_b,need_cash_a)
!инициализация кэша: однородное уравнение Гельмгольца через функции Бесселя \Delta\psi-k_helm^2\psi=0
use pgmod
integer(4) mode,der
logical need_cash_b(0:max_int) !интегралы по границе
logical need_cash_a(0:max_int) !интегралы по области
select case (mode)
case (1)
  need_cash_b(14:15)=.true.
case (2)
  select case (der)
  case (1)
    need_cash_b(16:18:2)=.true.
  case (2)
    need_cash_b(17:19:2)=.true.
  case default
    call gs_print_stop("Error init_need_cash_eq10")
  end select
case default
  call gs_print_stop("Error init_need_cash_eq10")
end select
return
need_cash_a(0)=need_cash_a(0)
end

subroutine init_need_cash_eq15(mode,der,need_cash_b,need_cash_a)
!инициализация кэша: уравнение Бринкмана \Delta^2\psi-k_helm^2\Delta\psi=0 через единую функцию Грина, содержащую функции Бесселя
use pgmod
integer(4) mode,der
logical need_cash_b(0:max_int) !интегралы по границе
logical need_cash_a(0:max_int) !интегралы по области
select case (mode)
case (1)
  need_cash_b(0:1)=.true.
  need_cash_b(14:15)=.true.
case (2)
  select case (der)
  case (1)
    need_cash_b(4:6:2)=.true.
    need_cash_b(16:18:2)=.true.
  case (2)
    need_cash_b(5:7:2)=.true.
    need_cash_b(17:19:2)=.true.
  case default
    call gs_print_stop("Error init_need_cash_eq15")
  end select
case default
  call gs_print_stop("Error init_need_cash_eq15")
end select
return
need_cash_a(0)=need_cash_a(0)
end

subroutine init_need_cash_eq21(mode,der,need_cash_b,need_cash_a)
!инициализация кэша: неоднородное бигармоническое уравнение (для главной функции)
use pgmod
integer(4) mode,der
logical need_cash_b(0:max_int) !интегралы по границе
logical need_cash_a(0:max_int) !интегралы по области
select case (mode)
case (1)
  need_cash_b(0:3)=.true.
  need_cash_a(2)=.true.
case (2)
  select case (der)
  case (1)
    need_cash_b(4:10:2)=.true.
    need_cash_a(8)=.true.
  case (2)
    need_cash_b(5:11:2)=.true.
    need_cash_a(9)=.true.
  case default
    call gs_print_stop("Error init_need_cash_eq21")
  end select
case default
  call gs_print_stop("Error init_need_cash_eq21")
end select
end

subroutine init_need_cash(mode,der,is_ba,i,init_need_cash_1)
!инициализация кэша: выбор расположения контрольных точек
!заменяет init_need_cash2 для уравнений второго порядка
use pgmod
integer(4) mode,der,is_ba,i
external init_need_cash_1
type(TArea), pointer :: a
a=>gs%a(i)
select case (is_ba)
case (1) !контрольные точки на границе
  call init_need_cash_1(mode,der,a%need_cash_bb,a%need_cash_ab)
case (2) !контрольные точки в области
  call init_need_cash_1(mode,der,a%need_cash_ba,a%need_cash_aa)
case default
  call gs_print_stop("Error init_need_cash")
end select
end

subroutine init_need_cash2(mode,der,is_ba,i,init_need_cash_1,init_need_cash_2)
!инициализация кэша: выбор главной или второй функции
use pgmod
integer(4) mode,der,is_ba,i
external init_need_cash_1,init_need_cash_2
select case (mode)
case (1,2) !главная функция
  call init_need_cash(mode,der,is_ba,i,init_need_cash_1)
case (3,4) !вторая функция
  call init_need_cash(mode-2,der,is_ba,i,init_need_cash_2)
case default
  call gs_print_stop("Error init_need_cash2")
end select
end

subroutine pg_init_need_cash_for_func(mode,der,is_ba,clear_need_cash)
!dec$ attributes dllexport:: pg_init_need_cash_for_func
!инициализаия кэша для вычисления функций для текущей области
use pgmod
integer(4) mode !1 - для вычисления главной функции
                !2 - для вычисления производной главной функции
                !3 - для вычисления второй функции
                !4 - для вычисления производной второй функции
integer(4) der  !тип производной
                !1 - d/dx
                !2 - d/dy
integer(4) is_ba  !1 - на границе
                  !2 - в области
logical clear_need_cash !очистить массивы. true для первого вызова
call init_need_cash_for_func(mode,der,is_ba,gsarea%i,clear_need_cash)
end

subroutine init_need_cash_for_func(mode,der,is_ba,ia,clear_need_cash)
!инициализаия кэша для вычисления функций
use pgmod
integer(4) mode !1 - для вычисления главной функции
                !2 - для вычисления производной главной функции
                !3 - для вычисления второй функции
                !4 - для вычисления производной второй функции
integer(4) der  !тип производной
                !1 - d/dx
                !2 - d/dy
integer(4) is_ba  !1 - на границе
                  !2 - в области
integer(4) ia     !индекс области
logical clear_need_cash !очистить массивы. true для первого вызова
integer(4) type_eq,type_eq2
type(TArea), pointer :: a
external init_need_cash_eq1,init_need_cash_eq3,init_need_cash_eq4,init_need_cash_eq9,init_need_cash_eq10,init_need_cash_eq15,init_need_cash_eq21
a=>gs%a(ia)
type_eq=a%type_eq(1)
type_eq2=0
if (a%nu==2) type_eq2=a%type_eq(2)
if (clear_need_cash) then
  a%need_cash_bb=.false.
  a%need_cash_ab=.false.
  a%need_cash_aa=.false.
  a%need_cash_ba=.false.
endif
select case (type_eq)
case (1)
  call init_need_cash(mode,der,is_ba,ia,init_need_cash_eq1)
case (3)
  call init_need_cash2(mode,der,is_ba,ia,init_need_cash_eq3,init_need_cash_eq1)
case (4,5,6,16:20,26)
  call init_need_cash(mode,der,is_ba,ia,init_need_cash_eq4)
case (8,11)
  select case (type_eq2)
  case (2)
    call init_need_cash2(mode,der,is_ba,ia,init_need_cash_eq4,init_need_cash_eq1)
  case (7,14)
    call init_need_cash2(mode,der,is_ba,ia,init_need_cash_eq4,init_need_cash_eq4)
  case (12)
    call init_need_cash2(mode,der,is_ba,ia,init_need_cash_eq4,init_need_cash_eq9)
  case (13)
    call init_need_cash2(mode,der,is_ba,ia,init_need_cash_eq4,init_need_cash_eq10)
  endselect
case (9)
  call init_need_cash(mode,der,is_ba,ia,init_need_cash_eq9)
case (10)
  call init_need_cash(mode,der,is_ba,ia,init_need_cash_eq10)
case (15)
  call init_need_cash2(mode,der,is_ba,ia,init_need_cash_eq15,init_need_cash_eq10)
case (21,22,24,25)
  call init_need_cash2(mode,der,is_ba,ia,init_need_cash_eq21,init_need_cash_eq4)
endselect
end

subroutine init_need_cash_for_matrix(i)
!инициализаия кэша для составления матрицы
use pgmod
integer(4) i
integer(4) type_eq,type_eq2
type(TArea), pointer :: a
a=>gs%a(i)
type_eq=a%type_eq(1)
type_eq2=0
if (a%nu==2) type_eq2=a%type_eq(2)
call init_need_cash_for_func(1,0,1,i,.true.) !общий для всех
select case (type_eq)
case (1,4,9,10)
  !вынесли вверх как общий для всех
case (3,15,21)
  call init_need_cash_for_func(3,0,1,i,.false.)
case (5,6)
  call init_need_cash_for_func(1,0,2,i,.false.)
case (8,11)
  call init_need_cash_for_func(3,0,1,i,.false.)
  call init_need_cash_for_func(3,0,2,i,.false.)
case (16,17)
  call init_need_cash_for_func(2,1,2,i,.false.)
  call init_need_cash_for_func(2,2,2,i,.false.)
case (18:20)
  call init_need_cash_for_func(2,2,2,i,.false.)
case (22,24,25)
  call init_need_cash_for_func(3,0,1,i,.false.)
  if (a%eq_var(2)) call init_need_cash_for_func(1,0,2,i,.false.) !psi
  if (a%eq_var(3)) call init_need_cash_for_func(3,0,2,i,.false.) !eta
  if (a%eq_var(4)) call init_need_cash_for_func(2,1,2,i,.false.) !dpsi/dx
  if (a%eq_var(5)) call init_need_cash_for_func(2,2,2,i,.false.) !dpsi/dy
  if (a%eq_var(6)) call init_need_cash_for_func(4,1,2,i,.false.) !deta/dx
  if (a%eq_var(7)) call init_need_cash_for_func(4,2,2,i,.false.) !deta/dy
case (26)
  if (a%eq_var(2)) call init_need_cash_for_func(1,0,2,i,.false.) !psi
  if (a%eq_var(3)) call init_need_cash_for_func(2,1,2,i,.false.) !dpsi/dx
  if (a%eq_var(4)) call init_need_cash_for_func(2,2,2,i,.false.) !dpsi/dy
endselect
end

subroutine test_need_cash(need_cash_bb_,need_cash_ba_,need_cash_ab_,need_cash_aa_,need_cash_bb,need_cash_ba,need_cash_ab,need_cash_aa)
use pgmod
logical need_cash_bb(0:max_int),need_cash_ab(0:max_int),need_cash_aa(0:max_int),need_cash_ba(0:max_int)
logical need_cash_bb_(0:max_int),need_cash_ab_(0:max_int),need_cash_aa_(0:max_int),need_cash_ba_(0:max_int)
integer(4) i
do i=0,max_int
  if (need_cash_bb_(i)/=need_cash_bb(i)) call gs_print_stop('test_need_cash: bb '//itoa(i))
  if (need_cash_ba_(i)/=need_cash_ba(i)) call gs_print_stop('test_need_cash: ba '//itoa(i))
  if (need_cash_ab_(i)/=need_cash_ab(i)) call gs_print_stop('test_need_cash: ab '//itoa(i))
  if (need_cash_aa_(i)/=need_cash_aa(i)) call gs_print_stop('test_need_cash: aa '//itoa(i))
enddo
end

subroutine pg_initCash_area_ar(ia1,ia2)
!dec$ attributes dllexport:: pg_initCash_area_ar
!инициализировать кэш интегралов для областей в диапазоне
use pgmod
integer(4) ia1,ia2 !диапазон индексов областей
integer(4) i
gs_cash_is_solve=.false.
if (gs_cash_presolve) call init_gs_time
do i=ia1,ia2
  call initCash_area(i)
enddo
if (gs_cash_presolve.and.gs_cash_is_solve) call print_total_time
end

subroutine pg_initCash_area
!dec$ attributes dllexport:: pg_initCash_area
!инициализировать кэш интегралов для текущей области
use pgmod
call pg_initCash_area_ar(gsarea%i,gsarea%i)
end

subroutine init_type_solve(ia,type_rotate,ni,nf,type_solve,nf_ref,ci_ref)
use pgmod
integer(4) ia,i
integer(4) type_rotate !тип перемещения области
                       !-1 - нет ссылки
                       !0 - просто перенос
                       !1 - поворот на 180 град
                       !2 - поворот на 90
                       !3 - поворот на 270
                       !4 - произвольный поворот
integer(4) nf(max_int+1),ni,nf1,nf2,nf_ref(max_int+1)
integer(4) type_solve(max_int+1) !способ вычисления интеграла
                                 !0 - вычислить стандартным способом
                                 !>0 - ускоренное вычисление
                                 !1 - взять со знаком "+" интеграл nf_ref(i)
                                 !2 - взять со знаком "-" интеграл nf_ref(i)
                                 !3 - вычислить по формуле поворота для d/dx
                                      !nf(i)*cosa-nf_ref(i)*sina   nf_ref(i) - содержит номер интеграла для d/dy
                                 !4 - вычислить по формуле поворота для d/dy
                                      !nf(i)*cosa+nf_ref(i)*sina   nf_ref(i) - содержит номер интеграла для d/dx
type(TCashIntegral) ci_ref
type(TArea), pointer :: a
do i=1,ni
  type_solve(i)=0
  nf_ref(i)=-1
enddo

!return  !чтобы отказаться от алгоритма ускоренного вычисления кэша, раскомментарить эту строку
         !в этом случае вычисление будет происходить по обычному через интегралы

a=>gs%a(ia)
if (type_rotate<0.or.a%oss_resolve) return
!!!ci_ref передается = null при type_rotate<1
!!!при компиляции в режиме Release оптимизация компилятора почемуто строит код так, что происходит обращение к ci_ref%solved
!!!там, где этого не должно происходить - это вызывает access violation
!!!хотя в режиме debug все нормально, до обращение к ci_ref%solved не доходит
!!!поэтому код для type_rotate<0 вынесли из цикла вперед отдельно
do i=1,ni
  nf1=nf(i)
  if (type_rotate==0.or.int_type(nf1)==0) then
    if (ci_ref%solved(nf1)) then
      type_solve(i)=1
      nf_ref(i)=nf1
    endif
  else
    if (int_type(nf1)==1) then
      nf2=nf1+1
    else
      nf2=nf1-1
    endif
    select case (type_rotate)
    case (1)
      if (ci_ref%solved(nf1)) then
        type_solve(i)=2
        nf_ref(i)=nf1
      endif
    case (2,3)
      if (ci_ref%solved(nf2)) then
        if ((int_type(nf1)==1.and.type_rotate==2).or.(int_type(nf1)==2.and.type_rotate==3)) then
          type_solve(i)=2
        else
          type_solve(i)=1
        endif
        nf_ref(i)=nf2
      endif
    case (4)
      if (ci_ref%solved(nf1).and.ci_ref%solved(nf2)) then
        if (int_type(nf1)==1) then
          type_solve(i)=3
        else
          type_solve(i)=4
        endif
        nf_ref(i)=nf2
      endif
    endselect
  endif
enddo
end

subroutine initCash_area(i)
use pgmod
integer(4) i !номер области
integer(4) type_eq,j,j1,k,k_,j_,i_,istep,ts
integer(8) n1,total_n,tn,tn1
type(TArea), pointer :: a
integer(4), parameter :: step_write=1000
integer(4) nf(max_int+1),ni,nf1,nf_ref(max_int+1),nf2
logical nf_dual(max_int+1)
integer(4) type_solve(max_int+1) !способ вычисления интеграла
                                 !0 - вычислить стандартным способом
                                 !>0 - ускоренное вычисление
                                 !1 - взять со знаком "+" интеграл nf_ref(i)
                                 !2 - взять со знаком "-" интеграл nf_ref(i)
                                 !3 - вычислить по формуле поворота для d/dx
                                      !nf(i)*cosa-nf_ref(i)*sina   nf_ref(i) - содержит номер интеграла для d/dy
                                 !4 - вычислить по формуле поворота для d/dy
                                      !nf(i)*cosa+nf_ref(i)*sina   nf_ref(i) - содержит номер интеграла для d/dx
type(TCashIntegral), pointer :: c
type(TBound), pointer :: b,b0
type(TCashIntegral), pointer :: ci_ref
logical only_numerical_int(max_int+1)
real(8) val,val2,calcintn,calcintm_,intd_z,xmc,ymc,calcint_ar_,intd_oss_z
character(80) str
logical use_dual,is_dual

if (.not.gs_use_cash) then
  call gs_print("Error initCash. Cash not use!!!")
  return
endif
use_dual=gs_use_dual_integral_solve.and.(.not.gs%const%use_numerical_int)
a=>gs%a(i)
type_eq=a%type_eq(1)
if (gs_cash_presolve) call bind_AreaConst(i)
!cash bb
call get_ni(a%need_cash_bb,ni)
if (ni>0) then
  total_n=0
  if (gs_cash_presolve) then
    do j=1,a%nb 
      b=>a%bnd(j)
      if (b%boundLineType/=1) cycle
      do j1=1,a%nb 
        b0=>a%bnd(j1)
        if (b0%boundLineType/=1) cycle
        total_n=total_n+b%npanel*b0%npanel
      enddo
    enddo
  endif
  tn=0
  tn1=0
  do j=1,a%nb !по панелям интегрирования
    b=>a%bnd(j)
    if (b%boundLineType/=1) cycle
    call init_bb_cash(i,j)
    !if (gs_cash_presolve.and.associated(b%cash,b%self_cash)) then
    if (gs_cash_presolve) then
      istep=0
      do j1=1,a%nb !по контрольным точкам
        b0=>a%bnd(j1)
        if (b0%boundLineType/=1) cycle
        c=>b%cash%bnd(j1)
        call get_ni_ki(a%need_cash_bb,ni,nf,only_numerical_int,c,nf_dual)
        n1=b%npanel*b0%npanel
        if (ni>0) then
          call gs_printInit_cash
          do k=1,ni
            call initCashIntegral_nf(c,nf(k),b%npanel,b0%npanel,b%cash%ia)
          enddo
          ci_ref=>null()
          if (a%type_rotate>=0) ci_ref=>a%a_ref%bnd(j)%cash%bnd(j1)
          call init_type_solve(i,a%type_rotate,ni,nf,type_solve,nf_ref,ci_ref)
          tn1=tn
          !$omp parallel do if (gs_use_parallel_matrix_build==1) reduction(+ : istep) private (k_,j_,i_,k,nf1,val,nf2,ts,val2,is_dual) shared (n1,b0,ni,nf,only_numerical_int,i,j,j1,tn,tn1,total_n,type_solve,nf_ref,a,use_dual,nf_dual) schedule (dynamic,1)
          do k_=0,n1-1
            istep=istep+1
            if (istep>=step_write) then
              !$omp critical (lock_screen)
              tn=tn+step_write
              istep=0
              write(str,"('cash area ',i0,' (bb): ',i0,'%')") i, tn*100/total_n
	            call write_by_time(trim(str))
              !$omp end critical (lock_screen)
            endif
            j_=k_/b0%npanel+1 !индекс панели интегрирования
            i_=mod(k_,b0%npanel)+1 !индекс панели с контрольной точкой
            do k=1,ni
              nf1=nf(k)
              is_dual=.false.
              ts=type_solve(k)
              if (ts==0) then
                if (only_numerical_int(k)) then
                  val=calcintn(j_,b0%xc(i_),b0%yc(i_),nf1,i,j)
                else
                  nf2=int_type_dual(nf1)
                  is_dual=use_dual.and.(nf2/=-1).and.nf_dual(k)
                  if (is_dual) then
                    if (nf2==-2) cycle
                    call calcintm_2(j_,i_,nf1,i,j,j1,val,val2)
                  else
                    val=calcintm_(j_,i_,nf1,i,j,j1)
                  endif
                endif
              else
                nf2=nf_ref(k)
                select case (ts)
                case (1)  !взять со знаком "+" интеграл nf_ref(i)
                  val=get_bb_cash_integral_2(a%a_ref,j,j_,j1,i_,nf2)
                case (2)  !взять со знаком "-" интеграл nf_ref(i)
                  val=-get_bb_cash_integral_2(a%a_ref,j,j_,j1,i_,nf2)
                case (3,4)  !nf(i)*cosa <+-> nf_ref(i)*sina   <+->= "-" для 3, "+" для 4
                  val=get_bb_cash_integral_2(a%a_ref,j,j_,j1,i_,nf1)
                  val2=get_bb_cash_integral_2(a%a_ref,j,j_,j1,i_,nf2)
                  if (ts==3) then
                    val=val*a%cosa-val2*a%sina
                  else
                    val=val*a%cosa+val2*a%sina
                  endif
                end select
              endif
              call set_bb_cash_integral_2(i,j,j_,j1,i_,nf1,val)
              if (is_dual) call set_bb_cash_integral_2(i,j,j_,j1,i_,nf2,val2)
            enddo
          enddo
          !omp end parallel do nowait
          do k=1,ni
            c%solved(nf(k))=.true.
          enddo
        endif
        tn=tn1+n1
      enddo
    endif
  enddo
endif
!cash ba
call get_ni(a%need_cash_ba,ni)
if (ni>0) then
  total_n=0
  if (gs_cash_presolve) then
    do j=1,a%nb 
      b=>a%bnd(j)
      if (b%boundLineType/=1) cycle
      total_n=total_n+b%npanel*a%a%ntr
    enddo
  endif
  tn=0
  tn1=0
  do j=1,a%nb !по панелям интегрирования
    b=>a%bnd(j)
    if (b%boundLineType/=1) cycle
    call init_ba_cash(i,j)
    !if (gs_cash_presolve.and.associated(b%cash,b%self_cash)) then
    if (gs_cash_presolve) then
      istep=0
      c=>b%cash%area
      call get_ni_ki(a%need_cash_ba,ni,nf,only_numerical_int,c,nf_dual)
      n1=b%npanel*a%a%ntr
      if (ni>0) then
        call gs_printInit_cash
        do k=1,ni
          call initCashIntegral_nf(c,nf(k),b%npanel,a%a%ntr,b%cash%ia)
        enddo
        ci_ref=>null()
        if (a%type_rotate>=0) ci_ref=>a%a_ref%bnd(j)%cash%area
        call init_type_solve(i,a%type_rotate,ni,nf,type_solve,nf_ref,ci_ref)
        tn1=tn
        !$omp parallel do if (gs_use_parallel_matrix_build==1) reduction(+ : istep) private (k_,j_,i_,k,nf1,val,nf2,ts,val2,is_dual) shared (n1,a,ni,nf,only_numerical_int,i,j,tn,tn1,total_n,type_solve,nf_ref,use_dual,nf_dual) schedule (dynamic,1)
        do k_=0,n1-1
          istep=istep+1
          if (istep>=step_write) then
            !$omp critical (lock_screen)
            tn=tn+step_write
            istep=0
            write(str,"('cash area ',i0,' (ba): ',i0,'%')") i, tn*100/total_n
	          call write_by_time(trim(str))
            !$omp end critical (lock_screen)
          endif
          j_=k_/a%a%ntr+1 !индекс панели интегрирования
          i_=mod(k_,a%a%ntr)+1 !индекс треугольника с контрольной точкой
          do k=1,ni
            nf1=nf(k)
            is_dual=.false.
            ts=type_solve(k)
            if (ts==0) then
              if (only_numerical_int(k)) then
                val=calcintn(j_,xmc(i_,i),ymc(i_,i),nf1,i,j)
              else
                nf2=int_type_dual(nf1)
                is_dual=use_dual.and.(nf2/=-1).and.nf_dual(k)
                if (is_dual) then
                  if (nf2==-2) cycle
                  call calcint_ar_2(j_,i_,nf1,i,j,val,val2)
                else
                  val=calcint_ar_(j_,i_,nf1,i,j)
                endif
              endif
            else
              nf2=nf_ref(k)
              select case (ts)
              case (1)  !взять со знаком "+" интеграл nf_ref(i)
                val=get_ba_cash_integral_2(a%a_ref,j,j_,i_,nf2)
              case (2)  !взять со знаком "-" интеграл nf_ref(i)
                val=-get_ba_cash_integral_2(a%a_ref,j,j_,i_,nf2)
              case (3,4)  !nf(i)*cosa <pm> nf_ref(i)*sina   <pm>= "-" для 3, "+" для 4
                val=get_ba_cash_integral_2(a%a_ref,j,j_,i_,nf1)
                val2=get_ba_cash_integral_2(a%a_ref,j,j_,i_,nf2)
                if (ts==3) then
                  val=val*a%cosa-val2*a%sina
                else
                  val=val*a%cosa+val2*a%sina
                endif
              end select
            endif
            call set_ba_cash_integral_2(i,j,j_,i_,nf1,val)
            if (is_dual) call set_ba_cash_integral_2(i,j,j_,i_,nf2,val2)
          enddo
        enddo
        !omp end parallel do nowait
        do k=1,ni
          c%solved(nf(k))=.true.
        enddo
      endif
      tn=tn1+n1
    endif
  enddo
endif
!cash ab
call get_ni(a%need_cash_ab,ni)
if (ni>0) then
  total_n=0
  if (gs_cash_presolve) then
    do j=1,a%nb 
      b=>a%bnd(j)
      if (b%boundLineType/=1) cycle
      total_n=total_n+b%npanel*a%a%ntr
    enddo
  endif
  tn=0
  tn1=0
  call init_ab_cash(i)
  !if (gs_cash_presolve.and.associated(a%a%cash,a%a%self_cash)) then
  if (gs_cash_presolve) then
    istep=0
    do j1=1,a%nb !по контрольным точкам
      b0=>a%bnd(j1)
      if (b0%boundLineType/=1) cycle
      c=>a%a%cash%bnd(j1)
      call get_ni_ki(a%need_cash_ab,ni,nf,only_numerical_int,c,nf_dual)
      n1=a%a%ntr*b0%npanel
      if (ni>0) then
        call gs_printInit_cash
        do k=1,ni
          call initCashIntegral_nf(c,nf(k),a%a%ntr,b0%npanel,a%a%cash%ia)
        enddo
        ci_ref=>null()
        if (a%type_rotate>=0) ci_ref=>a%a_ref%a%cash%bnd(j1)
        call init_type_solve(i,a%type_rotate,ni,nf,type_solve,nf_ref,ci_ref)
        tn1=tn
        !$omp parallel do if (gs_use_parallel_matrix_build==1) reduction(+ : istep) private (k_,j_,i_,k,nf1,val,nf2,ts,val2,is_dual) shared (n1,b0,ni,nf,i,j1,tn,tn1,total_n,type_solve,nf_ref,a,use_dual,nf_dual) schedule (dynamic,1)
        do k_=0,n1-1
          istep=istep+1
          if (istep>=step_write) then
            !$omp critical (lock_screen)
            tn=tn+step_write
            istep=0
            write(str,"('cash area ',i0,' (ab): ',i0,'%')") i, tn*100/total_n
	          call write_by_time(trim(str))
            !$omp end critical (lock_screen)
          endif
          j_=k_/b0%npanel+1 !индекс треугольника интегрирования
          i_=mod(k_,b0%npanel)+1 !индекс панели с контрольной точкой
          do k=1,ni
            nf1=nf(k)
            is_dual=.false.
            ts=type_solve(k)
            if (ts==0) then
              if (type_eq.eq.20) then
                val=intd_oss_z(b0%zc(i_),i,j_,nf1)
              else
                nf2=int_type_dual(nf1)
                is_dual=use_dual.and.(nf2/=-1).and.nf_dual(k)
                if (is_dual) then
                  if (nf2==-2) cycle
                  call intd_z_2(b0%zc(i_),i,j_,nf1,val,val2)
                else
                  val=intd_z(b0%zc(i_),i,j_,nf1)
                endif
              endif
            else
              nf2=nf_ref(k)
              select case (ts)
              case (1)  !взять со знаком "+" интеграл nf_ref(i)
                val=get_ab_cash_integral_2(a%a_ref,j_,j1,i_,nf2)
              case (2)  !взять со знаком "-" интеграл nf_ref(i)
                val=-get_ab_cash_integral_2(a%a_ref,j_,j1,i_,nf2)
              case (3,4)  !nf(i)*cosa <pm> nf_ref(i)*sina   <pm>= "-" для 3, "+" для 4
                val=get_ab_cash_integral_2(a%a_ref,j_,j1,i_,nf1)
                val2=get_ab_cash_integral_2(a%a_ref,j_,j1,i_,nf2)
                if (ts==3) then
                  val=val*a%cosa-val2*a%sina
                else
                  val=val*a%cosa+val2*a%sina
                endif
              end select
            endif
            call set_ab_cash_integral_2(i,j_,j1,i_,nf1,val)
            if (is_dual) call set_ab_cash_integral_2(i,j_,j1,i_,nf2,val2)
          enddo
        enddo
        !omp end parallel do nowait
        do k=1,ni
          c%solved(nf(k))=.true.
        enddo
      endif
      tn=tn1+n1
    enddo
  endif
endif
!cash aa
call get_ni(a%need_cash_aa,ni)
if (ni>0) then
  total_n=a%a%ntr**2
  tn=0
  call init_aa_cash(i)
  !if (gs_cash_presolve.and.associated(a%a%cash,a%a%self_cash)) then
  if (gs_cash_presolve) then
    istep=0
    c=>a%a%cash%area
    call get_ni_ki(a%need_cash_aa,ni,nf,only_numerical_int,c,nf_dual)
    if (ni>0) then
      call gs_printInit_cash
      n1=a%a%ntr**2
      do k=1,ni
        call initCashIntegral_nf(c,nf(k),a%a%ntr,a%a%ntr,a%a%cash%ia)
      enddo
      ci_ref=>null()
      if (a%type_rotate>=0) ci_ref=>a%a_ref%a%cash%area
      call init_type_solve(i,a%type_rotate,ni,nf,type_solve,nf_ref,ci_ref)
      !$omp parallel do if (gs_use_parallel_matrix_build==1) reduction(+ : istep) private (k_,j_,i_,k,nf1,val,nf2,ts,val2,is_dual) shared (n1,a,ni,nf,i,tn,tn1,total_n,type_solve,nf_ref,use_dual,nf_dual) schedule (dynamic,1)
      do k_=0,n1-1
        istep=istep+1
        if (istep>=step_write) then
          !$omp critical (lock_screen)
          tn=tn+step_write
          istep=0
          write(str,"('cash area ',i0,' (aa): ',i0,'%')") i, tn*100/total_n
	        call write_by_time(trim(str))
          !$omp end critical (lock_screen)
        endif
        j_=k_/a%a%ntr+1 !индекс треугольника интегрирования
        i_=mod(k_,a%a%ntr)+1 !индекс треугольника с контрольной точкой
        do k=1,ni
          nf1=nf(k)
          is_dual=.false.
          ts=type_solve(k)
          if (ts==0) then
            if (type_eq.eq.20) then
              val=intd_oss_z(a%a%zmc(i_),i,j_,nf1)
            else
              nf2=int_type_dual(nf1)
              is_dual=use_dual.and.(nf2/=-1).and.nf_dual(k)
              if (is_dual) then
                if (nf2==-2) cycle
                call intd_z_2(a%a%zmc(i_),i,j_,nf1,val,val2)
              else
                val=intd_z(a%a%zmc(i_),i,j_,nf1)
              endif
            endif
          else
            nf2=nf_ref(k)
            select case (ts)
            case (1)  !взять со знаком "+" интеграл nf_ref(i)
              val=get_aa_cash_integral_2(a%a_ref,j_,i_,nf2)
            case (2)  !взять со знаком "-" интеграл nf_ref(i)
              val=-get_aa_cash_integral_2(a%a_ref,j_,i_,nf2)
            case (3,4)  !nf(i)*cosa <pm> nf_ref(i)*sina   <pm>= "-" для 3, "+" для 4
              val=get_aa_cash_integral_2(a%a_ref,j_,i_,nf1)
              val2=get_aa_cash_integral_2(a%a_ref,j_,i_,nf2)
              if (ts==3) then
                val=val*a%cosa-val2*a%sina
              else
                val=val*a%cosa+val2*a%sina
              endif
            end select
          endif
          call set_aa_cash_integral_2(i,j_,i_,nf1,val)
          if (is_dual) call set_aa_cash_integral_2(i,j_,i_,nf2,val2)
        enddo
      enddo
      !omp end parallel do nowait
      do k=1,ni
        c%solved(nf(k))=.true.
      enddo
    endif
  endif
endif
a%need_cash_bb=.false.
a%need_cash_ab=.false.
a%need_cash_aa=.false.
a%need_cash_ba=.false.
end

subroutine get_ni(need_cash,ni)
use pgmod
integer(4) ni,j
logical need_cash(0:max_int)
ni=0
do j=0,max_int
  if (need_cash(j)) then
    ni=ni+1
  endif
enddo
end

subroutine get_ni_ki(need_cash,ni,nf,only_numerical_int,c,nf_dual)
use pgmod
integer(4) nf(max_int+1),ni,j
logical need_cash(0:max_int),only_numerical_int(max_int+1)
logical nf_dual(max_int+1)
type(TCashIntegral) c
ni=0
nf=0
nf_dual=.false.
do j=0,max_int
  if (need_cash(j).and.(.not.c%solved(j))) then
    ni=ni+1
    nf(ni)=j
  endif
enddo
only_numerical_int=.false.
do j=1,ni
  only_numerical_int(j)=nf(j)>=12.and.nf(j)<=19
  if (j<ni) then
    if (int_type_dual(nf(j))>=0.and.nf(j+1)==int_type_dual(nf(j))) nf_dual(j:j+1)=.true.
  endif
enddo
end

subroutine get_ind_costvalind(ind,j,ia)
use pgmod
integer(4) ind,j,ia
if (gs%constvalinds(ind)==0) then
  if (gs%constvala(ind)/=0.and.gs%constvala(ind)/=ia) then
    gs%constvalinds_fict_i=gs%constvalinds_fict_i+1
    gs%constvalinds(ind)=-gs%constvalinds_fict_i
    gs%constvalinds_fict(ind)=gs%constvalinds(ind)
  else
    j=j+1
    gs%constvalinds(ind)=j
  endif
endif
end

subroutine get_ind_area(ia)
!построение матрицы СЛАУ для задач включая Гельмгольца
use pgmod
integer(4) i,j,k,i1,i2,i3,i0,nu,j1,type_eq,type_eq2
integer(4) ia !номер области
type(TArea), pointer :: a
type(TBoundline), pointer :: bl
type(TBoundline2), pointer :: bl2
type(TBound), pointer :: b
type(TBoundline_GU), pointer :: gu
type(TBounLineFuncApprox), pointer :: ga,ga_ref
type(TAreaValue), pointer :: av,av2
integer(4), pointer :: arind(:)

a=>gs%a(ia)
j=a%m%ix-1 !номер текущей неизвестной
!nu=0
do i0=1,a%nb
  b=>a%bnd(i0)
  k=0 !индекс панели для массивов границы
  do i1=1,b%nline
    select case (b%boundLineType)
    case (1)
      bl=>b%line(i1)
	    do i2=1,bl%npanel
	      k=k+1
	      b%psiind(k,:)=int4_empty
        do j1=1,2 !для учета дополниетльных ГУ вместо точек коллокаций
	        do i3=1,a%nu
            if (j1==1) then
	            gu=>bl%gu(i3)
            else
              if (.not.allocated(bl.gu_add)) exit
              gu=>bl%gu_add(i3)
            endif
            select case (gu%bndu)
            case (1)
              b%psiind(k,gu%bndf)=0
              b%psiom(k,gu%bndf)=gu%bndval(i2,1)
            case (2)
              call get_ind_costvalind(gu%constvalind,j,ia)
	  	    	  b%psiind(k,gu%bndf)=gs%constvalinds(gu%constvalind)
            case (3)
              b%psiind(k,gu%bndf)=0
              gs%m%nx_fict=gs%m%nx_fict+1
            case (4)
              b%psiind(k,gu%bndf)=0
              b%psiom(k,gu%bndf)=gu%constval
            case (5)
              do i=1,gu%guglob%n
                call get_ind_costvalind(gu%guglob%gu_indsi(i),j,ia)
              enddo
              b%psiind(k,gu%bndf)=0
              gs%m%nx_fict=gs%m%nx_fict+1
            end select
          enddo
        enddo
	      do i3=1,a%umax
	        if (b%psiind(k,i3)==int4_empty) then
	  	      j=j+1
	  	      b%psiind(k,i3)=j
	  	    endif
	      enddo
	    enddo
    case (2)
      bl2=>b%line2(i1)
      do i2=1,a%umax
        ga=>bl2%ga(i2)
        if (ga%ga_ref/=0) cycle
        do i3=1,ga%n
		      if (ga%ind(i3)<0) then
            j=j+1
            ga%ind(i3)=j
		      endif
        enddo
      enddo
      !ссылки на другую аппроксимацию
      do i2=1,a%umax
        ga=>bl2%ga(i2)
        if (ga%ga_ref==0) cycle
        ga_ref=>bl2%ga(ga%ga_ref)
        if (ga_ref%ga_ref/=0) call gs_print_stop('Error ga_ref in get_ind_area!')
        do i3=1,ga%n
		      if (ga%ind(i3)<0) then
            ga%ind(i3)=ga_ref%ind(i3)
		      endif
        enddo
      enddo
    end select
  enddo
enddo
nu=a%cpp%ncp
a%m%nub=nu
type_eq=a%type_eq(1)
type_eq2=0
if (a%nu==2) type_eq2=a%type_eq(2)
if (a%haveAreaEq) then
  i0=j
  a%a%areaind=0
  do i1=1,a%umaxtr
    arind=>a%a%areaind(:,i1)
    av=>null()
    av2=>null()
    if (gs_test_areaval_eq_0) then
      select case (type_eq)
      case (5,6)
        av=>a%a%areaval(2)
      case (8,11)
        av=>a%a%areaval(3)
        if (type_eq==8.and.type_eq2==14) av2=>a%a%areaval(2)
      case (16,17)
        av=>a%a%areaval(3+i1)
      case (18:20)
        av=>a%a%areaval(4)
      case (22,24:26)
        av=>a%a%areaval(a%a%eq_ind(i1)+1)
      endselect
    endif
    !пропускаем неизвестные, перерд которыми стоят нулевые коэффициенты
    if (associated(av).and.associated(av2)) then
      do i=1,a%a%ntr
        if (av%v(i)==0.and.av2%v(i)==0) cycle  !это не проверял!!!
        j=j+1
        arind(i)=j
      enddo
    elseif (associated(av)) then
      do i=1,a%a%ntr
        if (av%v(i)==0) cycle  
        j=j+1
        arind(i)=j
      enddo
    else
      do i=1,a%a%ntr
        j=j+1
        arind(i)=j 
      enddo
    endif
  enddo
  nu=nu+j-i0
  if (a%haveAreaEq.and.(.not.a%a%cppa%inited)) call allocate_area_collocate_points(ia,nu-a%m%nub)
endif
do i=1,gs%constvaln
  if (gs%constvalinds(i)<=0.and.gs%constvala(i)==ia) then
    j=j+1
    gs%constvalinds(i)=j
  endif
enddo
a%m%nu=nu
a%m%nx=j-a%m%ix+1
end

subroutine get_ind_area2(ia,j)
!учет фиктивных переменных
use pgmod
integer(4) i,j,k,i1,i2,i3,i0,i4,k1,kk,j1
integer(4) ia !номер области
type(TArea), pointer :: a,a2
type(TBoundline), pointer :: bl,bl2
type(TBound), pointer :: b,b2
type(TBoundline_GU), pointer :: gu
type(TBoungline_GU_gen), pointer :: gug
type(TBoungline_GU_genGlobal), pointer :: guglob
type(fict_var), pointer :: fv

a=>gs%a(ia)
do i0=1,a%nb
  b=>a%bnd(i0)
  do i1=1,b%nline
    if (b%boundLineType==1) then
      bl=>b%line(i1)
      do j1=1,2
        do i3=1,a%nu
          if (j1==1) then
	          gu=>bl%gu(i3)
          else
            if (.not.allocated(bl.gu_add)) exit
            gu=>bl%gu_add(i3)
          endif
          selectcase (gu%bndu)
          case(3)
            k=bl%i_begin-1 !индекс панели для массивов границы
            gug=>gu%bndg
            kk=0
            if (gug%ia>0) then
              a2=>gs%a(gug%ia)
              if (gug%ibnd>0) then
                b2=>a2%bnd(gug%ibnd)
                if (gug%ibndl>0) then
                  bl2=>b2%line(gug%ibndl)
                  if (gug%is_direct) then
                    k1=bl2%i_begin-1
                    kk=1
                  else
                    k1=bl2%i_end+1
                    kk=-1
                  endif
                else
                  call gs_print_stop('Error in get_ind_area2!')
                endif
              else
                call gs_print_stop('Error in get_ind_area2!')
              endif
            else
              !второй области может не быть, если задествованы переменные только одной области
            endif
	          do i2=1,bl%npanel
              k=k+1
              if (kk.ne.0) k1=k1+kk
              j=j+1
              b%psiind(k,gu%bndf)=j
              i=j-gs%m%nx !индекс в массиве фиктивных переменных mfvar
              fv=>gs%m%fvar(i)
              fv%mode=1
              fv%gu_gen=>gu%bndg
              fv%i=i2
              allocate(fv%psiind(gug%n))
              allocate(fv%psiom(gug%n))
              do i4=1,gug%n
                if (gug%second_bnd(i4)) then
                  fv%psiind(i4)=b2%psiind(k1,gug%bndf(i4))
                  fv%psiom(i4)=b2%psiom(k1,gug%bndf(i4))
                else
                  fv%psiind(i4)=b%psiind(k,gug%bndf(i4))
                  fv%psiom(i4)=b%psiom(k,gug%bndf(i4))
                endif
              enddo
            enddo
          case (5)
            k=bl%i_begin-1 !индекс панели для массивов границы
            guglob=>gu%guglob
            do i2=1,bl%npanel
              k=k+1
              j=j+1
              b%psiind(k,gu%bndf)=j
              i=j-gs%m%nx !индекс в массиве фиктивных переменных mfvar
              fv=>gs%m%fvar(i)
              fv%mode=2
              fv%gu_genGlob=>gu%guglob
              fv%i=i2
              allocate(fv%psiind(guglob%n))
              do i4=1,guglob%n
                fv%psiind(i4)=gs%constvalinds(guglob%gu_indsi(i4))
              enddo
            enddo
          endselect
        enddo
      enddo
    endif
  enddo
enddo
end

subroutine get_ind_area3
!переделываем фиктивные индексы глобальных неизвестных
use pgmod
integer(4) i,ia,j,k,l
type(TArea), pointer :: a
type(TBound), pointer :: b
do i=1,gs%constvaln
  if (gs%constvalinds_fict(i)==0) cycle
  do ia=1,gs%na
    a=>gs%a(ia)
    do j=1,a%nb
      b=>a%bnd(j)
      do k=1,b%npanel
        do l=1,a%umax
          if (b%psiind(k,l)==gs%constvalinds_fict(i)) b%psiind(k,l)=gs%constvalinds(i)
        enddo
      enddo
    enddo
  enddo
enddo
end

subroutine init_col_info(a,ci,ind,buff_c,n,need_fict)
use pgmod
Type(TArea) a
type(col_info) ci
integer(4) n,i
logical ind(n)
integer(4) buff_c(n)
logical need_fict
ci%ncol=0
ci%icol_area=0
do i=1,gs%m%nx
  if (ind(i).and.(i<a%m%ix.or.i>a%m%ix2)) then
    ci%ncol=ci%ncol+1
    if (i<a%m%ix) ci%icol_area=ci%ncol
    buff_c(ci%ncol)=i
  endif
enddo
if (ci%ncol>0) then
  allocate(ci%col(ci%ncol))
  ci%col(1:ci%ncol)=buff_c(1:ci%ncol)
endif
ci%ncol_fict=0
if (need_fict) then
  do i=gs%m%ix_fict,gs%m%nx_all
    if (ind(i)) then
      ci%ncol_fict=ci%ncol_fict+1
      buff_c(ci%ncol_fict)=i
    endif
  enddo
  if (ci%ncol_fict>0) then
    allocate(ci%col_fict(ci%ncol_fict))
    ci%col_fict(1:ci%ncol_fict)=buff_c(1:ci%ncol_fict)
  endif
endif
end

subroutine prepare_global_sparsem
use pgmod
integer(4) ia,i0,i,j,k,i1,i2,i3,nclos,n
Type(TArea), pointer :: a
Type(TBound), pointer :: b
type(TBoundline2), pointer :: bl2
logical, allocatable :: ind(:)
integer(4), allocatable :: buff_c(:)
type(fict_var), pointer :: fv
type(TBounLineFuncApprox), pointer :: ga
type(sparse_matrix_row), pointer :: mrow
allocate(ind(gs%m%nx_all))
allocate(buff_c(gs%m%nx_all))
do ia=1,gs%na 
  a=>gs%a(ia)
  ind=.false.
  do i0=1,a%nb
    b=>a%bnd(i0)
    select case (b%boundLineType)
    case (1)
      do i=1,b%npanel
        do j=1,a%umax
          k=b%psiind(i,j)
	        if (k>0) then
            ind(k)=.true.
            if (k>=gs%m%ix_fict) then
              k=k-gs%m%nx !индекс в массиве фиктивных переменных mfvar
              fv=>gs%m%fvar(k)
              selectcase (fv%mode)
              case (1)
                n=fv%gu_gen%n
              case (2)
                n=fv%gu_genGlob%n
              endselect
              do i1=1,n
                k=fv%psiind(i1)
                if (k>0) ind(k)=.true.
              enddo
            endif
          endif
        enddo
      enddo
    case (2)
      do i1=1,b%nline
        bl2=>b%line2(i1)
        do i2=1,a%umax
          ga=>bl2%ga(i2)
          do i3=1,ga%n
		        k=ga%ind(i3)
            if (k>0) ind(k)=.true.
          enddo
        enddo
      enddo
    end select
  enddo
  if (a%haveAreaEq) then
    do i=1,a%a%ntr
      do i1=1,a%umaxtr
        k=a%a%areaind(i,i1)
        if (k>0) ind(k)=.true.
      enddo
    enddo
  endif
  if (gs%m%matrix_store_type==1.and.use_global_sparsem) then
    a%m%area_sm_count=1 !для учета обязательного диагонального элемента
    do i=1,gs%m%nx
      if (ind(i)) a%m%area_sm_count=a%m%area_sm_count+1
    enddo
  endif
  call init_col_info(a,a%m%self_ci,ind,buff_c,gs%m%nx_all,.true.)
enddo
deallocate(ind,buff_c)
if (gs%m%matrix_store_type==1.and.use_global_sparsem) then
  !подсчитываем, сколько нужно памяти
  forall (i=1:gs%m%nu) gs%m%sparse%mrows(i)%sm_count=0
  k=1
  do ia=1,gs%na 
    a=>gs%a(ia)
    do i=1,a%cpp%ncp
      mrow=>gs%m%sparse%mrows(a%cpp%j(i))
      mrow%sm_ibegin=k
      mrow%sm_count=a%m%area_sm_count
      k=k+a%m%area_sm_count
    enddo
    if (a%haveAreaEq) then
      do i=1,a%a%cppa%ncp
        mrow=>gs%m%sparse%mrows(a%a%cppa%j(i))
        mrow%sm_ibegin=k
        mrow%sm_count=a%m%area_sm_count
        k=k+a%m%area_sm_count
      enddo
    endif
  enddo
  nclos=gs%m%count_closing
  if (nclos==0) nclos=gs%m%nx
  do i=1,gs%m%nu  !оставшиеся уравнения для замыкания
    mrow=>gs%m%sparse%mrows(i)
    if (mrow%sm_count==0) then
      mrow%sm_ibegin=k
      mrow%sm_count=nclos
      k=k+nclos
    endif
  enddo
  gs%m%nc_global=k-1
  allocate(gs%m%sparse%m(gs%m%nc_global))
  allocate(gs%m%sparse%c(gs%m%nc_global))
  call sparse_memory(0)
endif
end

subroutine sparse_memory(mode)
use pgmod
integer(4) mode
real(8) nc
if (mode==0) then
  nc=gs%m%nc_global
else
  nc=gs%m%sparse%nc
endif
gs%m%sparse_percent=nc*100/gs%m%nx**2
gs%m%memory_mb=(nc*12+(gs%m%nu+1)*4)/1024**2
end

subroutine init_gs_time
use pgmod
real(8) rtc
gs_time0=rtc()
gs_time1=gs_time0-gs_time_interval*2
end

function test_gs_time
use pgmod
real(8) rtc,time
logical test_gs_time
time=rtc()
test_gs_time=time-gs_time1>gs_time_interval
end

function write_by_timef(str)
use pgmod
logical write_by_timef,test_gs_time
character(*) str
write_by_timef=test_gs_time()
if (write_by_timef) then
  call gs_print(str)
  call end_write_by_time
endif
end

subroutine end_write_by_time
use pgmod
real(8) rtc,time
time=rtc()
gs_time1=time
end

function write_by_timef2(str,test)
use pgmod
logical write_by_timef,test,write_by_timef2
character(*) str
if (test) then
  write_by_timef2=write_by_timef(str)
else
  write_by_timef2=.true.
  call gs_print(str)
endif
end

subroutine write_by_time(str)
use pgmod
character(*) str
logical write_by_timef,res
res=write_by_timef(str)
end

function get_total_time
use pgmod
real(8) get_total_time,rtc
get_total_time=rtc()-gs_time0
end

subroutine print_total_time
use pgmod
real(8) get_total_time
call gs_print('time=  '//ftoa(get_total_time()))
end

subroutine init_cpp_j
!инициализируем номера уравнений у контрольных точек
use pgmod
integer(4) ia,i,icp
Type(TArea), pointer :: a
do ia=1,gs%na 
  a=>gs%a(ia)
  i=a%m%iu-1  !номер уравнения
  do icp=1,a%cpp%ncp
    i=i+1
    a%cpp%j(icp)=i
  enddo
  if (a%haveAreaEq) then
    do icp=1,a%a%cppa%ncp
      i=i+1
      a%a%cppa%j(icp)=i
    enddo
  endif
enddo
end

subroutine get_matrix_area(ia,func_inf)
!построение матрицы СЛАУ для задач включая Гельмгольца
use pgmod
!use omp_lib
integer(4) icp
integer(4) ia !номер области
real(8) func_inf
external func_inf
type(TArea), pointer :: a
type(areatype), pointer :: aa
character(75) str
logical test_gs_time
a=>gs%a(ia)
a%m%ci=>a%m%self_ci
call bind_AreaConst(ia)
aa=>a%a
!цикл по граничным элементам
!call omp_set_num_threads(8)
!$omp parallel do if (gs_use_parallel_matrix_build==1) private (icp,str) shared (ia,a,gs) schedule (dynamic,1)
do icp=1,a%cpp%ncp
  !$omp critical (lock_screen)
  if (test_gs_time()) then
    write(str,"('Matrix [domain ', i0, ' (', I0, '), bound ', i0, ' (', I0, ')]. Bound El.Eq ', I0, ' (', I0, ')')") ia, gs%na, a%cpp%ibnd(icp), a%nb, icp,a%cpp%ncp
	  call gs_print(trim(str))
    call end_write_by_time
  endif
  !$omp end critical (lock_screen)
  call get_matrix_area_bp(ia,icp,func_inf)
enddo
!omp end parallel do nowait
if (a%haveAreaEq) then
  select case (aa%areaEq)  
  case (2,5:7,12:14,16:20,22,24:26)
    !цикл по площадным элементам
    !$omp parallel do if (gs_use_parallel_matrix_build==1) private (icp,str) shared (ia,a,gs) schedule (dynamic,1)
    do icp=1,aa%cppa%ncp
      !$omp critical (lock_screen)
      if (test_gs_time()) then
        write(str,"('Matrix [domain ', i0, ' (', I0, ')]. Area El.Eq ', I0, ' (', I0, ')')") ia, gs%na, icp, aa%cppa%ncp
	      call gs_print(trim(str))
        call end_write_by_time
      endif
      !$omp end critical (lock_screen)
  	call get_matrix_area_ap(ia,icp,func_inf)
    enddo
    !omp end parallel do nowait
  end select
endif
end

subroutine get_matrix_area_bp(ia,icp,func_inf)
!построение матрицы СЛАУ для задач включая Гельмгольца
use pgmod
integer(4) i,j,k,i1,j1,k1,i0,i2,j0,j2,i3,type_eq,type_eq2,icp,omp_get_thread_num,ith
integer(4) ia !номер области
real(8) calcintm,intd_m,calcint_ar,calcintn2,calcintn_ar,fff,fff2,sss,sss2,kk
real(8) cp_calcintm,cp_calcintm1,cp_calcintm2
real(8) val,val2,k2,val3
real(8) intd_oss_m
logical qlog,qlog2
real(8) func_inf
external func_inf
type(TArea), pointer :: a
type(areatype), pointer :: aa
type(TBound), pointer :: b,b0
type(TBoundline2), pointer :: bndl2
type(TBounLineFuncApprox), pointer :: ga
type(TBoundline_Collocate), pointer :: gc
real(8), pointer :: m_eq(:)
logical type_eq22,type_eq22orNS
real(8), pointer :: m_b
type(TOMP_buffer), pointer :: buff
type(TAreaValue), pointer :: av,av2,av3
logical use_dual,have_fff2
ith=omp_get_thread_num()
buff=>gs%m%omp_buff%b(ith)
m_eq=>buff%eq
m_b=>buff%b
use_dual=gs_use_dual_integral_solve.and.(.not.gs%const%use_numerical_int)
a=>gs%a(ia)
aa=>a%a
i0=a%cpp%ibnd(icp) !индекс границы с контрольной точкой
b0=>a%bnd(i0)
i2=a%cpp%iu(icp) !уравнение - соотношение Грина (главное или второе)
i=a%cpp%j(icp)   !индекс уравнения в матрице
buff%i=i
type_eq=a%type_eq(i2)
type_eq22orNS=(type_eq==22).or.(type_eq==24).or.(type_eq==25)
if (i2==1.and.a%nu==2) type_eq2=a%type_eq(2)
k=0
if ((type_eq.eq.2).or.(type_eq.eq.7).or.(type_eq.eq.12).or.(type_eq.eq.13).or.(type_eq.eq.14).or.(type_eq.eq.23)) k=2   
!m_eq=d0 !перенесли в add_equation_to_matrix_thread
m_b=d0
if (b0%boundLineType==2) then
  !контрольные точки на круговых частичках
  i1=a%cpp%ibndl(icp)
  gc=>b0%line2(i1)%gc(i2)
  k1=a%cpp%i(icp)

  fff=func_inf(dreal(gc%z(k1)),dimag(gc%z(k1)),i2,0,ia) !поведение на бесконечности
  if (fff/=d0) m_b=m_b-pi2*fff

  ga=>b0%line2(i1)%ga(1+k)
  if (ga%typea==1) then
    do i3=1,ga%n
      j1=ga%ind(i3)
      val=ga%val(i3)
      if (j1.ne.0.or.val.ne.d0) then
        select case (i3)
        case (1)
          kk=d1
        case (2)
          kk=dsin(gc%g(k1))
        case (3)
          kk=dcos(gc%g(k1))
        end select
        kk=kk*ga%kk(i3)
        if (j1.eq.0) then
          m_b=m_b+val*kk*pi
        else
          m_eq(j1)=m_eq(j1)-pi*kk
        endif
      endif
    enddo
  endif
  do j0=1,a%nb
    b=>a%bnd(j0)
    if (b%boundLineType==2) then
      !интегрирование по кругу
      do j=1,b%nline
        bndl2=>b%line2(j)
        !psi или eta
        ga=>bndl2%ga(1+k)
        if (ga%typea==1) then
          have_fff2=.false.
          do i3=1,ga%n !индекс коэффициента в разложении Фурье
            j1=ga%ind(i3)
            val=ga%val(i3)
            if (j1.ne.0.or.val.ne.d0) then
              select case (i3)
              case (1)
                fff=cp_calcintm2(j,i1,i2,k1,1,ia,j0,i0)
              case (2)
                if (use_dual) then
                  call cp_calcintm2_2(j,i1,i2,k1,22,ia,j0,i0,fff,fff2)
                  have_fff2=.true.
                else
                  fff=cp_calcintm2(j,i1,i2,k1,22,ia,j0,i0)
                endif
              case (3)
                if (use_dual.and.have_fff2) then
                  fff=fff2
                else
                  fff=cp_calcintm2(j,i1,i2,k1,23,ia,j0,i0)
                endif
              case default
                !!!сюда не должны приходить
                call gs_print_stop("Wrong parameters for boundLineType=2 in get_matrix_area")
              end select
              fff=fff*ga%kk(i3)
              if (j1.eq.0) then
                m_b=m_b-val*fff
              else
                m_eq(j1)=m_eq(j1)+fff
              endif
            endif
          enddo
        endif
        !dpsi/dn или deta/dn
        ga=>bndl2%ga(2+k)
        if (ga%typea==1) then
          have_fff2=.false.
          do i3=1,ga%n !индекс коэффициента в разложении Фурье
            j1=ga%ind(i3)
            val=ga%val(i3)
            if (j1.ne.0.or.val.ne.d0) then
              select case (i3)
              case (1)
                fff=cp_calcintm2(j,i1,i2,k1,0,ia,j0,i0)
              case (2)
                if (use_dual) then
                  call cp_calcintm2_2(j,i1,i2,k1,20,ia,j0,i0,fff,fff2)
                  have_fff2=.true.
                else
                  fff=cp_calcintm2(j,i1,i2,k1,20,ia,j0,i0)
                endif
              case (3)
                if (use_dual.and.have_fff2) then
                  fff=fff2
                else
                  fff=cp_calcintm2(j,i1,i2,k1,21,ia,j0,i0)
                endif
              case default
                !!!сюда не должны приходить
                call gs_print_stop("Wrong parameters for boundLineType=2 in get_matrix_area")
              end select
              fff=fff*ga%kk(i3)
              if (j1.eq.0) then
                m_b=m_b+val*fff
              else
                m_eq(j1)=m_eq(j1)-fff
              endif
            endif
          enddo
        endif
      enddo 
      !бигармоническое уравнение
      if (type_eq.eq.3) then
        do j=1,b%nline
          bndl2=>b%line2(j)
          !eta=-om
          ga=>bndl2%ga(3)
          if (ga%typea==1) then
            have_fff2=.false.
            do i3=1,ga%n !индекс коэффициента в разложении Фурье
              j1=ga%ind(i3)
              val=ga%val(i3)
              if (j1.ne.0.or.val.ne.d0) then
                select case (i3)
                case (1)
                  fff=cp_calcintm2(j,i1,i2,k1,3,ia,j0,i0)
                case (2)
                  if (use_dual) then
                    call cp_calcintm2_2(j,i1,i2,k1,26,ia,j0,i0,fff,fff2)
                    have_fff2=.true.
                  else
                    fff=cp_calcintm2(j,i1,i2,k1,26,ia,j0,i0)
                  endif
                case (3)
                  if (use_dual.and.have_fff2) then
                    fff=fff2
                  else
                    fff=cp_calcintm2(j,i1,i2,k1,27,ia,j0,i0)
                  endif
                case default
                  !!!сюда не должны приходить
                  call gs_print_stop("Wrong parameters for boundLineType=2 in get_matrix_area")
                end select
                fff=fff*ga%kk(i3)
                if (j1.eq.0) then
                  m_b=m_b-val*fff
                else
                  m_eq(j1)=m_eq(j1)+fff
                endif
              endif
            enddo
          endif
          !deta/dn=-dom/dn
          ga=>bndl2%ga(4)
          if (ga%typea==1) then
            have_fff2=.false.
            do i3=1,ga%n !индекс коэффициента в разложении Фурье
              j1=ga%ind(i3)
              val=ga%val(i3)
              if (j1.ne.0.or.val.ne.d0) then
                select case (i3)
                case (1)
                  fff=cp_calcintm2(j,i1,i2,k1,2,ia,j0,i0)
                case (2)
                  if (use_dual) then
                    call cp_calcintm2_2(j,i1,i2,k1,24,ia,j0,i0,fff,fff2)
                    have_fff2=.true.
                  else
                    fff=cp_calcintm2(j,i1,i2,k1,24,ia,j0,i0)
                  endif
                case (3)
                  if (use_dual.and.have_fff2) then
                    fff=fff2
                  else
                    fff=cp_calcintm2(j,i1,i2,k1,25,ia,j0,i0)
                  endif
                case default
                  !!!сюда не должны приходить
                  call gs_print_stop("Wrong parameters for boundLineType=2 in get_matrix_area")
                end select
                fff=fff*ga%kk(i3)
                if (j1.eq.0) then
                  m_b=m_b+val*fff
                else
                  m_eq(j1)=m_eq(j1)-fff
                endif
              endif
            enddo
          endif
        enddo
      endif
    else
      !интегрирование по панелям
      !b%boundLineType==1
      do j=1,b%npanel
        j1=b%psiind(j,1+k)
        val=b%psiom(j,1+k)
        if (j1.ne.0.or.val.ne.d0) then
          fff=cp_calcintm1(j,i1,i2,k1,1,ia,j0,i0)
          if (j1.eq.0) then
            m_b=m_b-val*fff
          else
            m_eq(j1)=m_eq(j1)+fff
          endif
        endif
        j1=b%psiind(j,2+k)
        val=b%psiom(j,2+k)
        if (j1.ne.0.or.val.ne.d0) then
          fff=cp_calcintm1(j,i1,i2,k1,0,ia,j0,i0)
          if (j1.eq.0) then
            m_b=m_b+val*fff
          else
            m_eq(j1)=m_eq(j1)-fff
          endif
        endif
      enddo
      !бигармоническое уравнение
      if (type_eq.eq.3) then
        do j=1,b%npanel
          j1=b%psiind(j,3)
          val=b%psiom(j,3)
          if (j1.eq.0) then
            if (val.ne.d0) m_b=m_b-val*cp_calcintm1(j,i1,i2,k1,3,ia,j0,i0)  
          else
            m_eq(j1)=m_eq(j1)+cp_calcintm1(j,i1,i2,k1,3,ia,j0,i0)  
          endif
          j1=b%psiind(j,4)
          val=b%psiom(j,4)
          if (j1.eq.0) then
            if (val.ne.d0) m_b=m_b+val*cp_calcintm1(j,i1,i2,k1,2,ia,j0,i0)
          else
            m_eq(j1)=m_eq(j1)-cp_calcintm1(j,i1,i2,k1,2,ia,j0,i0)
          endif
        enddo
      endif
    endif
  enddo
elseif (b0%boundLineType==1) then
  !контрольная точка на обычной панели
  i1=a%cpp%i(icp) !индекс панели с контрольной точкой
  j1=b0%psiind(i1,1+k)

  fff=func_inf(b0%xc(i1),b0%yc(i1),i2,0,ia) !поведение на бесконечности
  if (fff/=d0) m_b=m_b-pi2*fff

  val=b0%psiom(i1,1+k)
  if (j1.eq.0) then
    if (val/=d0) m_b=m_b+val*pi
  else
    m_eq(j1)=m_eq(j1)-pi
  endif
  do j0=1,a%nb
    b=>a%bnd(j0)
    if (b%boundLineType==2) then
      do j=1,b%nline
        bndl2=>b%line2(j)
        !psi или eta
        ga=>bndl2%ga(1+k)
        if (ga%typea==1) then
          have_fff2=.false.
          do k1=1,ga%n !индекс коэффициента в разложении Фурье
            j1=ga%ind(k1)
            val=ga%val(k1)
            if (j1.ne.0.or.val.ne.d0) then
              select case (k1)
              case (1)
                fff=cp_calcintm(j,i1,1,ia,j0,i0)
              case (2)
                if (use_dual) then
                  call cp_calcintm_2(j,i1,22,ia,j0,i0,fff,fff2)
                  have_fff2=.true.
                else
                  fff=cp_calcintm(j,i1,22,ia,j0,i0)
                endif
              case (3)
                if (use_dual.and.have_fff2) then
                  fff=fff2
                else
                  fff=cp_calcintm(j,i1,23,ia,j0,i0)
                endif
              case  default
                !!!сюда не должны приходить
                call gs_print_stop("Wrong parameters for boundLineType=2 in get_matrix_area")
              end select
              fff=fff*ga%kk(k1)
              if (j1.eq.0) then
                m_b=m_b-val*fff
              else
                m_eq(j1)=m_eq(j1)+fff
              endif
            endif
          enddo
        endif
        !dpsi/dn или deta/dn
        ga=>bndl2%ga(2+k)
        if (ga%typea==1) then
          have_fff2=.false.
          do k1=1,ga%n !индекс коэффициента в разложении Фурье
            j1=ga%ind(k1)
            val=ga%val(k1)
            if (j1.ne.0.or.val.ne.d0) then
              select case (k1)
              case (1)
                fff=cp_calcintm(j,i1,0,ia,j0,i0)
              case (2)
                if (use_dual) then
                  call cp_calcintm_2(j,i1,20,ia,j0,i0,fff,fff2)
                  have_fff2=.true.
                else
                  fff=cp_calcintm(j,i1,20,ia,j0,i0)
                endif
              case (3)
                if (use_dual.and.have_fff2) then
                  fff=fff2
                else
                  fff=cp_calcintm(j,i1,21,ia,j0,i0)
                endif
              case default
                !!!сюда не должны приходить
                call gs_print_stop("Wrong parameters for boundLineType=2 in get_matrix_area")
              end select
              fff=fff*ga%kk(k1)
              if (j1.eq.0) then
                m_b=m_b+val*fff
              else
                m_eq(j1)=m_eq(j1)-fff
              endif
            endif
          enddo
        endif
      enddo
    else
      !b%boundLineType==1
      do j=1,b%npanel
        j1=b%psiind(j,1+k)
        val=b%psiom(j,1+k)
        if (j1.ne.0.or.val.ne.d0) then
          select case (type_eq)
          case (9,12)
            fff=calcintn2(j,i1,13,ia,j0,i0)
          case (10,13)
            fff=calcintn2(j,i1,15,ia,j0,i0)
          case default
            fff=calcintm(j,i1,1,ia,j0,i0)
          end select
          if (j1.eq.0) then
            m_b=m_b-val*fff
          else
            m_eq(j1)=m_eq(j1)+fff
          endif
        endif
        j1=b%psiind(j,2+k)
        val=b%psiom(j,2+k)
        if (j1.ne.0.or.val.ne.d0) then
          select case (type_eq)
          case (9,12)
            fff=calcintn2(j,i1,12,ia,j0,i0)
          case (10,13)
            fff=calcintn2(j,i1,14,ia,j0,i0)
          case default
            fff=calcintm(j,i1,0,ia,j0,i0)
          end select
          if (j1.eq.0) then
            m_b=m_b+val*fff
          else
            m_eq(j1)=m_eq(j1)-fff
          endif
        endif
      enddo
    endif
  enddo
  !бигармоническое уравнение
  if ((type_eq.eq.3).or.(type_eq==21).or.type_eq22orNS) then
    do j0=1,a%nb
      b=>a%bnd(j0)
      if (b%boundLineType==2) then
        do j=1,b%nline
          bndl2=>b%line2(j)
          !eta=-om
          ga=>bndl2%ga(3)
          if (ga%typea==1) then
            have_fff2=.false.
            do k1=1,ga%n !индекс коэффициента в разложении Фурье
              j1=ga%ind(k1)
              val=ga%val(k1)
              if (j1.ne.0.or.val.ne.d0) then
                select case (k1)
                case (1)
                  fff=cp_calcintm(j,i1,3,ia,j0,i0)
                case (2)
                  if (use_dual) then
                    call cp_calcintm_2(j,i1,26,ia,j0,i0,fff,fff2)
                    have_fff2=.true.
                  else
                    fff=cp_calcintm(j,i1,26,ia,j0,i0)
                  endif
                case (3)
                  if (use_dual.and.have_fff2) then
                    fff=fff2
                  else
                    fff=cp_calcintm(j,i1,27,ia,j0,i0)
                  endif
                case default
                  !!!сюда не должны приходить
                  call gs_print_stop("Wrong parameters for boundLineType=2 in get_matrix_area")
                end select
                fff=fff*ga%kk(k1)
                if (j1.eq.0) then
                  m_b=m_b-val*fff
                else
                  m_eq(j1)=m_eq(j1)+fff
                endif
              endif
            enddo
          endif
          !deta/dn=-dom/dn
          ga=>bndl2%ga(4)
          if (ga%typea==1) then
            have_fff2=.false.
            do k1=1,ga%n !индекс коэффициента в разложении Фурье
              j1=ga%ind(k1)
              val=ga%val(k1)
              if (j1.ne.0.or.val.ne.d0) then
                select case (k1)
                case (1)
                  fff=cp_calcintm(j,i1,2,ia,j0,i0)
                case (2)
                  if (use_dual) then
                    call cp_calcintm_2(j,i1,24,ia,j0,i0,fff,fff2)
                    have_fff2=.true.
                  else
                    fff=cp_calcintm(j,i1,24,ia,j0,i0)
                  endif
                case (3)
                  if (use_dual.and.have_fff2) then
                    fff=fff2
                  else
                    fff=cp_calcintm(j,i1,25,ia,j0,i0)
                  endif
                case default
                  !!!сюда не должны приходить
                  call gs_print_stop("Wrong parameters for boundLineType=2 in get_matrix_area")
                end select
                fff=fff*ga%kk(k1)
                if (j1.eq.0) then
                  m_b=m_b+val*fff
                else
                  m_eq(j1)=m_eq(j1)-fff
                endif
              endif
            enddo
          endif
        enddo
      else
        !b%boundLineType==1
        do j=1,b%npanel
          j1=b%psiind(j,3)
          val=b%psiom(j,3)
          if (j1.ne.0.or.val.ne.d0) then
            fff=calcintm(j,i1,3,ia,j0,i0)
            if (j1.eq.0) then
              m_b=m_b-val*fff
            else
              m_eq(j1)=m_eq(j1)+fff
            endif
          endif
          j1=b%psiind(j,4)
          val=b%psiom(j,4)
          if (j1.ne.0.or.val.ne.d0) then
            fff=calcintm(j,i1,2,ia,j0,i0)
            if (j1.eq.0) then
              m_b=m_b+val*fff
            else
              m_eq(j1)=m_eq(j1)-fff
            endif
          endif
        enddo
      endif
    enddo
  endif
  !уравнение Бринкмана через единую функцию Грина
  if (type_eq.eq.15) then
    k2=d1/a%const%k_helm**2
    do j0=1,a%nb
      b=>a%bnd(j0)
      do j=1,b%npanel
        j1=b%psiind(j,3)
        val=b%psiom(j,3)
        if (j1.ne.0.or.val.ne.d0) then
          fff=calcintn2(j,i1,15,ia,j0,i0)
          fff=fff-calcintm(j,i1,1,ia,j0,i0)
          fff=fff*k2
          if (j1.eq.0) then
            m_b=m_b-val*fff
          else
            m_eq(j1)=m_eq(j1)+fff
          endif
        endif
        j1=b%psiind(j,4)
        val=b%psiom(j,4)
        if (j1.ne.0.or.val.ne.d0) then
          fff=calcintn2(j,i1,14,ia,j0,i0)
          fff=fff-calcintm(j,i1,0,ia,j0,i0)
          fff=fff*k2
          if (j1.eq.0) then
            m_b=m_b+val*fff
          else
            m_eq(j1)=m_eq(j1)-fff
          endif
        endif
      enddo
    enddo
  endif
  !неоднородные уравнения
  type_eq22=.false.
  if ((type_eq.eq.23).or.(type_eq.eq.26)) type_eq22=a%eq_var(1)
  if ((type_eq.eq.4).or.(type_eq.eq.6).or.(type_eq.eq.7).or.(type_eq22)) then
    av=>aa%areaval(1)
    do j=1,aa%ntr
      val=av%v(j)
      if (val.ne.d0) m_b=m_b-val*intd_m(ia,j,i1,i0,0) !intd_an(b0%zc(i1),zmz(1,j,ia),zmz(2,j,ia),zmz(3,j,ia),0)
    enddo
  endif
  !неоднородное бигармоничское уравнение
  type_eq22=.false.
  if (type_eq22orNS) type_eq22=a%eq_var(1)
  if ((type_eq.eq.21).or.(type_eq22)) then
    av=>aa%areaval(1)
    do j=1,aa%ntr
      val=av%v(j)
      if (val.ne.d0) m_b=m_b-val*intd_m(ia,j,i1,i0,2) !intd_an(b0%zc(i1),zmz(1,j,ia),zmz(2,j,ia),zmz(3,j,ia),2)
    enddo
  endif
  !Уравнение Гельмгольца
  select case (type_eq)
  case (5,6,14)
    av=>aa%areaval(2)
    do j=1,aa%ntr
      val=av%v(j)
      if (val.ne.d0) then
        j1=aa%areaind(j,1)
        m_eq(j1)=m_eq(j1)-val*intd_m(ia,j,i1,i0,0) !intd_an(b0%zc(i1),zmz(1,j,ia),zmz(2,j,ia),zmz(3,j,ia),0)
      endif
    enddo
  end select
  !система вида \Delta\psi=v(x,y)\omega  и уравнение на \omega (\omega  в области неизвестна)
  if (type_eq.eq.8) then
    av=>aa%areaval(3)
    do j=1,aa%ntr
      val=av%v(j)
      if (val.ne.d0) then
        j1=aa%areaind(j,1)
        m_eq(j1)=m_eq(j1)+val*intd_m(ia,j,i1,i0,0) !intd_an(b0%zc(i1),zmz(1,j,ia),zmz(2,j,ia),zmz(3,j,ia),0)
      endif
    enddo
  endif
  !система вида \Delta\psi=v(x,y)\omega  и уравнение на \omega
  if (type_eq.eq.11) then 
    av=>aa%areaval(3)
    do j0=1,a%nb
      b=>a%bnd(j0)
      do j=1,b%npanel
        sss=d0
        sss2=d0
        j1=b%psiind(j,3)
        j2=b%psiind(j,4)
        val=b%psiom(j,3)
        val2=b%psiom(j,4)
        qlog=.not.(j1.eq.0.and.val.eq.d0)
        qlog2=.not.(j2.eq.0.and.val2.eq.d0)
        if (qlog.or.qlog2) then
          do k1=1,aa%ntr
            fff=av%v(k1)
            if (fff.ne.d0) fff=fff*intd_m(ia,k1,i1,i0,0) !intd_an(b0%zc(i1),zmz(1,k1,ia),zmz(2,k1,ia),zmz(3,k1,ia),0)
            if (fff.ne.d0) then
              fff2=fff
              select case (type_eq2)
              case (2)
                if (qlog) fff=fff*calcint_ar(j,k1,1,ia,j0) !calcint(j,xmc(k1,ia),ymc(k1,ia),1,ia,j0)
                if (qlog2) fff2=fff2*calcint_ar(j,k1,0,ia,j0) !calcint(j,xmc(k1,ia),ymc(k1,ia),0,ia,j0)
              case (12)
                if (qlog) fff=fff*calcintn_ar(j,k1,13,ia,j0) !calcintn(j,xmc(k1,ia),ymc(k1,ia),13,ia,j0)
                if (qlog2) fff2=fff2*calcintn_ar(j,k1,12,ia,j0) !calcintn(j,xmc(k1,ia),ymc(k1,ia),12,ia,j0)
              case (13)
                if (qlog) fff=fff*calcintn_ar(j,k1,15,ia,j0) !calcintn(j,xmc(k1,ia),ymc(k1,ia),15,ia,j0)
                if (qlog2) fff2=fff2*calcintn_ar(j,k1,14,ia,j0) !calcintn(j,xmc(k1,ia),ymc(k1,ia),14,ia,j0)
              case default
                call gs_print_stop("Error for second equation in system in get_matrix_area")
              end select
              if (qlog) sss=sss+fff
              if (qlog2) sss2=sss2+fff2
            endif
          enddo
          if (qlog) then
            sss=sss/pi2 
            if (j1.eq.0) then
              m_b=m_b-val*sss
            else
              m_eq(j1)=m_eq(j1)+sss
            endif
          endif
          if (qlog2) then
            sss2=sss2/pi2
            if (j2.eq.0) then
              m_b=m_b+val2*sss2
            else
              m_eq(j2)=m_eq(j2)-sss2
            endif
          endif
        endif
      enddo
    enddo
  endif
  !уравнение переноса и осесимметричные уравнения
  select case (type_eq)
  case (16:20)
    av=>null()
    if ((type_eq.eq.17).or.(type_eq.eq.18)) av=>aa%areaval(1)
    av2=>aa%areaval(4)
    av3=>null()
    if ((type_eq.eq.16).or.(type_eq.eq.17)) av3=>aa%areaval(5)
    do j=1,aa%ntr
      val=d0
      val3=d0
      if ((type_eq.eq.17).or.(type_eq.eq.18)) val=av%v(j)
      val2=av2%v(j)
      if ((type_eq.eq.16).or.(type_eq.eq.17)) val3=av3%v(j)
      if ((val.ne.d0).or.(val2.ne.d0).or.(val3.ne.d0)) then
        if (type_eq.eq.20) then
          fff=intd_oss_m(ia,j,i1,i0,0) !intd_oss(b0%zc(i1),zmz(1,j,ia),zmz(2,j,ia),zmz(3,j,ia),0)
        else
          fff=intd_m(ia,j,i1,i0,0) !intd_an(b0%zc(i1),zmz(1,j,ia),zmz(2,j,ia),zmz(3,j,ia),0)
        endif
        if (val.ne.d0) m_b=m_b-val*fff
        if (val2.ne.d0) then
          j1=aa%areaind(j,1)
          m_eq(j1)=m_eq(j1)+val2*fff
        endif
        if (val3.ne.d0) then
          j1=aa%areaind(j,2)
          m_eq(j1)=m_eq(j1)+val3*fff
        endif
      endif
    enddo
  endselect
  !неоднородное бигармоническое уравнение с неизвестными в правой части
  if (type_eq22orNS.or.(type_eq.eq.23).or.(type_eq.eq.26)) then
    do j=1,aa%ntr
      fff=real8_inf
      do k1=2,a%n_eq_var
        if (a%eq_var(k1)) then
          val=aa%areaval(k1)%v(j)
          if (val.ne.d0) then
            if (fff==real8_inf) then
              selectcase (type_eq)
              case (23,26)
                fff=intd_m(ia,j,i1,i0,0)
              case default
                fff=intd_m(ia,j,i1,i0,2)
              endselect
            endif
            j1=aa%areaind(j,a%var_ind(k1))
            m_eq(j1)=m_eq(j1)+val*fff
          endif
        endif
      enddo
      !do k1=2,a%n_eq_var
      !  if (a%eq_var(k1)) then
      !    val=aa%areaval(j,k1)
      !    if (val.ne.d0) then
      !      if (type_eq.eq.23) then
      !        fff=intd_m(ia,j,i1,i0,0)
      !      else
      !        fff=intd_m(ia,j,i1,i0,2)
      !      endif
      !      j1=aa%areaind(j,a%var_ind(k1))
      !      m_eq(j1)=m_eq(j1)+val*fff
      !    endif
      !  endif
      !enddo
    enddo
  endif
endif
call update_fict_vars(m_eq,m_b,ia)
call add_equation_to_matrix_thread(ith,ia)
end

subroutine add_equation_to_matrix_thread(ith,ia)
use pgmod
integer(4) ith,ia
type(TOMP_buffer), pointer :: buff
buff=>gs%m%omp_buff%b(ith)
call add_equation_to_matrix(buff%i,buff%eq,buff%b,ia)
end

subroutine add_equation_to_matrix(i,m_eq,m_b,ia)
use pgmod
integer(4) ith,ia,omp_get_thread_num
real(8) m_eq(gs%m%nx)
real(8) m_b,maxm,tmax
type(sparse_matrix_row), pointer :: mrow
type(block_matrix), pointer :: b
real(8), pointer :: m(:)
integer(4), pointer :: c(:)
integer(4) i,j,i1,j2,idamax
type(matrix_mb), pointer :: mb
integer(4), pointer :: buffc(:)
ith=omp_get_thread_num()
buffc=>gs%m%omp_buff%b(ith)%buffc
mb=>gs%a(ia)%m
if (max_equation_coef_eq_1) then
  j=idamax(mb%nx,m_eq(mb%ix),1)
  maxm=dabs(m_eq(j))
  do j=1,mb%ci%ncol
    tmax=dabs(m_eq(mb%ci%col(j)))
    if (tmax>maxm) maxm=tmax
  enddo
  forall (j=mb%ix:mb%ix2) m_eq(j)=m_eq(j)/maxm
  do i1=1,mb%ci%ncol
    j=mb%ci%col(i1)
    m_eq(j)=m_eq(j)/maxm
  enddo
  m_b=m_b/maxm
endif
if (equation_coef_eq_0_eps>d0) then
  do j=mb%ix,mb%ix2
    if (dabs(m_eq(j))<equation_coef_eq_0_eps) m_eq(j)=d0
  enddo
  do i1=1,mb%ci%ncol
    j=mb%ci%col(i1)
    if (dabs(m_eq(j))<equation_coef_eq_0_eps) m_eq(j)=d0
  enddo
endif
select case (gs%m%matrix_store_type)
case (0)
  gs%m%m(i,mb%ix:mb%ix2)=m_eq(mb%ix:mb%ix2)
  do i1=1,mb%ci%ncol
    j=mb%ci%col(i1)
    gs%m%m(i,j)=m_eq(j)
  enddo
case (1)
  mrow=>gs%m%sparse%mrows(i)
  mrow%nc=0
  do j=1,mb%ci%icol_area
    i1=mb%ci%col(j)
    if (m_eq(i1)/=d0) then
      mrow%nc=mrow%nc+1
      buffc(mrow%nc)=i1
    endif
  enddo
  do j=mb%ix,mb%ix2
    if (m_eq(j)/=d0.or.j==i) then
      mrow%nc=mrow%nc+1
      buffc(mrow%nc)=j
    endif
  enddo
  do j=mb%ci%icol_area+1,mb%ci%ncol
    i1=mb%ci%col(j)
    if (m_eq(i1)/=d0) then
      mrow%nc=mrow%nc+1
      buffc(mrow%nc)=i1
    endif
  enddo
  if (use_global_sparsem) then
    if (mrow%nc>mrow%sm_count) call gs_print_stop("Error add_equation_to_matrix (nc>sm_count)")
    j2=mrow%sm_ibegin+mrow%sm_count-1
    m=>gs%m%sparse%m(mrow%sm_ibegin:j2)
    c=>gs%m%sparse%c(mrow%sm_ibegin:j2)
  else
    allocate(mrow%m(mrow%nc),mrow%c(mrow%nc))
    m=>mrow%m
    c=>mrow%c
  endif
  c(1:mrow%nc)=buffc(1:mrow%nc)
  do j=1,mrow%nc
    m(j)=m_eq(c(j))
  enddo
case (2)
  b=>gs%m%bm%blockm(ia)
  i1=i-b%iu+1 !локальный индекс уравнения в блоке
  b%m(i1,:)=m_eq(b%iu:b%iu2)
  mrow=>b%sparse%mrows(i1)
  do j=1,mb%ci%ncol
    j2=mb%ci%col(j)
    if (m_eq(j2)/=d0) then
      mrow%nc=mrow%nc+1
      buffc(mrow%nc)=j2
    endif
  enddo
  if (mrow%nc>0) then
    allocate(mrow%m(mrow%nc),mrow%c(mrow%nc))
    mrow%c=buffc(1:mrow%nc)
    do j=1,mrow%nc
      mrow%m(j)=m_eq(mrow%c(j))
    enddo
  endif
end select
gs%m%b(i)=m_b
!обнуляем за собой
m_eq(mb%ix:mb%ix2)=d0
do j=1,mb%ci%ncol
  m_eq(mb%ci%col(j))=d0
enddo
end

subroutine update_fict_vars(m_eq,m_b,ia)
!учет фиктивных неизвестных
use pgmod
real(8) m_eq(gs%m%nx_all),m_b
real(8) m_eqj,cc0,get_gug_c
integer(4) j,j1,j2,ia,i
type(fict_var), pointer :: fv
type(TBoungline_GU_gen), pointer :: gug
type(TBoungline_GU_genGlobal), pointer :: guglob
type(matrix_mb), pointer :: mb
mb=>gs%a(ia)%m
do i=1,mb%ci%ncol_fict
  j=mb%ci%col_fict(i)
  m_eqj=m_eq(j)
  if (m_eqj==d0) cycle
  fv=>gs%m%fvar(j-gs%m%ix_fict+1)
  selectcase (fv%mode)
  case (1)
    gug=>fv%gu_gen
    cc0=get_gug_c(gug%carr(0),fv%i)
    if (cc0.ne.d0) m_b=m_b-cc0*m_eqj
    do j1=1,gug%n
      j2=fv%psiind(j1)
      cc0=get_gug_c(gug%carr(j1),fv%i)
      if (j2==0) then
        m_b=m_b-cc0*m_eqj*fv%psiom(j1)
      elseif (j2>=gs%m%ix_fict) then
        call gs_print_stop("Error update_fict_vars. Fictive var depend from Fictive var")
        !!!переменные, стоящие слева (gu%bndf) - являются фиктивными переменными
        !!!переменные, стоящие справа (gu%bndg%bndf) - являются реальными переменными
      else
        m_eq(j2)=m_eq(j2)+cc0*m_eqj
      endif
    enddo
  case (2)
    guglob=>fv%gu_genGlob
    cc0=get_gug_c(guglob%carr(0),fv%i)
    if (cc0.ne.d0) m_b=m_b-cc0*m_eqj
    do j1=1,guglob%n
      j2=fv%psiind(j1)
      cc0=get_gug_c(guglob%carr(j1),fv%i)
      if (j2==0) then
        call gs_print_stop("Error update_fict_vars. Global_ind=0")
      elseif (j2>=gs%m%ix_fict) then
        call gs_print_stop("Error update_fict_vars. Fictive var depend from Fictive var")
        !!!переменные, стоящие слева (gu%bndf) - являются фиктивными переменными
        !!!переменные, стоящие справа (gu%bndg%bndf) - являются реальными переменными
      else
        m_eq(j2)=m_eq(j2)+cc0*m_eqj
      endif
    enddo
  endselect
enddo
do j=1,mb%ci%ncol_fict
  m_eq(mb%ci%col_fict(j))=d0
enddo
end

subroutine get_matrix_area_ap(ia,icp,func_inf)
!построение матрицы СЛАУ для площадных элементов
use pgmod
integer(4) i,j,k,i1,j1,i2,j0,i3,icp,ith,omp_get_thread_num
integer(4) ia !номер области
real(8) intd_ar,calcint_ar,calcintn_ar,fff
real(8) val,val2,val3
real(8) intd_oss_ar
real(8) func_inf
external func_inf
type(TArea), pointer :: a
type(areatype), pointer :: aa
type(TBound), pointer :: b
integer(4) nf22_1(6),nf22_2(6),nf22_3(6),nf22_4(6),nf22_5(6),nf26_1(3),nf26_2(3),nf26_3(3)
data nf22_1 /1,1,6,7,6,7/
data nf22_2 /0,0,4,5,4,5/
data nf22_3 /3,-1,10,11,-1,-1/
data nf22_4 /2,-1,8,9,-1,-1/
data nf22_5 /2,0,8,9,4,5/
data nf26_1 /1,6,7/
data nf26_2 /0,4,5/
data nf26_3 /0,4,5/
logical type_eq22,areaEq22
integer(4) i_eq
real(8), pointer :: m_eq(:)
real(8), pointer :: m_b
type(TOMP_buffer), pointer :: buff
type(TAreaValue), pointer :: av,av2,av3
ith=omp_get_thread_num()
buff=>gs%m%omp_buff%b(ith)
m_eq=>buff%eq
m_b=>buff%b

a=>gs%a(ia)
aa=>a%a
i=aa%cppa%j(icp)
buff%i=i
i1=aa%cppa%itr(icp)
i2=aa%cppa%iu(icp)
i_eq=0
areaEq22=aa%areaEq.eq.22.or.aa%areaEq.eq.24.or.aa%areaEq.eq.25
if (areaEq22.or.aa%areaEq==26) i_eq=aa%eq_ind(i2)
!m_eq=d0 !перенесли в add_equation_to_matrix_thread
m_b=d0
j1=aa%areaind(i1,i2)
m_eq(j1)=m_eq(j1)-pi2
k=0
type_eq22=(areaEq22).and.((i_eq==2).or.(i_eq==5).or.(i_eq==6))
if ((aa%areaEq.eq.2).or.(aa%areaEq.eq.7).or.(aa%areaEq.eq.12).or.(aa%areaEq.eq.13).or.(aa%areaEq.eq.14).or.(type_eq22)) k=2
do j0=1,a%nb
  b=>a%bnd(j0)
  do j=1,b%npanel
    j1=b%psiind(j,1+k)
    val=b%psiom(j,1+k)
    if (j1.ne.0.or.val.ne.d0) then
      select case (aa%areaEq)
      case (12)
        fff=calcintn_ar(j,i1,13,ia,j0) !calcintn(j,xmc(i1,ia),ymc(i1,ia),13,ia,j0)
      case (13)
	      fff=calcintn_ar(j,i1,15,ia,j0) !calcintn(j,xmc(i1,ia),ymc(i1,ia),15,ia,j0)
      case (16,17)
  	    fff=calcint_ar(j,i1,5+i2,ia,j0) !calcint(j,xmc(i1,ia),ymc(i1,ia),5+i2,ia,j0) 
      case (18:20)
        fff=calcint_ar(j,i1,7,ia,j0) !calcint(j,xmc(i1,ia),ymc(i1,ia),7,ia,j0) 
      case (22,24,25)
        fff=calcint_ar(j,i1,nf22_1(i_eq),ia,j0)
      case (26)
        fff=calcint_ar(j,i1,nf26_1(i_eq),ia,j0)
      case default
        fff=calcint_ar(j,i1,1,ia,j0) !calcint(j,xmc(i1,ia),ymc(i1,ia),1,ia,j0) 
      end select
      if (j1.eq.0) then
        m_b=m_b-val*fff
      else
        m_eq(j1)=m_eq(j1)+fff
      endif
    endif
    j1=b%psiind(j,2+k)
    val=b%psiom(j,2+k)
    if (j1.ne.0.or.val.ne.d0) then
      select case (aa%areaEq)
      case (12)
        fff=calcintn_ar(j,i1,12,ia,j0) !calcintn(j,xmc(i1,ia),ymc(i1,ia),12,ia,j0)
      case (13)
        fff=calcintn_ar(j,i1,14,ia,j0) !calcintn(j,xmc(i1,ia),ymc(i1,ia),14,ia,j0)
      case (16,17)
   	    fff=calcint_ar(j,i1,3+i2,ia,j0) !calcint(j,xmc(i1,ia),ymc(i1,ia),3+i2,ia,j0) 
      case (18:20)
        fff=calcint_ar(j,i1,5,ia,j0) !calcint(j,xmc(i1,ia),ymc(i1,ia),5,ia,j0) 
      case (22,24,25)
        fff=calcint_ar(j,i1,nf22_2(i_eq),ia,j0)
      case (26)
        fff=calcint_ar(j,i1,nf26_2(i_eq),ia,j0)
      case default
        fff=calcint_ar(j,i1,0,ia,j0) !calcint(j,xmc(i1,ia),ymc(i1,ia),0,ia,j0) 
      end select
      if (j1.eq.0) then
        m_b=m_b+val*fff
      else
        m_eq(j1)=m_eq(j1)-fff
      endif
    endif
  enddo
enddo
if ((areaEq22).and.((i_eq==1).or.(i_eq==3).or.(i_eq==4))) then
  do j0=1,a%nb
    b=>a%bnd(j0)
    do j=1,b%npanel
      j1=b%psiind(j,3)
      val=b%psiom(j,3)
      if (j1.ne.0.or.val.ne.d0) then
        fff=calcint_ar(j,i1,nf22_3(i_eq),ia,j0)
        if (j1.eq.0) then
          m_b=m_b-val*fff
        else
          m_eq(j1)=m_eq(j1)+fff
        endif
      endif
      j1=b%psiind(j,4)
      val=b%psiom(j,4)
      if (j1.ne.0.or.val.ne.d0) then
        fff=calcint_ar(j,i1,nf22_4(i_eq),ia,j0)
        if (j1.eq.0) then
          m_b=m_b+val*fff
        else
          m_eq(j1)=m_eq(j1)-fff
        endif
      endif
    enddo
  enddo
endif
select case (aa%areaEq)
case (5,6,14)
  av=>aa%areaval(2)
  do j=1,aa%ntr
    val=av%v(j)
    if (val.ne.d0) then
      j1=aa%areaind(j,1)
      m_eq(j1)=m_eq(j1)-val*intd_ar(ia,j,i1,0) !intd_an(aa%zmc(i1),zmz(1,j,ia),zmz(2,j,ia),zmz(3,j,ia),0)
    endif
  enddo
end select
select case (aa%areaEq)
case (6,7)
  av=>aa%areaval(1)
  do j=1,aa%ntr
    val=av%v(j)
    if (val.ne.d0) m_b=m_b-val*intd_ar(ia,j,i1,0) !intd_an(aa%zmc(i1),zmz(1,j,ia),zmz(2,j,ia),zmz(3,j,ia),0)
  enddo
end select
!уравнение переноса и осесимметричные уравнения
select case (aa%areaEq)
case (16:20)
  av=>null()
  if ((aa%areaEq.eq.17).or.(aa%areaEq.eq.18)) av=>aa%areaval(1)
  av2=>aa%areaval(4)
  av3=>null()
  if ((aa%areaEq.eq.16).or.(aa%areaEq.eq.17)) av3=>aa%areaval(5)
  do j=1,aa%ntr
    val=d0
    val3=d0
    if ((aa%areaEq.eq.17).or.(aa%areaEq.eq.18)) val=av%v(j)
    val2=av2%v(j)
    if ((aa%areaEq.eq.16).or.(aa%areaEq.eq.17)) val3=av3%v(j)
    if ((val.ne.d0).or.(val2.ne.d0).or.(val3.ne.d0)) then
      if (aa%areaEq.eq.20) then
        fff=intd_oss_ar(ia,j,i1,5) !intd_oss(aa%zmc(i1),zmz(1,j,ia),zmz(2,j,ia),zmz(3,j,ia),5)
      else
        i3=i2
        if ((aa%areaEq.eq.18).or.(aa%areaEq.eq.19)) i3=2
        fff=intd_ar(ia,j,i1,i3+3) !intd_an(aa%zmc(i1),zmz(1,j,ia),zmz(2,j,ia),zmz(3,j,ia),i3+3)
      endif
      if (val.ne.0) m_b=m_b-val*fff
      if (val2.ne.0) then
        j1=aa%areaind(j,1)
        m_eq(j1)=m_eq(j1)+val2*fff
      endif
      if (val3.ne.0) then
        j1=aa%areaind(j,2)
        m_eq(j1)=m_eq(j1)+val3*fff
      endif
    endif
  enddo
end select
if (areaEq22.or.aa%areaEq==26) then
  do j=1,aa%ntr
    fff=real8_inf
    do j0=1,a%n_eq_var
      if (a%eq_var(j0)) then
        val=aa%areaval(j0)%v(j)
        if (val.ne.d0) then
          if (fff==real8_inf) then
            if (areaEq22) then
              fff=intd_ar(ia,j,i1,nf22_5(i_eq))
            else
              fff=intd_ar(ia,j,i1,nf26_3(i_eq))
            endif
          endif
          if (j0==1) then
            m_b=m_b-val*fff
          else
            j1=aa%areaind(j,a%var_ind(j0))
            m_eq(j1)=m_eq(j1)+val*fff
          endif
        endif
      endif
    enddo
    !do j0=1,a%n_eq_var
    !  val=aa%areaval(j,j0)
    !  if (a%eq_var(j0).and.val.ne.d0) then
    !    fff=intd_ar(ia,j,i1,nf22_5(i_eq))
    !    if (j0==1) then
    !      m_b=m_b-val*fff
    !    else
    !      j1=aa%areaind(j,a%var_ind(j0))
    !      m_eq(j1)=m_eq(j1)+val*fff
    !    endif
    !  endif
    !enddo
  enddo
endif
call update_fict_vars(m_eq,m_b,ia)
call add_equation_to_matrix_thread(ith,ia)
end

subroutine main_sparse
use pgmod
integer(4) i,j,j1,k,i1
type(sparse_matrix_row), pointer :: mrow
gs%m%sparse%nc=0
do i=1,gs%m%nu
  gs%m%sparse%nc=gs%m%sparse%nc+gs%m%sparse%mrows(i)%nc
enddo
if (.not.use_global_sparsem) then
  allocate(gs%m%sparse%m(gs%m%sparse%nc))
  allocate(gs%m%sparse%c(gs%m%sparse%nc))
endif
allocate(gs%m%sparse%r(gs%m%nu+1))
j=1
do i=1,gs%m%nu
  mrow=>gs%m%sparse%mrows(i)
  j1=j+mrow%nc-1
  if (use_global_sparsem) then
    k=mrow%sm_ibegin
    if (j<k) then
      do i1=j,j1
        gs%m%sparse%m(i1)=gs%m%sparse%m(k)
        gs%m%sparse%c(i1)=gs%m%sparse%c(k)
        k=k+1
      enddo
    endif
  else
    gs%m%sparse%m(j:j1)=mrow%m(1:mrow%nc)
    gs%m%sparse%c(j:j1)=mrow%c(1:mrow%nc)
  endif
  gs%m%sparse%r(i)=j
  j=j1+1
  call bm_null_mrow(mrow)
enddo
gs%m%sparse%r(gs%m%nu+1)=j
if (use_global_sparsem.and.j<=gs%m%nc_global) then
  gs%m%sparse%m(j:gs%m%nc_global)=d0
  gs%m%sparse%c(j:gs%m%nc_global)=0
endif
call sparse_memory(1)
end

subroutine main_block_sparse
use pgmod
integer(4) ia,i,j
integer(8) nc
type(block_matrix), pointer :: b
type(sparse_matrix_row), pointer :: mrow
gs%m%memory_mb=0
nc=0
do ia=1,gs%m%bm%nblock
  b=>gs%m%bm%blockm(ia)
  j=0
  do i=1,b%nu
    mrow=>b%sparse%mrows(i)
    j=j+mrow%nc
  enddo
  nc=nc+b%nx*b%nu+j
  gs%m%memory_mb=gs%m%memory_mb+b%nx*b%nu*8+j*12+(b%nu+1)*4
enddo
gs%m%sparse_percent=nc*100.0/gs%m%nx/gs%m%nu
gs%m%memory_mb=gs%m%memory_mb/1024**2
call bm_main_block_sparse(gs%m%bm)
end

subroutine pg_solve
!dec$ attributes dllexport:: pg_solve
!решение СЛАУ
use pgmod
integer(4), parameter :: mode=1 !0 - IMSL
                                !1 - MKL  
integer(4), parameter :: mode_sp=1 !0 - MKL DSS
                                   !1 - MKL PARDISO
                                   !2 - amgcl
logical write_by_timef2
external write_by_timef2
call pg_solve_begin
call pg_solve_iter
call pg_solve_end(.true.)
end

subroutine pg_solve_iter
!dec$ attributes dllexport:: pg_solve_iter
!решение СЛАУ с запоминанием 
use pgmod
integer(4), parameter :: mode=1 !0 - IMSL
                                !1 - MKL  
integer(4), parameter :: mode_sp=1 !0 - MKL DSS
                                   !1 - MKL PARDISO
                                   !2 - amgcl
logical write_by_timef2
external write_by_timef2

call test_closing_slau
if (gs%m%matrix_store_type==1) call main_sparse
if (gs%m%matrix_store_type==2) call main_block_sparse
call gs_print('System solving')
call gs_print('neq='//itoa(gs%m%nu))
!call test_matrix
!call update_matr
!if (gs%m%matrix_store_type==1) call sparse_to_m
!if (gs%i==2.and.gs%m%matrix_store_type==2) then
!  call block_diagonal_to_m
!  gs%m%matrix_store_type=0
!endif
!call test_matr_maxmincol
!call matrix_resort
!call m_to_sparse
!gs%m%matrix_store_type=0
!call test_m_diagonal
!call drw_matrix
!call write_matr_b
!call write_matrix_mm


!gs%m%find_norm=.true.

if (gs_write_matrix) call write_matr
call init_gs_time
gs_use_mkl_progress=.true.
gs_system_count=1
gs_system_i=1
select case (gs%m%matrix_store_type)
case (0)
  if (gs%m%nx/=gs%m%nu) then
    call calc_slau_ls(gs%m%m, gs%m%b, gs%m%nu, gs%m%nx, gs%m%nu, gs%m%res, mode,gs%m%find_norm,gs%m%norm)
  else
    call calc_slau(gs%m%m, gs%m%b, gs%m%nu, gs%m%nx, gs%m%res, mode, gs%m%find_norm,gs%m%norm)
  endif
case (1)
  call calc_slau_sparse(gs%m%sparse%m, gs%m%sparse%c, gs%m%sparse%r, gs%m%b, gs%m%nu, gs%m%nx, gs%m%sparse%nc, gs%m%res, mode_sp, gs%m%find_norm,gs%m%norm, gs_pardiso_eps)
case (2)
  gs_system_count=gs%m%bm%nblock
  !call bm_BiCGStab(gs%m%bm,gs%m%b,gs%m%res,gs_block_matrix_solver_maxIter,gs_block_matrix_solver_eps,gs_block_matrix_solver_maxIter_convergence,write_by_timef2)
  call bm_calc_slau_block_diagonal2(gs%m%bm,gs%m%b,gs%m%res,gs_block_matrix_solver_lambda,gs_block_matrix_solver_maxIter,gs_block_matrix_solver_eps,gs_block_matrix_solver_maxIter_convergence,gs_system_i,write_by_timef2)
end select
gs_use_mkl_progress=.false.
call print_total_time
if (gs_write_matrix) call write_res(gs%m%res,gs%m%nx)
end

subroutine pg_solve_begin
use pgmod
if (allocated(gs%m%res)) deallocate(gs%m%res)
allocate(gs%m%res(gs%m%nx))
if(gs%m%have_initial_approxsolv) then
  call get_set_res(.false.)
else
  gs%m%res=d0
endif
end

subroutine pg_solve_end(free_res)
use pgmod
logical free_res
call get_set_res(.true.)
call solve_postprocessor
if (free_res) then
  deallocate(gs%m%res)
endif
end

!subroutine init_lam2
!use pgmod
!type(TBound), pointer :: b 
!type(TArea), pointer :: a
!integer(4) i,j,k,i0,j0
!gs%m%is_lam2=.false.
!if (gs%i==1) return
!do i0=1,gs%na
!  a=>gs%a(i0)
!  do j0=1,a%nb
!    b=>a%bnd(j0) 
!	  if (b%boundLineType==1) then
!      do i=1,b%npanel
!        do j=1,a%umax
!          k=b%psiind(i,j)
!          if (k.ne.0.and.k<gs%m%ix_fict) then
!            gs%m%is_lam2(k)=mod(j,2)==0
!          endif
!        enddo
!      enddo
!    endif
!  enddo
!enddo
!end

subroutine get_set_res(is_get)
use pgmod
type(TBound), pointer :: b 
type(TArea), pointer :: a
type(TBoundline2), pointer :: bl2
type(TBounLineFuncApprox), pointer :: ga
integer(4) i,j,k,i0,j0,j1,k1
type(TBoungline_GU_gen), pointer :: gug
type(TBoungline_GU_genGlobal), pointer :: guglob
type(fict_var), pointer :: fv
logical is_get
integer(4), allocatable :: nn(:)
integer(4) nn0
real(8) cc0,get_gug_c
if (.not.is_get) then
  allocate(nn(gs%m%nx))
  nn=0
  nn0=0
  gs%m%res=d0
endif
do i0=1,gs%na
  a=>gs%a(i0)
  do j0=1,a%nb
    b=>a%bnd(j0) 
	  if (b%boundLineType==1) then
      do i=1,b%npanel
        do j=1,a%umax
          k=b%psiind(i,j)
          if (k.ne.0.and.k<gs%m%ix_fict) then
            if (is_get) then
              b%psiom(i,j)=gs%m%res(k)
            else
              gs%m%res(k)=gs%m%res(k)+b%psiom(i,j)
              nn(k)=nn(k)+1
            endif
          endif
        enddo
      enddo
	  elseif (b%boundLineType==2.and.is_get) then
	    do i=1,b%nline
	      bl2=>b%line2(i)
	  	  do j=1,a%umax
	  	    ga=>bl2%ga(j)
	  	    do j1=1,ga%n
	  	      k=ga%ind(j1)
            if (k.ne.0) then
              ga%val(j1)=gs%m%res(k)
            endif
	  	    enddo
	  	  enddo
	    enddo
	  endif
  enddo
enddo
if (is_get) then
  !переменные в области
  do i0=1,gs%na
    a=>gs%a(i0)
    if(allocated(a%a%areaind)) then
      do i=1,a%a%ntr
  	    do j=1,a%umaxtr
          k=a%a%areaind(i,j)
          if (k.ne.0) a%a%psiarea(i,j)=gs%m%res(k)
        enddo
      enddo
    endif
  enddo
  !фиктивные переменные на границе
  do i0=1,gs%na
    a=>gs%a(i0)
    do j0=1,a%nb
      b=>a%bnd(j0) 
  	  if (b%boundLineType==1) then
        do i=1,b%npanel
          do j=1,a%umax
            k=b%psiind(i,j)
            if (k.ne.0.and.k>=gs%m%ix_fict) then
              fv=>gs%m%fvar(k-gs%m%ix_fict+1)
              selectcase (fv%mode)
              case (1)
                gug=>fv%gu_gen
                b%psiom(i,j)=get_gug_c(gug%carr(0),fv%i)
                do j1=1,gug%n
                  k1=fv%psiind(j1)                
                  cc0=get_gug_c(gug%carr(j1),fv%i)
                  if (k1==0) then
                    b%psiom(i,j)=b%psiom(i,j)+cc0*fv%psiom(j1)
                  else
                    b%psiom(i,j)=b%psiom(i,j)+cc0*gs%m%res(k1)
                  endif
                enddo
              case (2)
                guglob=>fv%gu_genGlob
                b%psiom(i,j)=get_gug_c(guglob%carr(0),fv%i)
                do j1=1,guglob%n
                  k1=fv%psiind(j1)                
                  cc0=get_gug_c(guglob%carr(j1),fv%i)
                  if (k1==0) then
                    call gs_print_stop("Error get_set_res. Global_ind=0")
                  else
                    b%psiom(i,j)=b%psiom(i,j)+cc0*gs%m%res(k1)
                  endif
                enddo
              endselect
            endif
          enddo
        enddo
  	  endif
    enddo
  enddo
  do i=1,gs%constvaln
    gs%constvals(i)=gs%m%res(gs%constvalinds(i))
  enddo
else
  do i=1,gs%m%nx
    if (nn(i)>1) then
      gs%m%res(i)=gs%m%res(i)/nn(i)
    endif
  enddo
  deallocate(nn)
endif
end

subroutine pg_set_initial_approxsolv(f_ias)
!dec$ attributes dllexport:: pg_set_initial_approxsolv
!задать начальное приближение к решению
!начальное приближение задается только на границе - необходимо для случая декомпозиции области - блочно-диагонаьная матрица
use pgmod
real(8) f_ias !(x,y,nf,der,ia)
              !real(8) x,y    - координаты
              !integer(4) nf  - номер функции 1..4
              !integer(4) der - для производных 1 - d/dx, 2 - d/dy
              !integer(4) ia  - область
external f_ias
integer(4) i0,j0,i,j,k
type(TBound), pointer :: b 
type(TArea), pointer :: a
logical need_ias(4),need_der,is_der(4),need_der_x,need_der_y
real(8) dfx,dfy,scal_p
real(8), parameter :: eps=1.0d-8
complex(8) nn
is_der=.false.
is_der(2)=.true.
is_der(4)=.true.
do i0=1,gs%na
  a=>gs%a(i0)
  do j0=1,a%nb
    b=>a%bnd(j0) 
	  if (b%boundLineType==1) then
      do i=1,b%npanel
        need_der=.false.
        do j=1,a%umax
          k=b%psiind(i,j)
          need_ias(j)=k.ne.0.and.k<gs%m%ix_fict
          !need_ias(j)=.true.
          need_der=need_der.or.(is_der(j).and.need_ias(j))
        enddo
        if (need_der) then
          nn=-b%ett(i)*ii
          need_der_x=dabs(dreal(nn))>eps
          need_der_y=dabs(dimag(nn))>eps
        endif
        do j=1,a%umax
          if (need_ias(j)) then
            if (is_der(j)) then
              dfx=d0
              dfy=d0
              if (need_der_x) dfx=f_ias(b%xc(i),b%yc(i),j,1,i0)
              if (need_der_y) dfy=f_ias(b%xc(i),b%yc(i),j,2,i0)
              b%psiom(i,j)=scal_p(nn,dcmplx(dfx,dfy))
            else
              b%psiom(i,j)=f_ias(b%xc(i),b%yc(i),j,0,i0)
            endif
          endif
        enddo
      enddo
	  endif
  enddo
enddo
gs%m%have_initial_approxsolv=.true.
end


subroutine write_res(res,n)
use pgmod
integer(4) n,i
real(8) res(n)
OPEN (1,FILE='res.dat')
do i=1,n
  write(1,*) res(i)
enddo
close(1)
end

subroutine write_matr
use pgmod
integer(4) j,i
OPEN (1,FILE='matr.dat')
do i=1,gs%m%nu
  write(1,*) i
  do j=1,gs%m%nx
    write(1,*) j, gs%m%m(i,j)
  enddo
enddo
write(1,*) 'b'
do j=1,gs%m%nu
  write(1,*) j, gs%m%b(j)
enddo
close(1)
end

subroutine write_matr_b
use pgmod
integer(4) j
OPEN (1,FILE='b.dat')
do j=1,gs%m%nu
  write(1,*) j, gs%m%b(j)
enddo
close(1)
end

subroutine test_matrix
!разделение коэффициентов матрицы по порядку близости значений к нулю
use pgmod
integer(4), parameter :: n=25
real(8) maxm,t,eps,mmmr(n+1)
integer(4) i,j,k,mmm(n+1)
!mmm(i) - количество элементов в диапазоне [10^(-i+1);10^(-i)]
!mmm(n+1) - чистые нули
mmm=0
if (max_equation_coef_eq_1) then
  maxm=d1
else
  maxm=d0
  do i=1,gs%m%nu
    do j=1,gs%m%nx
      if (dabs(gs%m%m(i,j))>maxm) maxm=dabs(gs%m%m(i,j))
    enddo  
  enddo  
endif
eps=dexp(-n*dlog(10.0d0))
do i=1,gs%m%nu
  do j=1,gs%m%nx
    if (gs%m%m(i,j)==d0) then
      mmm(n+1)=mmm(n+1)+1
      cycle
    endif
    t=dabs(gs%m%m(i,j))/maxm
    if (t<eps) then
      k=n
    else
      t=dlog10(t)
      k=1-int(t)
    endif
    if (k<1) then
      k=1
    elseif (k>n) then
      k=n
    endif
    mmm(k)=mmm(k)+1
  enddo  
enddo  
k=gs%m%nu*gs%m%nx
forall (i=1:n+1) mmmr(i)=(mmm(i)+d0)/k
end

subroutine update_matr
!замена малых коэффициентов матриы на ноль
use pgmod
real(8) maxm,eps,t,t1
integer(4) i,j,k,k1
maxm=d0
do i=1,gs%m%nu
  do j=1,gs%m%nx
    if (dabs(gs%m%m(i,j))>maxm) maxm=dabs(gs%m%m(i,j))
  enddo  
enddo  
eps=1.0d-12
k=0
k1=0
do i=1,gs%m%nu
  do j=1,gs%m%nx
    t=dabs(gs%m%m(i,j))/maxm
    if (t<eps) then
      if (gs%m%m(i,j)==d0) k1=k1+1
      gs%m%m(i,j)=d0
      k=k+1
    endif
  enddo  
enddo  
t=(k+d0)/(gs%m%nu*gs%m%nx)
t1=(k1+d0)/(gs%m%nu*gs%m%nx)
end

integer function mkl_progress( thread, step, stage )
use pgmod
 integer*4 thread, step, i
 character*(*) stage
 character(75) str
 mkl_progress = 0
 if (.not.gs_use_mkl_progress) return
 select case (gs%m%matrix_store_type)
 case (0)
   i=step*100/gs%m%nu
 case (1)
   i=step
 case (2)
   i=step*100/gs%m%bm%blockm(gs_system_i)%nu
   i=(i+(gs_system_i-1)*100)/gs_system_count
 end select
 write(str,"('Solve system(', a6, '): ', I0, '%')") stage,i
 call write_by_time(trim(str))
 i=thread !чтобы не было предупреждения при компиляции
 return
end

subroutine sparse_to_m
use pgmod
integer(4) i,j
if (.not.allocated(gs%m%m)) allocate(gs%m%m(gs%m%nu,gs%m%nx))
gs%m%m=d0
do i=1,gs%m%nu
  do j=gs%m%sparse%r(i),gs%m%sparse%r(i+1)-1
    gs%m%m(i,gs%m%sparse%c(j))=gs%m%sparse%m(j)
  enddo
enddo
end

subroutine m_to_sparse
use pgmod
integer(4) i,j,k
gs%m%sparse%nc=0
do i=1,gs%m%nu
  do j=1,gs%m%nx
    if (gs%m%m(i,j)/=d0.or.i==j) gs%m%sparse%nc=gs%m%sparse%nc+1
  enddo
enddo
if (.not.allocated(gs%m%sparse%m)) then  
  allocate(gs%m%sparse%m(gs%m%sparse%nc))
  allocate(gs%m%sparse%c(gs%m%sparse%nc))
  allocate(gs%m%sparse%r(gs%m%nu+1))
endif
k=1
do i=1,gs%m%nu
  gs%m%sparse%r(i)=k
  do j=1,gs%m%nx
    if (gs%m%m(i,j)/=d0.or.i==j) then
      gs%m%sparse%m(k)=gs%m%m(i,j)
      gs%m%sparse%c(k)=j
      k=k+1
    endif
  enddo
enddo
gs%m%sparse%r(gs%m%nu+1)=k
end

subroutine block_diagonal_to_m
use pgmod
integer(4) i,j,ia,i1
type(block_matrix), pointer :: b
if (.not.allocated(gs%m%m)) allocate(gs%m%m(gs%m%nu,gs%m%nx))
gs%m%m=d0
do ia=1,gs%m%bm%nblock
  b=>gs%m%bm%blockm(ia)
  gs%m%m(b%iu:b%iu2,b%ix:b%ix2)=b%m
  do i=1,b%nu
    i1=b%iu+i-1
    do j=b%sparse%r(i),b%sparse%r(i+1)-1
      gs%m%m(i1,b%sparse%c(j))=b%sparse%m(j)
    enddo
  enddo
enddo
end

subroutine drw_matrix
use pgmod
integer(4) i,j,k
k=0
do i=1,gs%m%nu
  do j=1,gs%m%nx
    if (gs%m%m(i,j)/=d0) k=k+1
  enddo
enddo
OPEN (1,FILE='matrp.dat')
write(1,*) 'TITLE = "m"'
write(1,*) 'VARIABLES = "j", "i", "m"'
write(1,"('ZONE T=""m0"", I=', i0, ', F=POINT')") gs%m%nu*gs%m%nx-k
do i=1,gs%m%nu
  do j=1,gs%m%nx
    if (gs%m%m(i,j)==d0) WRITE(1,"(i9,i9,E13.6)") j,i,dabs(gs%m%m(i,j))
  enddo
enddo
write(1,"('ZONE T=""m1"", I=', i0, ', F=POINT')") k
do i=1,gs%m%nu
  do j=1,gs%m%nx
    if (gs%m%m(i,j)/=d0) WRITE(1,"(i0,' ',i0,' ',E13.6)") j,i,dabs(gs%m%m(i,j))
  enddo
enddo
close(1)
end

subroutine test_m_diagonal
use pgmod
integer(4) i,j
j=0
do i=1,gs%m%nx
  if (gs%m%m(i,i)/=d0) j=j+1
enddo
write(*,"('Not 0:'i0,' from ',i0)") j,gs%m%nx
end

subroutine matrix_resort
use pgmod
integer(4) i,j,k,nb,mc(gs%m%nu),minc,cmne0count
logical badu(gs%m%nu)
real(8) tm(gs%m%nx),tb
badu=.false.
nb=0
do i=1,gs%m%nu
  !ищем уравнение с нулевым диагональным элементом
  if (gs%m%m(i,i)/=d0) cycle
  k=0
  !ищем уравнение, где этот элемент не нулевой
  do j=1,gs%m%nx
    if (gs%m%m(i,j)/=d0.and.gs%m%m(j,i)/=d0) then
      k=j
      exit
    endif
  enddo
  if (k==0) then
    !если такое не найдено, то помечаем текущее уравнение как плохое
    badu(i)=.true.
    nb=nb+1
  else
    !если найдено - меняем местами
    tm=gs%m%m(k,:)
    gs%m%m(k,:)=gs%m%m(i,:)
    gs%m%m(i,:)=tm
    tb=gs%m%b(k)
    gs%m%b(k)=gs%m%b(i)
    gs%m%b(i)=tb
  endif
enddo
!обработка плохих уравнений
if (nb>0) then
  !количество ненулевых коэффициентов в уравнениях
  do i=1,gs%m%nu
    mc(i)=cmne0count(i)
  enddo
  do i=1,gs%m%nu
    if (.not.badu(i)) cycle
    !ищем уравнение, где нужный диагональный элемент не нулевой, а количество коэф минимальное
    k=0
    minc=gs%m%nx+1
    do j=1,gs%m%nu
      if (mc(j)<minc.and.gs%m%m(j,i)/=d0) then
        minc=mc(j)
        k=j
      endif
    enddo
    !добавляем это уравнение к текущему
    gs%m%m(i,:)=gs%m%m(i,:)+gs%m%m(k,:)
    gs%m%b(i)=gs%m%b(i)+gs%m%b(k)
    mc(i)=cmne0count(i)
    badu(i)=.false.
  enddo
endif
end

function cmne0count(i)
use pgmod
integer(4) cmne0count,i,j,k
k=0
do j=1,gs%m%nx
  if (gs%m%m(i,j)/=d0) k=k+1
enddo
cmne0count=k
end

subroutine test_matr_maxmincol
use pgmod
integer(4) i,j,mc(gs%m%nx),k
!logical eq_eps
mc=0
k=0
do i=1,gs%m%nx
  do j=1,gs%m%nu
    !количество 1
    !if (eq_eps(d1,dabs(gs%m%m(j,i)),1.0d-8)) then
    !  mc(i)=mc(i)+1
    !  k=k+1
    !endif
    !количество ненулевых элементов
    if (gs%m%m(j,i)/=d0) then
      mc(i)=mc(i)+1
      k=k+1
    endif
  enddo
enddo
end

subroutine write_matrix_mm
!сохранение маткивы в формате MatrixMarket
!для работы этой функции надо раскомментарить call mmwrite
!и добавить к проекту файл common_fort/mmio.f
use pgmod
integer(4) ival(1),i,j,k
complex cval(1)
integer(4), allocatable :: indx(:)
allocate(indx(gs%m%sparse%nc))
ival=0
cval=cmplx(0.0,0.0)
k=0
do i=1,gs%m%nu
  do j=gs%m%sparse%r(i),gs%m%sparse%r(i+1)-1
    k=k+1
    indx(k)=i
  enddo
enddo
OPEN (1,FILE='m.dat')
!call mmwrite(1,'coordinate','real','general',gs%m%nu,gs%m%nx,gs%m%sparse%nc,indx,gs%m%sparse%c,ival,gs%m%sparse%m,cval)
close(1)
OPEN (1,FILE='b.dat')
do i=1,gs%m%nu
  write(1,*) gs%m%b(i)
enddo
close(1)
deallocate(indx)
end

subroutine test_closing_slau
use pgmod
integer(4) i,nnt
nnt=0
do i=1,gs%na
  nnt=nnt+gs%a(i)%m%nnt
enddo
if (nnt/=gs%m%nu) call gs_print_stop("Error closing matrix nnt/=nu!!! (test_closing_slau)")
end
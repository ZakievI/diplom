module pgmod
    use gen_mod
    use func_mod
    use slau_block

    real(8), parameter :: ds0numint=1.0d-3
    integer(4), parameter :: max_int=19 !максимальное количество типов интегралов для кэша
    !real(8), parameter :: real8_inf=1.0d0/0.0d0
    real(8), parameter :: real8_inf=Z'0FFFFFFFFFFFFFFF'
    integer(4), parameter :: int4_empty=Z'0FFFFFFF'
    integer(4), parameter :: gs_log_file_i=999
    integer(4), parameter :: int_type(0:max_int) = [0,0,0,0,1,2,1,2,1,2,1,2,0,0,0,0,1,2,1,2] !0 - обычный интеграл, 1 - d/dx, 2 - d/dy
    integer(4), parameter :: int_type_dn(0:max_int) = [0,1,0,1,0,0,1,1,0,0,1,1,0,1,0,1,0,0,1,1] !0 - обычный интеграл, 1 - интеграл с d/dn
    integer(4), parameter :: int_type_area(0:max_int) = [1,0,1,0,1,1,0,0,1,1,0,0,0,0,0,0,0,0,0,0] !1 - интеграл, который можно вычислить по области, 0 - нельзя
    integer(4), parameter :: int_type_dual(0:max_int) = [-1,-1,-1,-1,5,-2,7,-2,9,-2,11,-2,-1,-1,-1,-1,-1,-1,-1,-1] !номер парного интеграла, -1 - нет парного, -2 текущий интеграл является парным к другому интегралу
    
    type TAreaValue
      real(8), allocatable :: v(:) !(ntr) значение в треугольнике
    endtype
    
    type TAreaValue_c
      real(8) c !постоянное значение
      real(8), allocatable :: v(:) !переменное значение
    endtype
    
    type intf_struct
      integer(4) n
      real(8), allocatable :: x0(:),y0(:),s0(:)
      real(8), pointer :: x(:),y(:),s(:)
      real(8) xder,yder
    endtype
    

  type TCashValue
    real(8), allocatable :: i(:,:) !первый индекс - панель (треугольник) интегрирования
                                   !второй индекс - контрольная точка
  endtype
  
  type TCashIntegral
    logical, allocatable :: inited(:) !(0:max_int) true - массив intvals(k)%i инициализирован
                                      !индекс - тип интеграла
    logical, allocatable :: solved(:) !(0:max_int) true - массив intvals(k)%i вычислен в режиме gs_cash_presolve
                                      !индекс - тип интеграла
    type(TCashValue), allocatable :: intvals(:) !(0:max_int) кэш интегралов
  endtype
  
  type TCash 
    !кэш интегралов, принадлежит границе или области, по которому ведется интегрирование
    logical bnd_inited
    logical area_inited
    type(TCashIntegral), allocatable :: bnd(:) !(nb) кэш интегралов с контрольными точками на границах
    type(TCashIntegral) area !кэш интегралов с контрольными точками в треугольниках
    integer(4) ia            !индекс области, которой принадлежит кэш
  endtype
  
  type TBoungline_GU_gen
    logical is_direct  !прямое или обратное соответствие панелей на участках границы
    integer(4) ia      !индекс второй области
    integer(4) ibnd    !индекс второй границы
    integer(4) ibndl   !индекс второго участка границы
    integer(4) n       !количество неизвестных в правой части
    logical, allocatable :: second_bnd(:) !(n) текущая (false) или другая (true) граница
    integer(4), allocatable :: bndf(:)   !(n) тип функции
                                          !1 - psi
                                          !2 - dpsi/dn
                                          !3 - omega
                                          !4 - domega/dn
    type(TAreaValue_c), allocatable :: carr(:) !(0:n) массивы коэффициентов при неизвестных. 0 - для свободного члена с0
  endtype
  
  type TBoungline_GU_genGlobal
    !граничное условие в зависимости от глобальных констант с перемеными коэффициентами
    integer(4) n       !количество неизвестных в правой части
    type(TAreaValue_c), allocatable :: carr(:) !(0:n) массивы коэффициентов при неизвестных. 0 - для свободного члена с0
    integer(4), allocatable :: gu_indsi(:)  !(n) индексы глобальных неизвестных в gs%constvalinds
  endtype
  
  type TBound_near   !информация о точках для вычисления значений функций вблизи границы
    real(8) x,y      !координаты ближайшей точки на границе
    real(8) dist        !расстоние до точки
    real(8) panelL      !длина ближайшей панели
    logical in_domain   !точка находится в области 
  endtype
  
  type TBound_info !информация о точке на границе
    integer(4) ibnd !номер границы
    real(8) s       !дуговая абсцисса
    real(8) tt      !угол наклона касательной к границе  
    integer(4) in_bound !0 - область, 1 - граница, 2 - окрестность границы (расстояние меньше длины панели)
    type(TBound_near) bn
  endtype

	type TBoundline_GU !граничное условие
	  integer(4) bndu  ! тип граничного условия
		  !0 - ничего не задано 
      !1 - задана функция
      !2 - величина постоянная на всем участке, но заранее неизвестная
      !3 - общее условие - линейная зависимость от других неизвестных (возможно на другом участке, другой границы, другой области)
      !4 - задана функция - const
      !5 - ГУ в зависимсоти от глобальных констант 
    integer(4) bndf  ! функция, для которой задается ГУ
      !1 - psi
      !2 - dpsi/dn
      !3 - omega
      !4 - domega/dn
	  integer(4) constvalind !для bndu=2 - индекс неизвестной в constvals, constvalinds
    type(TBoungline_GU_gen) bndg !для общего ГУ bndu=3
	  real(8), allocatable :: bndval(:,:) !(npanel,nval) значение в граничном условии. nval = 1
    real(8) constval !значение функции для bndu=4
    type(TBoundline), pointer :: bndl !ссылка на участок границы
    type(TBoungline_GU_genGlobal) guglob !для ГУ в зависимсоти от глобальных констант bndu=5
	endtype

	type TBounLineFuncApprox !аппроксимация функции на участке границы
	  !аппроксимируемая функция определяется индексом в массиве TBoundline2%ga
	      !1 - psi
          !2 - dpsi/dn
          !3 - omega
          !4 - domega/dn
	  integer(4) typea ! тип аппроксимации  
	      !1 - кусок ряда фурье (g=2*pi*s, s - локальная относительная дуговая абсцисса: 0 - начало, 1 - конец)
	  integer(4) bnda  !подтиптип аппроксимации  
	      !1 - C0 (заданная)
          !2 - С0
          !3 - С1*sin(g)+С2*cos(g)
          !4 - С0+С1*sin(g)+С2*cos(g)
	  real(8), allocatable :: valkk(:)  !значение коэффициента с учетом умножения на kk
	  real(8), allocatable :: val(:)  !значение коэффициента
	  integer(4), allocatable :: ind(:) !индексы неизвестных для bnda>1
	  real(8), allocatable :: kk(:)  !масштабный множитель для коэффициентов
	  integer(4) n !количество неизвестных коэффициентов
    integer(4) ga_ref !ссылка на индекс родительской аппроксимации (когда dom/dn зависит от om)
	  real(8) cc !амплитуда 
	  real(8) delta !смещение оси (аргумент вектора скорости)    
	endtype

	type TBoundline_Collocate !точки коллокации в граничном условие для аппроксимированной функции
	  !номер уравнения определяется индексом в массиве TBoundline2%gc
	  real(8), pointer :: s_gu(:) !локальная относительная дуговая абсцисса точек, в которых необходимо выполнить уравнение
	  integer(4) n !число точек каллокаций в s_gu
	  real(8), allocatable :: g(:) !s_gu*pi2
	  complex(8),allocatable :: z(:) !коорддинаты точек коллокации
	  complex(8),allocatable :: ztc(:) !вектора от центра круга до точек коллокации
  endtype
  
  type TBoundline_geomdetale
    integer(4) mode !0 - не задано
                    !1 - массив координат
                    !2 - прямолинейный отрезок
                    !3 - дуга окружности
    real(8) s1,s2   !дуговые абсциссы начальной и конечной точек
    real(8) x,y     !первая точка (mode=1), центр окружности (2)
    real(8) x2,y2   !конечная точка (1)
    real(8) r,gam1,gam2 !радиус и углы (2)
    integer(4) i_begin, i_end !индексы первой и последней точки в массивах x,y у TBound (в начале у TBoundline)
    integer(4) ibndl    !индекс TBoundline
    real(8) panelLmax   !максимальная длина панели на участке
  endtype

  type Connect_info
    logical is_internal   !внутренняя граница между сабдоменами
    integer(4) ia         !номер области
    integer(4) ibnd       !номер границы
    integer(4) ibndl      !номер участка
  endtype
    
  type TBoundline !участок линии границы
	    integer(4) i !номер
        !геометрические характеристики
		integer(4) i_begin,i_end       !индекс первой и последней панели участка в массивах границы
    integer(4) npanel    !число панелей
    real(8), allocatable :: x(:), y(:)   !(npanel+1) координаты концов панелей        
        !граничные условия
		type(TBoundline_GU), allocatable :: gu(:) !(nu) граничные условия
    type(TBoundline_GU), allocatable :: gu_add(:) !(nu) альтернативные граничные условия, заменяющие точки коллокации и уравнения в них
		logical, allocatable :: use_gu(:,:)   !(npanel,nu) - true в тех панелях, где нужно записать уравнение
    logical, allocatable :: cp_line(:)   !(nu) - true для того уравнения, для которого нужно создать точки коллокации
    integer(4) gu_mode  !код граничного условия - определяется самим пользователем
    integer(4) igd_begin !индекс начального geomdetale в bnd%geom_detale0
    integer(4) igd_end !индекс конечного geomdetale в bnd%geom_detale0
    logical gu_inited  !граничные условия проинициализрованы. используется для pg_init_gu_multiarea
    type(Connect_info) ci  !информация о соседней границе
  endtype

	type TBoundline2 !участок линии границы
	    integer(4) i !номер
        !геометрические характеристики
		complex(8) zc !центр частицы
		real(8) r !радиус частицы
		real(8) g0 !начальный угол (отсчитывается от оси x против часовой стрелки)
		integer(4) dir !направление обхода окружности: 1 - против часовой стрелки (для внутренности), -1 - по часовой стрелке (для внешности)
		type(TBoundline_Collocate), allocatable :: gc(:) !(nu) точки коллокации граничного условия для аппроксимированной функции
		!граничные условия
		type(TBounLineFuncApprox), allocatable :: ga(:) !(umax) аппроксимации для функций
    integer(4) type_ga !0 - 6 неизвестных, 1 - 4 неизвестных (dom/dn выражается через om)
    real(8) m1 !локальное значение концентрации = 1-m
    real(8) u !локальное значение скорости фильтрации
    real(8) tet !угол наклона локальной скорости фильтрации
  endtype

	type TBound_Collocate_points !точки коллокации в граничном условие для аппроксимированной функции
	  !номер уравнения определяется индексом в массиве TArea%cpp
	  integer(4) ncp !число точек коллокации
	  integer(4), allocatable :: i(:) !(ncp) индекс точки (boundLineType=1 - индекс панели, boundLineType=2 - индекс точки)
	  integer(4), allocatable :: ibndl(:) !(ncp) индекс участка границы (для boundLineType /= 1)
	  integer(4), allocatable :: ibnd(:) !(ncp) индекс границы 
	  integer(4), allocatable :: iu(:) !(ncp) уравнение (1..nu)
	  integer(4), allocatable :: j(:) !(ncp) индекс уравнения в матрице
    logical inited
	endtype

	type TArea_Collocate_points !точки коллокации в граничном условие для аппроксимированной функции
	  !номер уравнения определяется индексом в массиве TArea%cppa
	  integer(4) ncp !число точек коллокации
	  integer(4), allocatable :: iu(:) !(ncp) номер уравнения в треугольнике
	  integer(4), allocatable :: itr(:) !(ncp) индекс треугольника
	  integer(4), allocatable :: j(:) !(ncp) индекс уравнения
    logical inited
	endtype

    type TBound  !линия граница
	    integer(4) i !номер
	    integer(4) nline !число участков границы
		integer(4) boundLineType !тип

		!***boundLineType=1 - граница с набором пнелей,
	    type(TBoundline), allocatable :: line(:)  !участки границы
      type(TBoundline_geomdetale), allocatable :: geom_detale0(:) !информация о геометрии участков границы (исходные)
      type(TBoundline_geomdetale), allocatable :: geom_detale(:) !информация о геометрии участков границы - после объединения, отсортированные
      integer(4) ngeom_detale0 !индекс последнего инициализированного элемента в geom_detale0
      integer(4) ngeom_detale !количество элементов в geom_detale
		!геометрические характеристики
	    integer(4) npanel    !число панелей
        real(8), allocatable :: x(:), y(:)   !(npanel+1) координаты концов панелей
        real(8), allocatable :: xc(:), yc(:) !(npanel)   координаты середин панелей
        real(8), allocatable :: s(:)         !(npanel+1) дуговые абсциссы концов
        real(8), allocatable :: sc(:)        !(npanel)   дуговые абсциссы середин
        real(8), allocatable :: l(:)         !(npanel)   длины панелей
        complex(8), allocatable :: z(:)      !(npanel+1) комплексные координаты концов панелей
        complex(8), allocatable :: zc(:)     !(npanel) комплексные координаты середин панелей
        complex(8), allocatable :: ett(:)    !(npanel) cdexp(ii*tt)
	  
	    !неизвестные
        real(8), allocatable :: psiom(:,:)     !(npanel,umax) - значения искомых функций в контрольных точках
        integer(4), allocatable :: psiind(:,:) !(npanel,umax) - индексы неизвестных
                             !1-й индекс - номер контрольной точки
                             !2-й индекс - 1- psi, 2- dpsi/dn, 3- om, 4- dom/dn
		real(8), allocatable :: bb_psioml(:,:) !(npanel+1,umax)
		real(8), allocatable :: cc_psioml(:,:,:) !(4,npanel+1,umax)

		!***boundLineType=2 - круговая частица с аппроксимацией искомых функций
		type(TBoundline2), allocatable :: line2(:)  !участки границы
    
    !кэш интегралов
    type(TCash) self_cash !собственный кэш
    type(TCash), pointer :: cash !ссылка на кэш другой границы или собственный кэш
    logical skip_getfun !пропустить границу при вычислении функции
  endtype
    
    type areapart
	    !геометрия участок
      integer(4) i !номер
      integer(4) n_begin !первый индекс в массиве zm области
      integer(4) ntr_begin !первый индекс в массиве trm области
      integer(4) n_end !последний индекс в массиве zm области
      integer(4) ntr_end !последний индекс в массиве trm области
      integer(4) npe !число углов в элементе (3,4) - для выделения памяти
      integer(4) n !число узлов
      integer(4) ntr !количество треугольников в сетке
      complex(8), allocatable :: zm(:) !(n) !координаты узлов сетки в области
      integer(4), allocatable :: trm(:,:) !(npe,ntr) !индексы вершин треугольной сетки
      !ссылка на другой участок
      type(areapart), pointer :: apref
      complex(8) a_ref,b_ref,c_ref
    endtype
    
    type areatype
	  !геометрия
      integer(4) npart !количество участков
      type(areapart), allocatable :: part(:) !(npart) участки области
      logical geom_inited
      integer(4) npe !максимальное число углов в элементе (3,4) - для выделения памяти
      integer(4), allocatable :: npe_ar(:) !(ntr) число углов в элементе (3,4)
      integer(4) n !число узлов
      complex(8), allocatable :: zm(:) !(n) !координаты узлов сетки в области
      integer(4), allocatable :: trm(:,:) !(npe,ntr) !индексы вершин треугольной сетки
      integer(4) ntr !количество треугольников в сетке
      complex(8), allocatable :: zmc(:) !(ntr) !координаты центров треугольников
	  !значения в треугольниках
      type(TAreaValue), allocatable :: areaval(:) !(k) !значения в треугольниках для задачи пуассона и т.п. (k=1,2,3,4 в зависимости от задачи)
                                           !1 - значение заданной функции в правой части неоднородного уравнения
    									   !2 - коэффициент (может быть переменным) при искомой функции в уравнении Гельмгольца
    									   !3 - значение функции v в первом уравнении системы вида \Delta\psi=v(x,y)\omega
										   !4 - значение функции vx в уравнении переноса или коэф при df/dy в частном случае осесимметричного уравнения (eq=18,19)
										                                                !или коэф при 1/y*df/dy в частном случае осесимметричного уравнения (eq=20) 
										   !5 - значение функции vy в уравнении переноса
                       
                       !для type_eq=22,24,25 (1 - неоднородное ур, 2 - коэф при psi, 3 - коэф при eta, 4 - dpsi/dx, 5 - dpsi/dy, 6 - deta/dx, 7 - deta/dy)
                       !для type_eq=26 (1 - неоднородное ур, 2 - коэф при psi, 3 - dpsi/dx, 4 - dpsi/dy)
      type(TAreaValue), allocatable :: arval_funp(:) !(nu) !значение функции в центре треугольников при вычислении решения уравнения Пуассона или Гельмгольца
                                              !второй индекс - номер искомой функции в системе (для уравнений 4 порядка или систем)
      !неизвестные (инициализируются только для задач с неизвестными в области)
      integer(4) areaEq                    !тип уравнения (для области)
      integer(4), allocatable :: areaind(:,:) !(ntr,umaxtr) !индекс неизвестной (номер столбца) в области в системе уравнений      
      real(8), allocatable :: psiarea(:,:) !(ntr,umaxtr) !значение искомой функции в области
	                                       !второй индекс=1,2
										   !1 - искомая функция в уравнениях типа гельмгольца
										   !1,2 - df/dx,df/dy в уравнении переноса
      integer(4), allocatable :: eq_ind(:) !(umaxtr) индекс функции, для которой строится текущее уравнение (TArea_Collocate_points%iu)
                       
      !кэш интегралов
      type(TCash) self_cash !собственный кэш
      type(TCash), pointer :: cash !ссылка на кэш другой границы или собственный кэш
      
      !точки коллокации 
		  type(TArea_Collocate_points) cppa  !точки коллокации - точки составления уравнения в области
    endtype
    
    type fict_var !данные для пересчета фиктивных неизвестных
      integer mode !1 - для gu_gen, 2 - для gu_genGlob
      type(TBoungline_GU_gen), pointer :: gu_gen !общее граничное условие для фиктивной неизвестной, где она выражена через реальные неизвестные
      type(TBoungline_GU_genGlobal), pointer :: gu_genGlob !общее граничное условие для фиктивной неизвестной, где она выражена через реальные неизвестные
      integer(4), allocatable :: psiind(:) !(gu_gen%n) индекс неизвестной
      real(8), allocatable :: psiom(:)     !(gu_gen%n) - значения функций, если они заданы
      integer(4) i !индекс для массива переменных коэффициентов (совпадает с номером панели на участке)
    endtype
  
  type col_info
    integer(4) ncol  !количество неизвестных в области вне блока
    integer(4) icol_area !индекс последнего элемента в массиве col, который стоит перед блоком
    integer(4), allocatable :: col(:) !(ncol)
    integer(4) ncol_fict
    integer(4), allocatable :: col_fict(:) !(ncol_fict)
  end type
    
  type matrix_mb
	  integer(4) ix !номер первого столбца неизвестных в глобальной матрице
    integer(4) ix2 !номер последнего столбца неизвестных в глобальной матрице
    integer(4) iu !номер первого уравнения в глобальной матрице
    integer(4) iu2 !номер последнего уравнения в глобальной матрице
    integer(4) nx !число неизвестных (столбцов)
    integer(4) nu !число уравнений (строк) составленных при пробегании по элементам в pg_get_matrix
    integer(4) nu_all !число уравнений, которое должно быть
	  integer(4) nub !количество уравнений на границе 
    integer(4) area_sm_count  !количество неизвестных в области (для matrix_main%sparse%m)
    type(col_info), pointer :: ci
    type(col_info) self_ci
    integer(4) nnt !счетчик добавленных уравнений в блок области
    integer(4) nu_add !число дополнительных уравнений для переопределенной матрицы
  endtype
  
  type TOMP_buffer
    real(8), allocatable :: eq(:) !gs%m%nx_all
    real(8) b
    integer(4) i  !номер текущего уравнения
    integer(4), allocatable :: buffc(:) !(nx) буффер для индексов
  endtype
  
  type TOMP_buffers
    integer(4) max_threads
    type(TOMP_buffer), allocatable :: b(:) !(max_threads)
  endtype
    
	type matrix_main
      integer(4) matrix_store_type  !0 - плотная глобальная матрица
                                    !1 - разреженная глобальная матрица
                                    !2 - блочно-диагональная плотно-разреженная матрица
      real(8), allocatable :: m(:,:)
      type(sparse_matrix) sparse
      real(8), allocatable :: b(:)
      integer(4) nx !число неизвестных (столбцов)
      integer(4) nu !число уравнений (строк)
      integer(4) ix_fict !номер столбца, с которого будут начинаться фиктивные переменные
      integer(4) nx_fict !число фиктивных переменных
      integer(4) nx_all  !общее число неизвестных c фиктивныи
      type(fict_var), allocatable :: fvar(:) !(nx_fict) фиктивные переменные
      real(8) norm     !норма невязки ||Ax-b|| для переопределенной системя
      logical find_norm !найти норма невязки
      real(8) sparse_percent !процент ненулевых (или диагональных элементов) от всех элементов
      real(8) memory_mb !размер в мегабайтах массивов, необходимых для хранения матрицы
      integer(4) allocated_matrix_store_type !-1 - если не инициализировано или =matrix_store_type, если инициализировано
      integer(4) count_closing  !количество элементов, выделяемых под замыкающие уравнения, если 0, то предполагаются все ненулевыми (=nx)
      integer(4) nc_global          !количество ненулевых коэффициентов (при предварительном выделении памяти для sparse)
      type(main_block_matrix) bm   !блочно-диагональная матрица для matrix_store_type=2
      logical have_initial_approxsolv
      integer(4) max_threads
      integer(4) max_ncol
      type(TOMP_buffers) omp_buff
      real(8), allocatable :: res(:)  !решение СЛАУ
    endtype
    
	type TConstArea
    real(8) k_helm     !коэффициент в задаче Гельмгольца с постоянным коэффициентом
		real(8) k_oss      !коэффициент в овесимметричной задаче (=1 для уравнения Лапласа, =-1 для уравнения для функции тока)
    real(8) Re         !число Рейнольдса для уравнений Навье-Стокса
    real(8) k          !проницаемость
    real(8) mub        !вязкость Бринкмана
	endtype

    type TArea !область
	    integer(4) i !номер
      integer(4) mode !признак области, задаваемый пользователем
	    integer(4), allocatable :: type_eq(:)  !(nu) типы уравнений в области
          !1 - гармоническая функция
		      !2 - гармоническая функция (второе уравнение в бигармоническом уравнении)
          !3 - бигармоническая функция
          !4 - уравнение Пуассона
          !5 - однородное уравнение Гельмгольца (через неизвестные в области)
          !6 - неоднородное уравнение Гельмгольца (через неизвестные в области)
          !7 - уравнение Пуассона (второе уравнение в неоднородном бигармоническом уравнении)
          !8 - система вида \Delta\psi=v(x,y)\omega  и уравнение на \omega. неизвестные только на границе и для \omega в области 
          !9 - однородное уравнение Гельмгольца через функции Бесселя \Delta\psi+k_helm^2\psi=0
          !10 - однородное уравнение Гельмгольца через функции Бесселя \Delta\psi-k_helm^2\psi=0
          !11 - система вида \Delta\psi=v(x,y)\omega  и уравнение на \omega. неизвестные только на границе
                !тип уравнения на \omega определяется значением type_eq(2)
          !12 - второе однородное уравнение Гельмгольца через функции Бесселя \Delta\psi+k_helm^2\psi=0 (в системе, также к которой сводится уравнение Бринкмана)
          !13 - второе однородное уравнение Гельмгольца через функции Бесселя \Delta\psi-k_helm^2\psi=0 (в системе, также к которой сводится уравнение Бринкмана)
          !14 - второе однородное уравнение Гельмгольца (через неизвестные в области)
		      !15 - уравнение Бринкмана \Delta^2\psi-k_helm^2\Delta\psi=0 через единую функцию Грина, содержащую функции Бесселя
		      !16 - однородное уравнение переноса \Delta C = vx*dC/dx + vy*dC/dy
		      !17 - неоднородное уравнение переноса \Delta C = vx*dC/dx + vy*dC/dy + f
		      !18 - неоднородное уравнение для осесимметричного случая \Delta\psi = vy*d\psi/dy + f
		      !19 - уравнение Лапласа для осесимметричного случая \Delta\psi = -(k_oss/y)*d\psi/dy  (k_oss=1 - ур-ие Лапласа, k_oss=-1 - ур-ие для функции тока)
		      !20 - уравнение Лапласа для осесимметричного случая \Delta\psi = -(k_oss/y)*d\psi/dy  (k_oss=1 - ур-ие Лапласа, k_oss=-1 - ур-ие для функции тока) - случай интегрирования G_1/y
          !21 - неоднородное бигармоническое уравнение (неизвестные только на границе) или уравнение Навье-Стокса
          !22 - неоднородное бигармоническое уравнение с неизвестными в правой части
                !\Delta^2\psi=f1+f2\psi+f3\eta+f4 d\psi/dx+f5 d\psi/dy+f6 d\eta/dx+f7 d\eta/dy
          !23 - второе уравнение для (22) неоднородного бигармонического уравнение с неизвестными в правой части
          !24 - уравнение Навье-Стокса - неоднородное бигармоническое с неизвестными d\psi/dx,d\psi/dy
          !25 - уравнение Навье-Стокса - неоднородное бигармоническое с неизвестными d\om/dx,d\om/dy
          !26 - уравнение Пуассона с неизвестными в правой части
                !\Delta\psi=f1+f2\psi+f3 d\psi/dx+f4 d\psi/dy
      integer(4) n_eq_var !количество неизвестных в правой части в уравнениях типа 22,24,25,26
      logical, allocatable :: eq_var(:) !(n_eq_var) true - неизвестная присутствует
      integer(4), allocatable :: var_ind(:) !(n_eq_var) второй индекс для неизвестной в areaind
      integer(4) nb  !число границ (связность области)
      type(TBound), allocatable :: bnd(:)
      integer(4) nu ! =1,2 порядок уравнения или системы (число уравнений второго порядка)
      integer(4) umax !=nu*2 - количество неизвестных в каждой точке границы
		  integer(4) umaxtr !количество неизвестных в каждом треугольнике
      type(matrix_mb) m  !матрица
      type(areatype) a   !область с треугольниками
		  logical haveAreaEq
		  type(TConstArea) const   !константы

      !точки коллокации 
		  type(TBound_Collocate_points) cpp  !точки коллокации - точки составления уравнения на границе

		  !real(8), pointer :: func_inf
		  !pointer (pfunc_inf,func_inf)
      real(8) cash_size_mb
      logical need_cash_bb(0:max_int),need_cash_ab(0:max_int),need_cash_aa(0:max_int),need_cash_ba(0:max_int)
      !информация о области, полученной копированием другой области
      real(8) dx,dy !смещение
      real(8) sina,cosa !поворот
      type(TArea), pointer :: a_ref !ссылка на область, из которой была скопирована геометрия (на одну из предыдущих областей, не будущую!!!)
      integer(4) type_rotate !тип перемещения области
                       !-1 - нет ссылки
                       !0 - просто перенос
                       !1 - поворот на 180 град
                       !2 - поворот на 90
                       !3 - поворот на 270
                       !4 - произвольный поворот
      logical oss_resolve !перевычислять кэш обычным способом
      integer(4) iorder !номер области при составлении графа для упорядочивания блоков областей в блочно-диагональной матрице
    endtype
    
	type TConstGlobal
		real(8) ds0numint  !шаг численного интегрирования
    logical use_numerical_int   !использовать численное интегрирование
		logical test_izlom !проверять на излом точки границы при построении сплайнов функций на грание
		real(8) angle_izlom !максимальный угол излома в градусах
		integer(4) di !сколько точек скраю не берем для производных, т.к. там содержатся большие прогрешности
    real(8) epsMinDistBound !минимальное расстояние от точки до границы
    real(8) epsMinDistBound2 !минимальное расстояние от точки до границы (если удовлетворяется условие с этим eps, то дальнейший поиск прекращается)
    real(8) epsMinDistBound3 !минимальное относительное расстояние от точки до границы (для определения in_bound=2)
  endtype
  
  type Tsndns_type
    integer(4) n !число точек коллокаций
    real(8) shift !сдвиг дуговой абсциссы первой точки
    logical use_shift2 !true - точки коллокации для второго уравнения сдвинуты относительно первого
    real(8), allocatable :: s_ndns(:,:) !дуговые абсциссы точек выполнения стандартных ГУ для boundLineType=2
                                        !первый индекс - уравнения
                                        !второй индекс - точка коллокации
  endtype
    
    type TGreenSolve   !панельный метод - одна задача
	    integer(4) i !номер
        integer(4) na     !количество областей
        type(TArea), allocatable :: a(:)  !области
        type(matrix_main) m  !общая матрица
        type(TConstGlobal) const   !константы
        real(8) cash_size_mb
        integer(4) constvaln 
        integer(4) constvalk 
		    integer(4), allocatable :: constvalinds(:)
		    real(8), allocatable :: constvals(:)
        integer(4), allocatable :: constvala(:) !индексы областей, к кторотым должны быть приписаны переменные (для matrix_store_type=2)
                                               !=0 если не инициализировано 
        integer(4), allocatable :: constvalinds_fict(:) !фиктивные (отрицательные) индексы для неизвестных, которые встретились в одной области, а должны быть отнесены к другой 
        integer(4) constvalinds_fict_i !счетчик для фиктивных индексов
    endtype
    
    type TGreenSolves   !панельный метод
	    integer(4) ngs  !число задач
      type(TGreenSolve), allocatable :: ggs(:)  !задачи
      integer(4) nsn !число вариантов sn
      type(Tsndns_type), allocatable :: sn(:) !дуговые абсциссы точек выполнения стандартных ГУ для boundLineType=2
    endtype
	
    type(TGreenSolves), target :: gsMain
	type(TGreenSolve), pointer :: gs  !текущая задача
	type(TArea), pointer :: gsarea    !текущая область текущей задачи
  type(areapart), pointer :: gsareapart    !текущий участок области
	type(TBound), pointer :: gsbnd    !текущая граница    
	type(TConstArea), pointer :: gsareaconst    !константы текущей области
	type(TBoundline), pointer :: gsbndl    !текущий участок границы
	type(TBoundline2), pointer :: gsbndl2    !текущий участок границы

	!для вывода на экран
  logical gs_use_mkl_progress
	real(8) gs_time0, gs_time_interval, gs_time1
    integer gs_use_parallel_matrix_build
    integer gs_use_parallel_get_fun
    integer gs_use_parallel_num_int
    logical gs_stop_on_int_dxdy_inbound
    logical gs_use_cash
    logical gs_cash_presolve              !предварительное вычисление кэша перед составлением матрицы
    logical gs_use_cash_lock_get
    logical gs_write_matrix
    logical gs_rebuild_matrix_for_iterate
    logical first_matrix_build
    logical use_global_sparsem   !использовать сразу глобальную разреженную матрицу
    logical max_equation_coef_eq_1 !приводить коэффициенты матрицы к максимальному по модулю значению =1
    real(8) equation_coef_eq_0_eps !проверка на машинный ноль при добавлении коэффициентов уравнения в матрицу (=d0 - не учитывать)
    real(8) gs_block_matrix_solver_eps  !eps для решателя блочно-диагональной матрицы (calc_slau_block_diagonal)
    integer(4) gs_block_matrix_solver_maxIter  !максимальное число итераций для решателя блочно-диагональной матрицы (calc_slau_block_diagonal)
    integer(4) gs_block_matrix_solver_maxIter_convergence !макимальное количество итераций, за которое не достигается уменьшение погрешности (calc_slau_block_diagonal)
    real(8) gs_block_matrix_solver_lambda  !параметр релаксации x=lam*x[i]+(1-lam)*x[i-1] (calc_slau_block_diagonal)
    integer(4) gs_pardiso_eps !(см. slau.f90 pardiso_eps) 
    integer(4) gs_system_count,gs_system_i !количество и текущая решаемая СЛАУ (для mkl_progreess)
    logical gs_write_log_file
    logical :: gs_used_pg_start=.false.
    logical :: gs_cash_is_solve !внутренняя переменная, используется при вычислении кэша
    logical gs_profiling        !для быстрого нахождения мест профилировки (после оптимизации - удалить везде)
    logical gs_test_areaval_eq_0 !проверять areaval чтобы не заводить лишних неизвестных
    logical gs_use_dual_integral_solve !вычислять два интеграла, являющиеся Re и Im одного комплексного выражения, за один раз при предварительном вычислении кэша
    logical gs_test_point_collenear_for_BesselG  !аналитическое вычисление интегралов с Бесселевыми функциями (для точек лежащих на одной прямой и т.п.)
    logical gs_test_point_near_bound !вычисление функций вблизи границы на расстоянии меньше длины панели
    
    
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
    integer(4) npanel !количество панелей в текущем участке
    integer(4) npanel_k  !количество панелей на сгущающемся участке
    real(8) x1,y1,x2,y2  !координаты концов отрезка
    real(8) ds1          !шаг в начале отрезка
    real(8) ds2          !шаг в конце (или в середине) отрезка
    real(8) kk           !коэффициент изменения длины (>1!!!)
    integer(4) mode      !тип сгущения
    					 !1 - -- --- ---- ---- ----
    					 !3 - -- --- --- --- --- -- - (получается гарантированное четное число панелей)
    integer(4) min_np    !минимальное число панелей
    					 !если получается число панелей,меньше минимального, 
    					 !то делается равномерное разбиение с минимальным числом панелей
    type(TAreaValue), target :: vx,vy
    end
    
    recursive subroutine init_areapart_geom_ref(aref,aa,bb,cc)
    !инициализировать участок ссылкой на другой участок и преобразованием координат
    !z1=A*z+B*conj(z)+C
    import :: areapart
    complex(8) aa,bb,cc
    type(areapart), target :: aref
    end
    
    subroutine init_boundline_geomlineds2_array(ds1,ds2,x1,y1,x2,y2,vx,vy,npanel,min_np)
    !инициализировать геометрию для участка границы в виде прямой линии с переменным ds
    !меняющимся от ds1 до ds2
    import :: TAreaValue
    integer(4) npanel !количество панелей в текущем участке
    real(8) x1,y1,x2,y2  !координаты концов отрезка
    real(8) ds1          !шаг в начале отрезка
    real(8) ds2          !шаг в конце отрезка
    integer(4) min_np    !минимальное число панелей (если 0 - не учитывается)
    type(TAreaValue), target :: vx,vy
    end
    
    function get_bb_cash_integral_2(a,knd,j,knd0,i,nf) result(res)
    !кэш интегрирование по панели контрольная точка на панели
    !ускоренный вариант без проверок
    Import :: TArea
    !j - номер тек панели, по которой идет интегрирование
    !i - номер контрольной точки (связан с номером уравнения)
    !nf - номер интеграла
    !knd0 - номера области и границы с контрольной точкой
    !a   - область где находится панель интегрирования
    !knd - номер границы, где находится панель интегрирования
    integer(4) knd,j,knd0,i,nf
    real(8) res
    !DEC$ IF DEFINED (DEBUG)
    type(TArea), target :: a
    !DEC$ ELSE
    type(TArea) a
    !DEC$ ENDIF
    end
    
    function get_ba_cash_integral_2(a,knd,j,i,nf) result (res)
    !кэш интегрирование по панели контрольная точка в треугольнике
    !ускоренный вариант без проверок
    Import :: TArea
    !j - номер тек панели, по которой идет интегрирование
    !i - номер контрольной точки (треугольника)
    !nf - номер интеграла
    !a   - область где находится панель интегрирования
    !knd - номер границы, где находится панель интегрирования
    integer(4) knd,j,i,nf
    real(8) res
    !DEC$ IF DEFINED (DEBUG)
    type(TArea), target :: a
    !DEC$ ELSE
    type(TArea) a
    !DEC$ ENDIF
    end
    
    function get_ab_cash_integral_2(a,j,knd0,i,nf) result(res)
    !кэш интегрирование по треугольнику контрольная точка на панели
    !ускоренный вариант без проверок
    Import :: TArea
    !j - номер тек треугольника, по которой идет интегрирование
    !i - номер контрольной точки (связан с номером уравнения)
    !nf - номер интеграла
    !knd0 - номер границы с контрольной точкой
    !a - область
    integer(4) j,knd0,i,nf
    real(8) res
    !DEC$ IF DEFINED (DEBUG)
    type(TArea), target :: a
    !DEC$ ELSE
    type(TArea) a
    !DEC$ ENDIF
    end
    
    function get_aa_cash_integral_2(a,j,i,nf) result (res)
    !кэш интегрирование по треугольнику контрольная точка на треугольнике
    !ускоренный вариант без проверок
    Import :: TArea
    !j - номер тек треугольника, по которой идет интегрирование
    !i - номер контрольной точки (треугольника)
    !nf - номер интеграла
    !a - область
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
    integer(4) mode !1 - заданы концы x1,y1,x2,y2,ds
                    !2 - горизонтальная линия x,y1,n
                    !3 - линия массивом координат x,y,n
    type(intf_struct), target :: ifs
    end
    
    subroutine init_boundline_gu_gen_var(gu,i,var)
    Import :: TBoundline_GU
    type(TBoundline_GU), target :: gu
    integer(4) i    !номер массива, =0 для c0
    real(8) var(gu%bndl%npanel)   !массив коэффициентов
    end
    
    subroutine init_boundline_gu_genGlobal_var(gu,i,ind,var)
    Import :: TBoundline_GU
    type(TBoundline_GU), target :: gu
    integer(4) i    !номер массива, =0 для c0
    integer(4) ind  !индекс глобальной переменной в gs%constvalinds
    real(8) var(gu%bndl%npanel)   !массив коэффициентов
    end
    
    subroutine init_boundline_gu_genGlobal_const(gu,i,ind,varconst)
    Import :: TBoundline_GU
    type(TBoundline_GU), target :: gu
    integer(4) i    !номер массива, =0 для c0
    integer(4) ind  !индекс глобальной переменной в gs%constvalinds
    real(8) varconst     !константа для массив коэффициентов
    end
    
    end interface
        
    end

subroutine pg_start
!dec$ attributes dllexport:: pg_start
!инициализация глобальных констант
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
equation_coef_eq_0_eps=1.0d-16  !=1.0d-16 точность double
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
!завершение
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
!инициализировать количество задач
use pgmod
integer(4) n !количество задач
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
!получить следующий индекс в массиве глобальных констант
use pgmod
integer(4) pg_get_next_constind
gs%constvalk=gs%constvalk+1
pg_get_next_constind=gs%constvalk
end

subroutine pg_set_constvala(ia,i1,i2)
!dec$ attributes dllexport:: pg_set_constvala
!инициалищировать индексы областей по диапазону
use pgmod
integer(4) ia !номер области
integer(4) i1,i2 !диапазон в массиве gs%constvala
gs%constvala(i1:i2)=ia
end

subroutine pg_allocate_constvalind(n)
!dec$ attributes dllexport:: pg_allocate_constvalind
!инициализировать низвестные, постоянные на нескольких панелях
use pgmod
integer(4) n !количество неизвестных
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
!выбрать текущую задачу
use pgmod
integer(4) k  !номер объекта задачи
gs=>gsMain%ggs(k)
gs%i=k
end

subroutine pg_bind_domain(k)
!dec$ attributes dllexport:: pg_bind_domain
!выбрать текущую область
use pgmod
integer(4) k  !номер области в текущей задаче
gsarea=>gs%a(k)
gsarea%i=k
end

subroutine pg_bind_areapart(k)
!dec$ attributes dllexport:: pg_bind_areapart
!выбрать текущий участок области
use pgmod
integer(4) k  !номер участка в текущей области
gsareapart=>gsarea%a%part(k)
gsareapart%i=k
end

subroutine pg_bind_bound(k) 
!dec$ attributes dllexport:: pg_bind_bound
!выбрать текущую границу
use pgmod
integer(4) k  !номер границы в текущей области
gsbnd=>gsarea%bnd(k)
gsbnd%i=k
end

subroutine pg_bind_boundline(k) 
!dec$ attributes dllexport:: pg_bind_boundline
!выбрать текущий участок границы
use pgmod
integer(4) k  !номер участка в текущей границе
if (gsbnd%boundLineType==1) then
  gsbndl=>gsbnd%line(k)
  gsbndl%i=k
elseif (gsbnd%boundLineType==2) then
  gsbndl2=>gsbnd%line2(k)
  gsbndl2%i=k
endif
end

subroutine bind_AreaConst(k)
!выбрать текущую область
use pgmod
integer(4) k  !номер области в текущей задаче задачи
gsareaconst=>gs%a(k)%const
end

subroutine initconst_GreenSolve
!инициализировать константы
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
!задание количества дополнительных уравнений для переопределенной СЛАУ
use pgmod
integer(4) n
gsarea%m%nu_add=n
end

subroutine pg_set_matrix_sparse_count_closing(n)
!dec$ attributes dllexport:: pg_set_matrix_sparse_count_closing
!задание количества элементов, выделяемых в разреженной СЛАУ для замыкающих уравнений
use pgmod
integer(4) n
gs%m%count_closing=n+1  !+1 для учета диагонального элемента
end

subroutine pg_set_areaconst(j,val)
!dec$ attributes dllexport:: pg_set_areaconst
!задание констант области
use pgmod
integer(4) j !тип константы
                    !1-k_helm
					!2-k_oss
          !3-Re
real(8) val !значение константы
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
!инициализировать количество областей в задаче и сделать задачу теущей
use pgmod
integer(4) n !количество областей
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
!освободить память во всех задачах
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
!освободить память в задаче
use pgmod
integer(4) i
call deallocate_mm(.true.) !должно быть перед областями, чтобы освободить матрицы областей
if (allocated(gs%a)) then
    do i=1,gs%na
        call deallocate_bounds(i)
    enddo
    deallocate(gs%a)
endif
call deallocate_domains_constval
end

subroutine deallocate_domains_constval
!освободить память gs%constval в задаче
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
!освобождение памяти глобальной матрицы
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
!обнулить данные, связанные с областью
use pgmod
integer(4) ia !номер области
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
!задать тип уравнения для области
use pgmod
integer(4) type_eq  !тип задачи в текущей области
integer(4) type_eq2  !тип второй задачи в текущей области
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
!задать тип уравнения для области (cсистма уравнений)
use pgmod
integer(4) type_eq  !тип задачи в текущей области
integer(4) type_eq2  !тип второй задачи в текущей области
if (type_eq.ne.8.and.type_eq.ne.11) call gs_print_stop('!!! type_eq=8,11')
if (type_eq==8) then
  if (type_eq2.ne.2.and.type_eq2.ne.7.and.type_eq2.ne.12.and.type_eq2.ne.13.and.type_eq2.ne.14) call gs_print_stop('!!! type_eq=2,7,12,13,14')
elseif (type_eq==11) then 
  if (type_eq2.ne.2.and.type_eq2.ne.12.and.type_eq2.ne.13) call gs_print_stop('!!! type_eq=2,12,13')
endif
call set_domain_equation(type_eq,type_eq2)
end 

subroutine set_domain_equation(type_eq,type_eq2)
!задать тип уравнения для области
use pgmod
integer(4) type_eq  !тип задачи в текущей области
integer(4) type_eq2  !тип второй задачи в текущей области
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
!задать массив gsarea%eq_var
use pgmod
integer(4) n !количество элементов в eq_var
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
!выделить память для области
use pgmod
integer(4) inb  !количество границ в текущей области
integer(4) i
gsarea%nb=inb
allocate(gsarea%bnd(inb))
do i=1,inb
  call nullbound(i)
enddo
end

subroutine deallocate_bounds(ia)
!освободить память для области
use pgmod
integer(4) ia !номер области
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
!обнулить данные для границы
use pgmod
integer(4) ibnd !номер гриницы
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
!выделить память для границы
!для boundLineType=1
use pgmod
integer(4) nline !число участков границы
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
!выделить память для границы
!для boundLineType=2
use pgmod
integer(4) nline !число участков границы
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
!освободить память для границы
use pgmod
integer(4) ia   !номер области
integer(4) ibnd !номер гриницы
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
!выделить память для границы (геометрия)
!для BoundLineType=1
use pgmod
type(TBound) b
integer(4) inpanel !количество панелей в текущей границе
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
!выделить память для границы (точки коллокации)
use pgmod
integer(4) incp !число точек коллокации
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
!выделить память для области (точки коллокации)
use pgmod
integer(4) incp !число точек коллокации
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
!обнулить данные для участка границы
use pgmod
integer(4) ibndl !номер участка границы
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
!задать пользовательский код граничного условия
use pgmod
integer(4) ibndl !номер участка границы
integer(4) imode !пользовательский код граничного условия
type(TBoundline), pointer :: b
b=>gsbnd%line(ibndl)
b%gu_mode=imode
end

subroutine nullboundline2(ibndl)
!обнулить данные для участка границы
use pgmod
integer(4) ibndl !номер участка границы
type(TBoundline2), pointer :: b
b=>gsbnd%line2(ibndl)
!!!пока нечего
end

subroutine allocate_boundline_geom(ibndl,inpanel)
!выделить память для участка границы (геометрия)
!для BoundLineType=1
use pgmod
integer(4) ibndl !номер участка границы
integer(4) inpanel !количество панелей в текущей границе
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
!выделить память для участка границы (геометрия)
!для BoundLineType=2
use pgmod
integer(4) ibndl !номер участка границы
type(TBoundline2), pointer :: b
b=>gsbnd%line2(ibndl)
allocate(b%gc(gsarea%nu))
end

subroutine pg_allocate_bound_gu
!dec$ attributes dllexport:: pg_allocate_bound_gu
!выделить память для границы (граничные условия)
!для BoundLineType=1
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
!выделить память для участка границы (альтернативные граничные условия)
!для BoundLineType=1
use pgmod
call pg_allocate_boundline_gu_add2(.false.)
end

subroutine pg_allocate_boundline_gu_add2(cp_line)
!dec$ attributes dllexport:: pg_allocate_boundline_gu_add2
!выделить память для участка границы (альтернативные граничные условия)
!для BoundLineType=1
use pgmod
logical cp_line !true - если на линии создавать точки коллокации при составлении матрицы
integer(4) i
allocate(gsbndl%gu_add(gsarea%nu))
gsbndl%cp_line=cp_line
do i=1,gsarea%nu
  call nullgu(0,i,.true.)
enddo
end

subroutine pg_allocate_bound_gu2
!dec$ attributes dllexport:: pg_allocate_bound_gu2
!выделить память для границы (граничные условия)
!для BoundLineType=2
use pgmod
integer(4) j
type(TBoundline2), pointer :: b
do j=1,gsbnd%nline
  b=>gsbnd%line2(j)
  allocate(b%ga(gsarea%umax)) 
enddo
end

subroutine deallocate_boundline(ia,ibnd,ibndl)
!освободить память для границы
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
!освободить память для границы
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
!освободить память для границы
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
!освободить память для точек коллокации
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
!освободить память для точек коллокации
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
!обнулить данные, связанные с областью
use pgmod
integer(4) igu !номер граничного условия
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
!выделить память для граничного условия
use pgmod
integer(4) igu !номер граничного условия
integer(4) ival !размерность массива значений (пока = 1)
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
!выделить память для граничного условия
use pgmod
integer(4) n !размерность массивов
type(TBoundline_GU) gu !граничное условие
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
!выделить память для граничного условия
use pgmod
integer(4) n !размерность массивов
type(TBoundline_GU) gu !граничное условие
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
!освободить память для граничного условия
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
!освободить память для граничного условия
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
!выделить память для use_gu на текущем участке границы
use pgmod
allocate(gsbndl%use_gu(gsbndl%npanel,gsarea%nu))
gsbndl%use_gu=.true.
end

subroutine pg_allocate_and_set_ga2(iga,n)
!dec$ attributes dllexport:: pg_allocate_and_set_ga2
!выделить память для коэффициентов аппроксимационной функции
!проставить индексы по умолчанию (-1 - нужно индексировать в get_ind_area)
!для BoundLineType=2
use pgmod
integer(4) iga !номер аппрксимационной функции
integer(4) n
type(TBounLineFuncApprox), pointer :: ga
ga=>gsbndl2%ga(iga)
call pg_allocate_and_set_ga2_base(ga,n)
end

subroutine pg_allocate_and_set_ga2_base(ga,n)
!dec$ attributes dllexport:: pg_allocate_and_set_ga2_base
!выделить память для коэффициентов аппроксимационной функции
!проставить индексы по умолчанию (-1 - нужно индексировать в get_ind_area)
!для BoundLineType=2
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
!проставить индексы по умолчанию (-1 - нужно индексировать в get_ind_area)
!для BoundLineType=2
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
!выделить память для точек коллокации
!для BoundLineType=2
use pgmod
integer(4) ibndl
integer(4) igc !номер аппрксимационной функции
integer(4) n   !количество точек коллокации
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
!освободить память для граничного условия
use pgmod
integer(4) ia !номер области
integer(4) ibnd  !номер границы
integer(4) ibndl !номер участка границы
integer(4) igu !номер граничного условия
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
!освободить память для граничного условия
use pgmod
integer(4) ia !номер области
integer(4) ibnd  !номер границы
integer(4) ibndl !номер участка границы
integer(4) igu !номер граничного условия
type(TBoundline_Collocate), pointer :: gc
gc=>gs%a(ia)%bnd(ibnd)%line2(ibndl)%gc(igu)
if (allocated(gc%g)) deallocate(gc%g)
if (allocated(gc%z)) deallocate(gc%z)
if (allocated(gc%ztc)) deallocate(gc%ztc)
end

subroutine deallocate_ga(ia,ibnd,ibndl,iga)
!освободить память для граничного условия
use pgmod
integer(4) ia !номер области
integer(4) ibnd  !номер границы
integer(4) ibndl !номер участка границы
integer(4) iga !номер граничного условия
type(TBounLineFuncApprox), pointer :: ga
ga=>gs%a(ia)%bnd(ibnd)%line2(ibndl)%ga(iga)
call pg_deallocate_ga_base(ga)
end

subroutine pg_deallocate_ga_base(ga)
!dec$ attributes dllexport:: pg_deallocate_ga_base
!освободить память для граничного условия
use pgmod
type(TBounLineFuncApprox) ga
if (allocated(ga%val)) deallocate(ga%val)
if (allocated(ga%ind)) deallocate(ga%ind)
if (allocated(ga%kk)) deallocate(ga%kk)
if (allocated(ga%valkk)) deallocate(ga%valkk)
end

subroutine pg_allocate_area(inpart)
!dec$ attributes dllexport:: pg_allocate_area
!инициализировать количество участков в текущей геометрии
use pgmod
integer(4) inpart   !количество участков
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
!выделить память для геометрии сетки в области
use pgmod
integer(4) in   !количество точек
integer(4) intr !количество треугольников
integer(4) inpe !количество вершин в одном элементе (3 - треугольник, 4 - четырехугольник)
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
!обнулить данные, связанные с областью
use pgmod
integer(4) ia !номер участка области
type(areapart), pointer :: a
a=>gsarea%a%part(ia)
a%n=0
a%ntr=0
a%npe=0
end

subroutine pg_allocate_areapart_geom(in,intr,inpe)
!dec$ attributes dllexport:: pg_allocate_areapart_geom
!выделить память для геометрии сетки в участке области
use pgmod
integer(4) in   !количество точек
integer(4) intr !количество треугольников
integer(4) inpe !количество вершин в одном элементе (3 - треугольник, 4 - четырехугольник)
call pg_allocate_areapart_geom_n(in)
call pg_allocate_areapart_geom_tr(intr,inpe)
end

subroutine pg_allocate_areapart_geom_n(in)
!dec$ attributes dllexport:: pg_allocate_areapart_geom_n
!выделить память для геометрии сетки в участке области
!только узлы
use pgmod
integer(4) in   !количество точек
type(areapart), pointer :: a
a=>gsareapart
a%n=in
allocate(a%zm(in))
a%zm=c0
end

subroutine pg_allocate_areapart_geom_tr(intr,inpe)
!dec$ attributes dllexport:: pg_allocate_areapart_geom_tr
!выделить память для геометрии сетки в участке области
!только элементы
use pgmod
integer(4) intr !количество треугольников
integer(4) inpe !количество вершин в одном элементе (3 - треугольник, 4 - четырехугольник)
type(areapart), pointer :: a
a=>gsareapart
a%ntr=intr
a%npe=inpe
allocate(a%trm(inpe,intr))
a%trm=0
end

subroutine pg_init_areapart_geom_xy(x,y)
!dec$ attributes dllexport:: pg_init_areapart_geom_xy
!инициализировать геометрию сетки в участке области (координаты узлов)
use pgmod
real(8) x(gsareapart%n),y(gsareapart%n) !координаты
gsareapart%zm=dcmplx(x,y)
end

subroutine pg_init_areapart_geom_tr(tr)
!dec$ attributes dllexport:: pg_init_areapart_geom_tr
!инициализировать геометрию сетки в участке области (треугольники или четырехугольники)
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
!выделить память для условий в ячейках (аналог ГУ на границе) в области
use pgmod
integer(4) ku   !тип основного уравнения в области
integer(4) ku2   !тип второго уравнения в области
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
!неоднородные уравнения
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
!уравнения в неизвестными в области
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
!особободить память в области
use pgmod
integer(4) ia   !номер области
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
!особободить память в области
use pgmod
integer(4) ia   !номер области
integer(4) ipart   !номер участка области
type(areapart), pointer :: a
a=>gs%a(ia)%a%part(ipart)
if (allocated(a%zm)) deallocate(a%zm)
if (allocated(a%trm)) deallocate(a%trm)
end

subroutine init_mesh_gu
!задать тип уравнения в ячейках области
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
!задать тип уравнения в ячейках области
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
!значения в треугольниках для задачи пуассона, Гельмгольца и.т.п
!постоянное значение во всех треугольниках
use pgmod
real(8) val !значения в треугольниках 
integer(4) k !тип значения, см. описание areatype%areaval
real(8) g(gsarea%a%ntr) !массив значений
g=val
call pg_init_area_gu(g,k)
end

subroutine pg_init_areapart_gu_const(val,k)
!dec$ attributes dllexport:: pg_init_areapart_gu_const
!значения в треугольниках для задачи пуассона, Гельмгольца и.т.п
!постоянное значение во всех треугольниках
!для участка
use pgmod
real(8) val !значения в треугольниках 
integer(4) k !тип значения, см. описание areatype%areaval
real(8) g(gsareapart%ntr) !массив значений
g=val
call pg_init_areapart_gu(g,k)
end

subroutine pg_init_area_gu_oss
!dec$ attributes dllexport:: pg_init_area_gu_oss
!значения в треугольниках для осесимметричных задач 19,20
use pgmod
integer(4) i
integer(4) k !тип значения, см. описание areatype%areaval
real(8) g(gsarea%a%ntr) !массив значений
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
!значения в треугольниках для задачи пуассона, Гельмгольца и.т.п
use pgmod
integer(4) k !тип значения, см. описание areatype%areaval
real(8) g(gsarea%a%ntr) !массив значений
integer(4) i
type(areapart), pointer :: a
do i=1,gsarea%a%npart
  a=>gsarea%a%part(i)
  call init_areapart_gu(i,g(a%ntr_begin:a%ntr_end),k)
enddo
end

subroutine pg_init_areapart_gu(g,k)
!dec$ attributes dllexport:: pg_init_areapart_gu
!значения в треугольниках для задачи пуассона, Гельмгольца и.т.п
!для участка
use pgmod
integer(4) k !тип значения, см. описание areatype%areaval
real(8) g(gsareapart%ntr) !массив значений
call init_areapart_gu(gsareapart%i,g,k)
end

subroutine pg_init_area_gu_test(g,k,err,ierr,erre)
!dec$ attributes dllexport:: pg_init_area_gu_test
!значения в треугольниках для задачи пуассона, Гельмгольца и.т.п
use pgmod
integer(4) k !тип значения, см. описание areatype%areaval
real(8) g(gsarea%a%ntr) !массив значений
real(8) err !абсолютная пограшность
real(8) erre !относительная погрешность
integer(4) ierr !индекс ячейки, на которой достигнута абсолютная погрешность
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
!значения в треугольниках для задачи пуассона, Гельмгольца и.т.п
!для участка
use pgmod
integer(4) k !тип значения, см. описание areatype%areaval
real(8) g(gsareapart%ntr) !массив значений
real(8) err !абсолютная пограшность
real(8) erre !относительная погрешность
integer(4) ierr !индекс ячейки (с 1 для участка), на которой достигнута абсолютная погрешность
err=d0
erre=d0
ierr=-1
call init_areapart_gu_test(gsareapart%i,g,k,.true.,err,ierr,erre)
erre=erre/gsareapart%ntr
ierr=ierr-gsareapart%ntr_begin+1
end

function test_area_gu(k)
!значения в треугольниках для задачи пуассона, Гельмгольца и.т.п
use pgmod
logical test_area_gu
integer(4) k !тип значения, см. описание areatype%areaval
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
integer(4) ipart !номер участка
integer(4) k !тип значения, см. описание areatype%areaval
real(8) g(gsarea%a%part(ipart)%ntr) !массив значений
real(8) err,erre
integer(4) ierr
call init_areapart_gu_test(ipart,g,k,.false.,err,ierr,erre)
end

subroutine init_areapart_gu_test(ipart,g,k,need_test,err,ierr,erre)
use pgmod
integer(4) ipart !номер участка
integer(4) k !тип значения, см. описание areatype%areaval
real(8) g(gsarea%a%part(ipart)%ntr) !массив значений
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
integer(4) n !число точек коллокации
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
integer(4) n !число точек коллокации
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
!кэш интегрирование по панели контрольная точка на панели
use pgmod
!j - номер тек панели, по которой идет интегрирование
!i - номер контрольной точки (связан с номером уравнения)
!nf - номер интеграла
!knd0 - номера области и границы с контрольной точкой
!ia, knd - номера области и границы, где находится панель интегрирования
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
!кэш интегрирование по панели контрольная точка на панели
use pgmod
!j - номер тек панели, по которой идет интегрирование
!i - номер контрольной точки (связан с номером уравнения)
!nf - номер интеграла
!knd0 - номера области и границы с контрольной точкой
!ia, knd - номера области и границы, где находится панель интегрирования
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
!кэш интегрирование по панели контрольная точка на панели
!ускоренный вариант без проверок
use pgmod
!j - номер тек панели, по которой идет интегрирование
!i - номер контрольной точки (связан с номером уравнения)
!nf - номер интеграла
!knd0 - номера области и границы с контрольной точкой
!a   - область где находится панель интегрирования
!knd - номер границы, где находится панель интегрирования
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
!кэш интегрирование по панели контрольная точка на панели
use pgmod
!j - номер тек панели, по которой идет интегрирование
!i - номер контрольной точки (связан с номером уравнения)
!nf - номер интеграла
!knd0 - номера области и границы с контрольной точкой
!ia, knd - номера области и границы, где находится панель интегрирования
integer(4) ia,knd,j,knd0,i,nf
real(8) val
!$omp critical (lock_cash)
call set_bb_cash_integral_(ia,knd,j,knd0,i,nf,val)
!$omp end critical (lock_cash)
end

subroutine set_bb_cash_integral_(ia,knd,j,knd0,i,nf,val)
!кэш интегрирование по панели контрольная точка на панели
use pgmod
!j - номер тек панели, по которой идет интегрирование
!i - номер контрольной точки (связан с номером уравнения)
!nf - номер интеграла
!knd0 - номера области и границы с контрольной точкой
!ia, knd - номера области и границы, где находится панель интегрирования
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
!кэш интегрирование по панели контрольная точка на панели
!ускоренный вариант без проверок
use pgmod
!j - номер тек панели, по которой идет интегрирование
!i - номер контрольной точки (связан с номером уравнения)
!nf - номер интеграла
!knd0 - номера области и границы с контрольной точкой
!ia, knd - номера области и границы, где находится панель интегрирования
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
!кэш интегрирование по панели контрольная точка в треугольнике
use pgmod
!j - номер тек панели, по которой идет интегрирование
!i - номер контрольной точки (треугольника)
!nf - номер интеграла
!ia, knd - номера области и границы, где находится панель интегрирования
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
!кэш интегрирование по панели контрольная точка в треугольнике
use pgmod
!j - номер тек панели, по которой идет интегрирование
!i - номер контрольной точки (треугольника)
!nf - номер интеграла
!ia, knd - номера области и границы, где находится панель интегрирования
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
!кэш интегрирование по панели контрольная точка в треугольнике
!ускоренный вариант без проверок
use pgmod
!j - номер тек панели, по которой идет интегрирование
!i - номер контрольной точки (треугольника)
!nf - номер интеграла
!a   - область где находится панель интегрирования
!knd - номер границы, где находится панель интегрирования
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
!кэш интегрирование по панели контрольная точка в треугольнике
use pgmod
!j - номер тек панели, по которой идет интегрирование
!i - номер контрольной точки (треугольника)
!nf - номер интеграла
!ia, knd - номера области и границы, где находится панель интегрирования
integer(4) ia,knd,j,i,nf
real(8) val
!$omp critical (lock_cash)
call set_ba_cash_integral_(ia,knd,j,i,nf,val)
!$omp end critical (lock_cash)
end

subroutine set_ba_cash_integral_(ia,knd,j,i,nf,val)
!кэш интегрирование по панели контрольная точка в треугольнике
use pgmod
!j - номер тек панели, по которой идет интегрирование
!i - номер контрольной точки (треугольника)
!nf - номер интеграла
!ia, knd - номера области и границы, где находится панель интегрирования
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
!кэш интегрирование по панели контрольная точка в треугольнике
!ускоренный вариант без проверок
use pgmod
!j - номер тек панели, по которой идет интегрирование
!i - номер контрольной точки (треугольника)
!nf - номер интеграла
!ia, knd - номера области и границы, где находится панель интегрирования
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
!кэш интегрирование по треугольнику контрольная точка на панели
use pgmod
!j - номер тек треугольника, по которой идет интегрирование
!i - номер контрольной точки (связан с номером уравнения)
!nf - номер интеграла
!knd0 - номера области и границы с контрольной точкой
!ia - номера области
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
!кэш интегрирование по треугольнику контрольная точка на панели
use pgmod
!j - номер тек треугольника, по которой идет интегрирование
!i - номер контрольной точки (связан с номером уравнения)
!nf - номер интеграла
!knd0 - номера области и границы с контрольной точкой
!ia - номера области
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
!кэш интегрирование по треугольнику контрольная точка на панели
!ускоренный вариант без проверок
use pgmod
!j - номер тек треугольника, по которой идет интегрирование
!i - номер контрольной точки (связан с номером уравнения)
!nf - номер интеграла
!knd0 - номер границы с контрольной точкой
!a - область
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
!кэш интегрирование по треугольнику контрольная точка на панели
use pgmod
!j - номер тек треугольника, по которой идет интегрирование
!i - номер контрольной точки (связан с номером уравнения)
!nf - номер интеграла
!knd0 - номера области и границы с контрольной точкой
!ia - номера области
integer(4) ia,j,knd0,i,nf
real(8) val
!$omp critical (lock_cash)
call set_ab_cash_integral_(ia,j,knd0,i,nf,val)
!$omp end critical (lock_cash)
end

subroutine set_ab_cash_integral_(ia,j,knd0,i,nf,val)
!кэш интегрирование по треугольнику контрольная точка на панели
use pgmod
!j - номер тек треугольника, по которой идет интегрирование
!i - номер контрольной точки (связан с номером уравнения)
!nf - номер интеграла
!knd0 - номера области и границы с контрольной точкой
!ia - номера области
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
!кэш интегрирование по треугольнику контрольная точка на панели
!ускоренный вариант без проверок
use pgmod
!j - номер тек треугольника, по которой идет интегрирование
!i - номер контрольной точки (связан с номером уравнения)
!nf - номер интеграла
!knd0 - номера области и границы с контрольной точкой
!ia - номера области
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
!кэш интегрирование по треугольнику контрольная точка на треугольнике
use pgmod
!j - номер тек треугольника, по которой идет интегрирование
!i - номер контрольной точки (треугольника)
!nf - номер интеграла
!ia - номера области
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
!кэш интегрирование по треугольнику контрольная точка на треугольнике
use pgmod
!j - номер тек треугольника, по которой идет интегрирование
!i - номер контрольной точки (треугольника)
!nf - номер интеграла
!ia - номера области
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
!кэш интегрирование по треугольнику контрольная точка на треугольнике
!ускоренный вариант без проверок
use pgmod
!j - номер тек треугольника, по которой идет интегрирование
!i - номер контрольной точки (треугольника)
!nf - номер интеграла
!a - область
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
!кэш интегрирование по треугольнику контрольная точка на треугольнике
use pgmod
!j - номер тек треугольника, по которой идет интегрирование
!i - номер контрольной точки (треугольника)
!nf - номер интеграла
!ia - номера области
integer(4) ia,j,i,nf
real(8) val
!$omp critical (lock_cash)
call set_aa_cash_integral_(ia,j,i,nf,val)
!$omp end critical (lock_cash)
end

subroutine set_aa_cash_integral_(ia,j,i,nf,val)
!кэш интегрирование по треугольнику контрольная точка на треугольнике
use pgmod
!j - номер тек треугольника, по которой идет интегрирование
!i - номер контрольной точки (треугольника)
!nf - номер интеграла
!ia - номера области
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
!кэш интегрирование по треугольнику контрольная точка на треугольнике
!ускоренный вариант без проверок
use pgmod
!j - номер тек треугольника, по которой идет интегрирование
!i - номер контрольной точки (треугольника)
!nf - номер интеграла
!ia - номера области
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
!установить ссылку на cash другой области, из которой текущая область получена копированием границы и сетки в области со смещением и поворотом
use pgmod
integer(4) ip !индекс задачи с областью, на которую будет ссылка
integer(4) ia !индекс области, на которую будет ссылка
integer(4) i
type(TArea), pointer :: a2
type(TBound), pointer :: b,b2
if (gsarea%type_rotate<0) then
  call gs_print("Error pg_set_cash_ref. Set cash ref for an uncopied domain (ia="//itoa(ia)//")!")
  call gs_print_stop("Use pg_copy_subdomain_dxdy_rotate");
else if (gsarea%type_rotate==0.and.(.not.gsarea%oss_resolve)) then
  !перенос (кроме осесимметричных задач с переносом по y)
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
!очистить свойство skip_getfun у всех границ текущей области
use pgmod
integer(4) i
do i=1,gsarea%nb
  gsarea%bnd(i)%skip_getfun=.false.
enddo
end

function get_gug_c(av,k) result(cc0)
!получить i-й коэффициент (0 - с0) для k-й панели
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
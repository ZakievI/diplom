subroutine solve !построение графика числа T = N_impact/N от числа стокса
    use mod
    implicit none
    integer(4), parameter::n = 4
    integer(4):: i, num, f_t(3)
    integer(4) function_impact_testing
    real(8), allocatable :: y_out(:,:) 
    real(8):: y(n), dlt, stocs, y_t(3)
    real(8):: pg_get_fun_xy
    external fcn, fcn_s, function_impact_testing
    !open (1, file='traektorie.dat')
    !write(1,*) 'title = "traektorie"'
    !write(1,*) 'variables = "x", "y"'
    open (1, file='grafik.dat')
    write(1,*) 'title = "traektorie"'
    write(1,*) 'variables = "x", "y"'
    write(1,"('ZONE T=""area',i0,'"", I=', i0, ', F=POINT')") 1, 10
    print "(4x, 'istep' , 5x, 'time', 9x, 'x', 11x, 'y')"
    num             = 5000
    dlt             = d0
    stocs           = -1d0
    allocate(y_out(num,2))
    !$omp parallel do if (gs_use_parallel_build_grafic == 1) private(i, n1, ido, s, y, Curve_tempr, param)
    do while (stocs <= 2)
        st              = 10**(stocs)
        !write(*,"('st=')")
        print * ,'st=', stocs
        !st              = 1d0
        y_t(1)          = eps
        y_t(3)          = d1 + dlt + eps
        f_t(2)          = 1
        do i=1,3,2
            y(1)        = -L1/2
            y(2)        = y_t(i)
            y(3)        = pg_get_fun_xy(y(1),y(2),2,d0,d1,0)
            y(4)        = -pg_get_fun_xy(y(1),y(2),2,d1,d0,0)
            call impact_test(n, y , y_out, num, dlt)
            f_t(i) = function_impact_testing(n, y, dlt)
        end do
        do while(f_t(2) /= 0)
            y_t(2)      = (y_t(1)+y_t(3))/2
            y(1)        = -L1/2
            y(2)        = y_t(2)
            y(3)        = pg_get_fun_xy(y(1),y(2),2,d0,d1,0)
            y(4)        = -pg_get_fun_xy(y(1),y(2),2,d1,d0,0)
            call impact_test(n, y , y_out, num, dlt)
            f_t(2) = function_impact_testing(n, y, dlt)
            !print * , f_t
            !print * , y_t
            if ((f_t(1)*f_t(3)) < 0) then
                if (f_t(2) == 0) then
                    !print * , f_t
                    !print * , y_t
                    exit
                else if (f_t(2) == -1) then
                    y_t(3) = y_t(2)
                    f_t(3) = f_t(2)
                else
                    y_t(1) = y_t(2)
                    f_t(1) = f_t(2)
                end if
            else
                y_t(1) = y_t(1)**2
            end if 
            if (y_t(3)-y_t(1) < eps) then
                exit
            end if
        end do
        print *, y_t(2)
        write(1,"(E15.5, ' ', E15.5)") st, y_t(2)
        print *, stocs, ' ', y_t(2)
        stocs = stocs + d1/3
    end do
    !$OMP END PARALLEL DO
    deallocate(y_out)
    end subroutine solve
subroutine build_bound() !постройка границы в области после цилиндра
    use mod
    implicit none
    integer(4), parameter::n = 4
    integer(4):: i, num, f_t(3)
    integer(4) function_impact_testing, size_arr
    real(8), allocatable :: bound_first_particle(:,:)
    real(8), allocatable :: new_arr(:,:)
    real(8):: y(n), dlt, stocs, y_t(3)
    real(8):: pg_get_fun_xy
    external fcn, fcn_s, function_impact_testing
    num             = 5000
    dlt             = d0
    stocs           = -1d0
    allocate(bound_first_particle(num,2))
    st              = 10**(stocs)
    y_t(1)          = eps
    y_t(3)          = d1 + dlt + eps
    f_t(2)          = 1
    do i=1,3,2
        y(1)        = -L1/2
        y(2)        = y_t(i)
        y(3)        = pg_get_fun_xy(y(1),y(2),2,d0,d1,0)
        y(4)        = -pg_get_fun_xy(y(1),y(2),2,d1,d0,0)
        call impact_test(n, y , bound_first_particle, num, dlt)
        f_t(i) = function_impact_testing(n, y, dlt)
    end do
    do while(f_t(2) /= 0)
        y_t(2)      = (y_t(1)+y_t(3))/2
        y(1)        = -L1/2
        y(2)        = y_t(2)
        y(3)        = pg_get_fun_xy(y(1),y(2),2,d0,d1,0)
        y(4)        = -pg_get_fun_xy(y(1),y(2),2,d1,d0,0)
        call impact_test(n, y , bound_first_particle, num, dlt)
        f_t(2) = function_impact_testing(n, y, dlt)
        if ((f_t(1)*f_t(3)) < 0) then
            if (f_t(2) == 0) then
                exit
            else if (f_t(2) == -1) then
                y_t(3) = y_t(2)
                f_t(3) = f_t(2)
            else
                y_t(1) = y_t(2)
                f_t(1) = f_t(2)
            end if
        else
            y_t(1) = y_t(1)**2
        end if 
        if (y_t(3)-y_t(1) < eps) then
            exit
        end if
    end do
    y(1)        = -L1/2
    y(2)        = y_t(3)
    y(3)        = pg_get_fun_xy(y(1),y(2),2,d0,d1,0)
    y(4)        = -pg_get_fun_xy(y(1),y(2),2,d1,d0,0)
    call coordinate_first_particle(n, y, bound_first_particle, num, dlt, size_arr)
    allocate(new_arr(size_arr,2))
    do i= 1, size_arr
        new_arr(i,1) = bound_first_particle(i,1)
        new_arr(i,2) = bound_first_particle(i,2)
    end do
    deallocate(bound_first_particle)
    allocate(bound_first_particle(size_arr,2))
    bound_first_particle = new_arr
    deallocate(new_arr)
    end subroutine build_bound
subroutine fcn(n, t, y, yprime) !интегрирование по времени
    use mod
    integer(4) :: n
    real(8) t, y(n), yprime(n), u_x, u_y, pg_get_fun_xy
    real(8), parameter :: g = 9.81d0, v0 = 10.0d0
    u_x         = pg_get_fun_xy(y(1),y(2),2,d0,d1,0)
    u_y         = -pg_get_fun_xy(y(1),y(2),2,d1,d0,0)
    yprime(1)   = y(3) 
    yprime(2)   = y(4)
    yprime(3)   = (u_x-y(3))/st
    yprime(4)   = (u_y-y(4))/st
    end subroutine fcn
subroutine fcn_s(n, t, y, yprime) !интегрирование по длине дуги s
    use mod
    integer(4) :: n
    real(8) t, y(n), yprime(n), u_x, u_y, pg_get_fun_xy, V
    real(8), parameter :: g = 9.81d0, v0 = 10.0d0
    V           = dsqrt(y(3)**2+y(4)**2)
    u_x         = pg_get_fun_xy(y(1),y(2),2,d0,d1,0)
    u_y         = -pg_get_fun_xy(y(1),y(2),2,d1,d0,0)
    yprime(1)   = y(3)/V
    yprime(2)   = y(4)/V
    yprime(3)   = (u_x-y(3))/(st*V)
    yprime(4)   = (u_y-y(4))/(st*V)
    end subroutine fcn_s
function function_impact_testing(n, y, dlt)  !индикаторная функция
    use mod
    integer(4)::n
    real(8) y(n), dlt, tol
    integer(4) function_impact_testing
    function_impact_testing = 2
    tol = 0.001d0
    if (y(2)<d1+dlt-tol) then
        function_impact_testing = 1
    else if ((y(2)>= d1+dlt-tol).and.(y(2)<= d1+dlt+tol)) then
        function_impact_testing = 0
    else if (y(2)> d1+dlt+tol) then
        function_impact_testing = -1
    end if
    end function function_impact_testing
subroutine impact_test(n, y , y_out, num, dlt) !проверка на ударение чатcицы
    use mod
    integer(4) :: n, num, k1, ido
    integer(4), parameter::mxparm = 100
    real(8) :: y_out(num,2)
    real(8) :: y(n), s, d_s, tol, param(mxparm), dlt
    external fcn, fcn_s
    s           = d0                !��������� ������� ��� �������������� �� ����
    tol         = 0.001d0           !���������� ������
    param       = d0                !��������� �� ���������
    param(4)    = num               !������������ ���-�� ��������
    param(10)   = 1.0d0             !������������ ���������� �������
    d_s         = 0.05d0            !��� �� ����
    k1          = 2
    ido = 1
    do while ((sqrt(y(1)**2+y(2)**2)>d1+dlt).and.(-L1/2<=y(1)).and.(y(1)<d0).and.(y(2)<H1).and.(d0<y(2)).and.(k1<=num))
        call divprk(ido, n, fcn_s, s, s+d_s, tol, param, y)
        !print '(i6, 6f12.3)', k1, t, y
        !y_out(k1,1)=y(1)
        !y_out(k1,2)=y(2)
        k1 = k1 + 1
    end do
    call divprk(3, n, fcn_s, s, s+d_s, tol, param, y)
    end subroutine impact_test
subroutine coordinate_first_particle(n, y, arr_bound, num, dlt, size_arr) !координаты экстремальной частицы, прошедшей мимо цилиндрва  
    use mod 
    integer(4) :: n, num, k1, ido, size_arr
    integer(4), parameter::mxparm = 100
    real(8) :: arr_bound(num,2)
    real(8) :: y(n), s, d_s, tol, param(mxparm), dlt
    external fcn, fcn_s
    s           = d0                !��������� ������� ��� �������������� �� ����
    tol         = 0.001d0           !���������� ������
    param       = d0                !��������� �� ���������
    param(4)    = num               !������������ ���-�� ��������
    param(10)   = 1.0d0             !������������ ���������� �������
    d_s         = 0.05d0            !��� �� ����
    k1          = 2
    ido = 1
    arr_bound(1, 1) = d0
    arr_bound(1, 2) = d1
    do while ((y(1)<d0).and.(k1<=num))
        call divprk(ido, n, fcn_s, s, s+d_s, tol, param, y)
        !print '(i6, 6f12.3)', k1, t, y
        !y_out(k1,1)=y(1)
        !y_out(k1,2)=y(2)
        !k1 = k1 + 1
    end do
    do while ((y(1)<L1/2).and.(k1<=num))
        call divprk(ido, n, fcn_s, s, s+d_s, tol, param, y)
        print *, k1, y(1), y(2)
        arr_bound(k1,1)=y(1)
        arr_bound(k1,2)=y(2)
        k1 = k1 + 1
    end do
    size_arr = k1-1
    call divprk(3, n, fcn_s, s, s+d_s, tol, param, y)
    end subroutine coordinate_first_particle
function area_quadrilateral(z1, z2, z3, z4) !нахождение объема по 4-м точкам
    use mod
    real(8) area_quadrilateral
    complex(8) z1, z2, z3, z4
    ! area_quadrilateral = (dreal(p3-p1)*dimag(p4-2*p1+p2)+dimag(p3-p1)*dreal(p4-2*p1+p2))*d5
    area_quadrilateral = d5 * ABS(DIMAG(z1*CONJG(z2) + z2*CONJG(z3) + z3*CONJG(z4) + z4*CONJG(z1)))
    end function area_quadrilateral
function search_for_extreme_particles() ! поиск критической частицы, после которой след. частицы пролетаю мимо цилиндра
    use mod
    implicit none
    integer(4), parameter::n = 4
    integer(4) :: i, num, f_t(3)
    integer(4) function_impact_testing
    real(8), allocatable :: bound_first_particle(:,:)
    real(8) :: y(n), dlt, y_t(3)
    real(8) :: pg_get_fun_xy
    real(8) :: search_for_extreme_particles
    external fcn_s, function_impact_testing
    num             = 5000
    dlt             = d0
    allocate(bound_first_particle(num,2))
    y_t(1)          = eps
    y_t(3)          = d1 + dlt + eps
    f_t(2)          = 1
    do i=1,3,2
        y(1)        = -L1/2
        y(2)        = y_t(i)
        y(3)        = pg_get_fun_xy(y(1),y(2),2,d0,d1,0)
        y(4)        = -pg_get_fun_xy(y(1),y(2),2,d1,d0,0)
        call impact_test(n, y , bound_first_particle, num, dlt)
        f_t(i) = function_impact_testing(n, y, dlt)
    end do
    do while(f_t(2) /= 0)
        y_t(2)      = (y_t(1)+y_t(3))/2
        y(1)        = -L1/2
        y(2)        = y_t(2)
        y(3)        = pg_get_fun_xy(y(1),y(2),2,d0,d1,0)
        y(4)        = -pg_get_fun_xy(y(1),y(2),2,d1,d0,0)
        call impact_test(n, y , bound_first_particle, num, dlt)
        f_t(2) = function_impact_testing(n, y, dlt)
        if ((f_t(1)*f_t(3)) < 0) then
            if (f_t(2) == 0) then
                exit
            else if (f_t(2) == -1) then
                y_t(3) = y_t(2)
                f_t(3) = f_t(2)
            else
                y_t(1) = y_t(2)
                f_t(1) = f_t(2)
            end if
        else
            y_t(1) = y_t(1)**2
        end if 
        if (y_t(3)-y_t(1) < eps) then
            exit
        end if
    end do
    search_for_extreme_particles = y_t(2)
    deallocate(bound_first_particle)
    end function
subroutine draw_Curves() !вывод кривых
    use mod
    integer(4) :: i, l
    open (1, file='traektorie.dat')
    write(1,*) 'title = "traektorie"'
    write(1,*) 'variables = "x", "y", "t", "c"'
    do i = 1, size(Curves)
        write(1,"('ZONE T=""area',i0,'"", I=', i0, ', F=POINT')") i, Curves(i)%n
        do l = 1, Curves(i)%n
            write(1,"(E15.5, ' ', E15.5, ' ', E15.5, ' ', E15.5)") Curves(i)%x(l), Curves(i)%y(l), Curves(i)%t(l), Curves(i)%c(l)
        end do
    end do
    end subroutine
subroutine build_curve() !поиск кривых
    use mod
    integer(4) :: i, n1, ido
    integer(4), parameter :: n = 5
    real(8) area_quadrilateral, search_for_extreme_particles, S1
    integer(4), parameter :: mxparm = 50
    real(8) :: param(mxparm), d_s, s, y(n), dlt, tol, pg_get_fun_xy
    real(8), allocatable :: Curve_tempr(:,:)
    complex(8) :: p1, p2, p3, p4
    external fcn_s_t, area_quadrilateral, fcn_s_t_top, fcn_s_t_bottom, search_for_extreme_particles
    !p1              = (1d0, 1d0)
    !p2              = (3d0, 1d0)
    !p3              = (3d0, 3d0)
    !p4              = (1d0, 3d0)
    !S1              = area_quadrilateral(p1, p2, p3, p4)
    tol             = 0.0001d0
    d_s             = 0.01d0
    s               = d0
    ido             = 1
    dlt             = d0
    allocate(Curves(num_particle))
    ! cord_extreme_particles = search_for_extreme_particles()
    
    !$omp parallel do if (gs_use_parallel_build_cerves == 1) private(i, n1, ido, s, y, Curve_tempr, param)
    do i = 1, num_particle
        allocate(Curve_tempr(N_arr,4))
        write(*,"('I=',i0)") i
        n1                 = 1
        ido                = 1
        s                  = d0
        Curve_tempr(n1, 1) = -L1/2 
        if ((H1*((i)/(num_particle-d1))**3>cord_extreme_particles).and.(cord_extreme_particles>H1*((i-d1)/(num_particle-d1))**3)) then
            cord_extreme_particles = search_for_extreme_particles()
            Curve_tempr(n1, 2) = cord_extreme_particles
            index_extreme_particles = i
        else
            !Curve_tempr(n1, 2) = H1*(i - d1)/(num_particle - d1)
            Curve_tempr(n1, 2) = H1*((i-d1)/(num_particle-d1))**3
        end if
        Curve_tempr(n1, 3) = d0
        Curve_tempr(n1, 4) = s
        y(1)            = -L1/2
        y(2)            = Curve_tempr(n1, 2)
        if (i == 1) then
            y(3)            = pg_get_fun_xy(y(1),y(2) + eps,2,d0,d1,0)
            y(4)            = -pg_get_fun_xy(y(1),y(2) + eps,2,d1,d0,0)
        else if(i == num_particle)then 
            y(3)            = pg_get_fun_xy(y(1),y(2) - eps,2,d0,d1,0)
            y(4)            = -pg_get_fun_xy(y(1),y(2) -  eps,2,d1,d0,0)
        else 
            y(3)            = pg_get_fun_xy(y(1),y(2),2,d0,d1,0)
            y(4)            = -pg_get_fun_xy(y(1),y(2),2,d1,d0,0)
        end if
        y(5)            = d0
        param           = d0
        param(4)        = N_arr
        param(10)       = 1.0d0
        do while((sqrt(y(1)**2+y(2)**2)>d1+dlt).and.(-L1/2<=y(1)).and.(y(1)<=L1/2).and.(y(2)<=H1).and.(d0<=y(2)).and.(n1<N_arr))
            n1 = n1 + 1
            !write(*,"('N=',i0)") n1
            if (i == 1) then
                call divprk(ido, n, fcn_s_t_bottom, s, s+d_s, tol, param, y)
                y(2) = d0
            else if(i == num_particle)then 
                call divprk(ido, n, fcn_s_t_top, s, s+d_s, tol, param, y)
                y(2) = H1
            else
                call divprk(ido, n, fcn_s_t, s, s+d_s, tol, param, y)
            end if 
            Curve_tempr(n1, 1) = y(1)
            Curve_tempr(n1, 2) = y(2)
            Curve_tempr(n1, 3) = y(5)
            Curve_tempr(n1, 4) = s
        end do
        allocate(Curves(i)%x(n1))
        allocate(Curves(i)%y(n1))
        allocate(Curves(i)%c(n1))
        allocate(Curves(i)%t(n1))
        allocate(Curves(i)%s(n1))
        Curves(i)%x(1:n1) = Curve_tempr(1:n1, 1)
        Curves(i)%y(1:n1) = Curve_tempr(1:n1, 2)
        Curves(i)%t(1:n1) = Curve_tempr(1:n1, 3)
        Curves(i)%s(1:n1) = Curve_tempr(1:n1, 4)
        Curves(i)%n = n1
        call divprk(3, n, fcn_s_t, s, s+d_s, tol, param, y)
        deallocate(Curve_tempr)
    end do 
    !$OMP END PARALLEL DO
    call find_concentration()
    call draw_Curves()
    !call build_mesh_1()
    end subroutine build_curve
subroutine build_mesh_1
    use mod
    integer(4) :: i, l
    real(8), allocatable :: tempr_y_c(:), tempr_x_c(:)
    integer(4) :: num_of_partitions_by_x
    real(8) :: x_c
    N_part_2 = index_extreme_particles
    num_of_partitions_by_x = (N_part_1 + N_part_2 + N_part_3 + N_part_4)
    ! allocate(mesh%x_y_(num_particle * num_of_partitions_by_x, 2),mesh%t(num_particle * num_of_partitions_by_x, 1),mesh%c(num_particle * num_of_partitions_by_x, 1)&
    !         mesh%v(num_particle * num_of_partitions_by_x, 2),mesh%s(num_particle * num_of_partitions_by_x, 2))
    allocate(mesh)
    mesh%n_i = num_of_partitions_by_x
    mesh%n_j = num_particle
    allocate(mesh%x_y_(num_of_partitions_by_x * num_particle, 2))
    ! область 1 - x меняеться от -L1 до -1, y - от 0 до H1
    ! область 2 - x меняеться от -1 до 0, y - от Y_Last до H1
    ! область 3 - x меняеться от 0 до 1, y - от Y_Last до H1
    ! область 4 - x меняеться от 1 до L1, y - от 0 до H1
    
    ! область 1
    do i = 1, 4
        if (i == 1) then
            do l = 1, N_part_1
                x_c = -L1 + (L1 - d1)*(l-1)/(N_part_1-1)
                mesh%x_y_(l, 1) = x_c
            end do
            allocate(tempr_y_c(N_part_1),tempr_x_c(N_part_1))
            tempr_x_c = mesh%x_y_(1:N_part_1, 1)
            do l = 1, num_particle
                call dcsiez(Curves(l)%n, Curves(l)%x, Curves(l)%y, N_part_1, tempr_x_c, tempr_y_c )
                mesh%x_y_(1+num_of_partitions_by_x * (l-1):N_part_1 + num_of_partitions_by_x * (l-1),2) = tempr_y_c
                mesh%x_y_(1+num_of_partitions_by_x * (l-1):N_part_1 + num_of_partitions_by_x * (l-1),1) = mesh%x_y_(1:N_part_1, 1)
            end do
            deallocate(tempr_y_c,tempr_x_c)
        elseif (i == 2)then
            allocate(tempr_x_c(N_part_2 - 1))
            do l = 2, N_part_2
                tempr_x_c(l-1) = Curves(l)%x(Curves(l)%n)
            end do
        elseif (i == 3) then 
        else 
        end if
        ! end do
    end do 
    end subroutine build_mesh_1
subroutine find_concentration() !поиск концентрации
    use mod
    integer(4) :: i, n1, j
    real(8) V_0, area_quadrilateral 
    !real(8), allocatable :: value_(:), xvec(:), xdata(:), fdata(:)
    type(Curve) :: Curve_tempr_top, Curve_tempr_bottom
    complex(8) p1, p2, p3, p4
    external area_quadrilateral
    V_0 = d1
    do i = 1, num_particle
        if ((i == 1).or.(i == index_extreme_particles + 1)) then
            Curve_tempr_top%n = Curves(i)%n
            allocate(Curve_tempr_top%x(Curves(i)%n), Curve_tempr_top%y(Curves(i)%n))
            call dcsiez(Curves(i+1)%n, Curves(i+1)%t, Curves(i+1)%x, Curve_tempr_top%n, Curves(i)%t, Curve_tempr_top%x)!построение сплайна 
            call dcsiez(Curves(i+1)%n, Curves(i+1)%t, Curves(i+1)%y, Curve_tempr_top%n, Curves(i)%t, Curve_tempr_top%y)
            do j = 1, Curves(i)%n
                if (j == 1) then
                    p1 = cmplx(Curves(i)%x(j),Curves(i)%y(j))
                    p2 = cmplx(Curves(i)%x(j+1),Curves(i)%y(j+1))
                    p3 = cmplx(Curve_tempr_top%x(j+1),Curve_tempr_top%y(j+1))
                    p4 = cmplx(Curve_tempr_top%x(j),Curve_tempr_top%y(j))
                    V_0 = 4*area_quadrilateral(p1,p2,p3,p4)
                    Curves(i)%c(j)=d1
                elseif ( j == Curves(i)%n ) then
                    p1 = cmplx(Curves(i)%x(j-1),Curves(i)%y(j-1))
                    p2 = cmplx(Curves(i)%x(j),Curves(i)%y(j))
                    p3 = cmplx(Curve_tempr_top%x(j),Curve_tempr_top%y(j))
                    p4 = cmplx(Curve_tempr_top%x(j-1),Curve_tempr_top%y(j-1))
                    Curves(i)%c(j) = V_0/(4*area_quadrilateral(p1,p2,p3,p4))
                else
                    p1 = cmplx(Curves(i)%x(j-1),Curves(i)%y(j-1))
                    p2 = cmplx(Curves(i)%x(j+1),Curves(i)%y(j+1))
                    p3 = cmplx(Curve_tempr_top%x(j+1),Curve_tempr_top%y(j+1))
                    p4 = cmplx(Curve_tempr_top%x(j-1),Curve_tempr_top%y(j-1))
                    Curves(i)%c(j) = V_0/(2*area_quadrilateral(p1,p2,p3,p4))
                end if
            end do
            call finalize(Curve_tempr_top)
        ! else if (i == index_extreme_particles): 
        else if ((i == num_particle) .or. (i == index_extreme_particles)) then 
            Curve_tempr_bottom%n = Curves(i)%n
            allocate(Curve_tempr_bottom%x(Curves(i)%n), Curve_tempr_bottom%y(Curves(i)%n))
            call dcsiez(Curves(i-1)%n, Curves(i-1)%t, Curves(i-1)%x, Curve_tempr_bottom%n, Curves(i)%t, Curve_tempr_bottom%x)
            call dcsiez(Curves(i-1)%n, Curves(i-1)%t, Curves(i-1)%y, Curve_tempr_bottom%n, Curves(i)%t, Curve_tempr_bottom%y)
            do j = 1, Curves(i)%n
                if (j == 1) then
                    p1 = cmplx(Curve_tempr_bottom%x(j),Curve_tempr_bottom%y(j))
                    p2 = cmplx(Curve_tempr_bottom%x(j+1),Curve_tempr_bottom%y(j+1))
                    p3 = cmplx(Curves(i)%x(j+1),Curves(i)%y(j+1))
                    p4 = cmplx(Curves(i)%x(j),Curves(i)%y(j))
                    V_0 = 4*area_quadrilateral(p1,p2,p3,p4)
                    Curves(i)%c(j)=d1
                elseif ( j == Curves(i)%n ) then
                    p1 = cmplx(Curve_tempr_bottom%x(j-1),Curve_tempr_bottom%y(j-1))
                    p2 = cmplx(Curve_tempr_bottom%x(j),Curve_tempr_bottom%y(j))
                    p3 = cmplx(Curves(i)%x(j),Curves(i)%y(j))
                    p4 = cmplx(Curves(i)%x(j-1),Curves(i)%y(j-1))
                    Curves(i)%c(j) = V_0/(4*area_quadrilateral(p1,p2,p3,p4))
                else
                    p1 = cmplx(Curve_tempr_bottom%x(j-1),Curve_tempr_bottom%y(j-1))
                    p2 = cmplx(Curve_tempr_bottom%x(j+1),Curve_tempr_bottom%y(j+1))
                    p3 = cmplx(Curves(i)%x(j+1),Curves(i)%y(j+1))
                    p4 = cmplx(Curves(i)%x(j-1),Curves(i)%y(j-1))
                    Curves(i)%c(j) = V_0/(2*area_quadrilateral(p1,p2,p3,p4))
                end if
            end do
            call finalize(Curve_tempr_bottom)
        else 
            Curve_tempr_top%n = Curves(i)%n
            Curve_tempr_bottom%n = Curves(i)%n
            allocate(Curve_tempr_top%x(Curves(i)%n), Curve_tempr_top%y(Curves(i)%n))
            allocate(Curve_tempr_bottom%x(Curves(i)%n), Curve_tempr_bottom%y(Curves(i)%n))
            call dcsiez(Curves(i+1)%n, Curves(i+1)%t, Curves(i+1)%x, Curve_tempr_top%n, Curves(i)%t, Curve_tempr_top%x)
            call dcsiez(Curves(i+1)%n, Curves(i+1)%t, Curves(i+1)%y, Curve_tempr_top%n, Curves(i)%t, Curve_tempr_top%y)
            call dcsiez(Curves(i-1)%n, Curves(i-1)%t, Curves(i-1)%x, Curve_tempr_bottom%n, Curves(i)%t, Curve_tempr_bottom%x)
            call dcsiez(Curves(i-1)%n, Curves(i-1)%t, Curves(i-1)%y, Curve_tempr_bottom%n, Curves(i)%t, Curve_tempr_bottom%y)
            do j = 1, Curves(i)%n
                if (j == 1) then
                    p1 = cmplx(Curve_tempr_bottom%x(j),Curve_tempr_bottom%y(j))
                    p2 = cmplx(Curve_tempr_bottom%x(j+1),Curve_tempr_bottom%y(j+1))
                    p3 = cmplx(Curve_tempr_top%x(j+1),Curve_tempr_top%y(j+1))
                    p4 = cmplx(Curve_tempr_top%x(j),Curve_tempr_top%y(j))
                    V_0 = 2*area_quadrilateral(p1,p2,p3,p4)
                    Curves(i)%c(j)=d1
                elseif ( j == Curves(i)%n ) then
                    p1 = cmplx(Curve_tempr_bottom%x(j-1),Curve_tempr_bottom%y(j-1))
                    p2 = cmplx(Curve_tempr_bottom%x(j),Curve_tempr_bottom%y(j))
                    p3 = cmplx(Curve_tempr_top%x(j),Curve_tempr_top%y(j))
                    p4 = cmplx(Curve_tempr_top%x(j-1),Curve_tempr_top%y(j-1))
                    Curves(i)%c(j) = V_0/(2*area_quadrilateral(p1,p2,p3,p4))
                else
                    p1 = cmplx(Curve_tempr_bottom%x(j-1),Curve_tempr_bottom%y(j-1))
                    p2 = cmplx(Curve_tempr_bottom%x(j+1),Curve_tempr_bottom%y(j+1))
                    p3 = cmplx(Curve_tempr_top%x(j+1),Curve_tempr_top%y(j+1))
                    p4 = cmplx(Curve_tempr_top%x(j-1),Curve_tempr_top%y(j-1))
                    Curves(i)%c(j) = V_0/(area_quadrilateral(p1,p2,p3,p4))
                end if
            end do
            call finalize(Curve_tempr_top)
            call finalize(Curve_tempr_bottom)
            end if
        end do
    end subroutine
subroutine fcn_s_t(n, t, y, yprime) !интегрирование по длине дуги s с поиском времени
    use mod
    integer(4) :: n
    real(8) t, y(n), yprime(n), u_x, u_y, pg_get_fun_xy, V
    V           = dsqrt(y(3)**2+y(4)**2)
    u_x         = pg_get_fun_xy(y(1),y(2),2,d0,d1,0)
    u_y         = -pg_get_fun_xy(y(1),y(2),2,d1,d0,0)
    yprime(1)   = y(3)/V
    yprime(2)   = y(4)/V
    yprime(3)   = (u_x-y(3))/(st*V)
    yprime(4)   = (u_y-y(4))/(st*V)
    yprime(5)   = 1/V
    end subroutine fcn_s_t
subroutine fcn_s_t_top(n, t, y, yprime) !интегрирование по длине дуги s с поиском времени для верхней границы
    use mod
    integer(4) :: n
    real(8) t, y(n), yprime(n), u_x, pg_get_fun_xy, V
    V           = abs(y(3))
    u_x         = pg_get_fun_xy(y(1),y(2)-eps,2,d0,d1,0)
    yprime(1)   = y(3)/V
    yprime(2)   = d0
    yprime(3)   = (u_x-y(3))/(st*V)
    yprime(4)   = d0
    yprime(5)   = 1/V
    end subroutine fcn_s_t_top
subroutine fcn_s_t_bottom(n, t, y, yprime) !интегрирование по длине дуги s с поиском времени для нижней границы
    use mod
    integer(4) :: n
    real(8) t, y(n), yprime(n), u_x, pg_get_fun_xy, V
    V           = abs(y(3))
    u_x         = pg_get_fun_xy(y(1),y(2)+eps,2,d0,d1,0)
    yprime(1)   = y(3)/V
    yprime(2)   = d0
    yprime(3)   = (u_x-y(3))/(st*V)
    yprime(4)   = d0
    yprime(5)   = 1/V
    end subroutine fcn_s_t_bottom
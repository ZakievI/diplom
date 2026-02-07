subroutine solve !построение графика числа T = N_impact/N от числа стокса
    use mod
    implicit none
    integer(4), parameter::n = 4
    integer(4):: i, num, f_t(3)
    integer(4) function_impact_testing
    real(8), allocatable :: y_out(:,:) 
    real(8):: y(n), stocs, y_t(3)
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
    real(8):: y(n), stocs, y_t(3)
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
function function_impact_testing(n, y, dlt1)  !индикаторная функция
    use mod
    integer(4)::n
    real(8) y(n), tol, dlt1
    integer(4) function_impact_testing
    function_impact_testing = 2
    tol = 0.001d0
    if (y(2)<d1+dlt1-tol) then
        function_impact_testing = 1
    else if ((y(2)>= d1+dlt-tol).and.(y(2)<= d1+dlt1+tol)) then
        function_impact_testing = 0
    else if (y(2)> d1+dlt1+tol) then
        function_impact_testing = -1
    end if
    end function function_impact_testing
subroutine impact_test(n, y , y_out, num, dlt1) !проверка на ударение чатcицы
    use mod
    integer(4) :: n, num, k1, ido
    integer(4), parameter::mxparm = 100
    real(8) :: y_out(num,2)
    real(8) :: y(n), s, d_s, tol, param(mxparm), dlt1
    external fcn, fcn_s
    s           = d0                !��������� ������� ��� �������������� �� ����
    tol         = 0.001d0           !���������� ������
    param       = d0                !��������� �� ���������
    param(4)    = num               !������������ ���-�� ��������
    param(10)   = 1.0d0             !������������ ���������� �������
    d_s         = 0.05d0            !��� �� ����
    k1          = 2
    ido = 1
    do while ((sqrt(y(1)**2+y(2)**2)>d1+dlt1).and.(-L1/2<=y(1)).and.(y(1)<d0).and.(y(2)<H1).and.(d0<y(2)).and.(k1<=num))
        call divprk(ido, n, fcn_s, s, s+d_s, tol, param, y)
        !print '(i6, 6f12.3)', k1, t, y
        !y_out(k1,1)=y(1)
        !y_out(k1,2)=y(2)
        k1 = k1 + 1
    end do
    call divprk(3, n, fcn_s, s, s+d_s, tol, param, y)
    end subroutine impact_test
subroutine coordinate_first_particle(n, y, arr_bound, num, dlt1, size_arr) !координаты экстремальной частицы, прошедшей мимо цилиндрва  
    use mod 
    integer(4) :: n, num, k1, ido, size_arr
    integer(4), parameter::mxparm = 100
    real(8) :: arr_bound(num,2)
    real(8) :: y(n), s, d_s, tol, param(mxparm), dlt1
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
    real(8) :: y(n), y_t(3)
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
subroutine draw_Curves(par) !вывод кривых
    use mod
    !par=
    !   1-вывод траектории с концертацией и с элементами якобиана
    !   2-вывод траектории
    integer(4) :: i, l
    integer(4) :: par
    open (1, file='traektorie.dat')
    write(1,*) 'title = "traektorie"'
    select case (par)
    case (1)
        write(1,*) 'variables = "x", "y", "s", "consetr", "J_11", "J_12", "J_21", "J_22"'
    case (2)
        write(1,*) 'variables = "x", "y"'
    case DEFAULT
        write(1,*) 'variables = "x", "y"'
    end select
    do i = 1, size(Curves)
        write(1,"('ZONE T=""area',i0,'"", I=', i0, ', F=POINT')") i, Curves(i)%n
        do l = 1, Curves(i)%n
            select case (par)
            case (1)
                write(1,"(E15.5, ' ', E15.5, ' ', E15.5, ' ', E15.5, ' ', E15.5, ' ', E15.5, ' ', E15.5, ' ', E15.5)") &
                Curves(i)%x(l), Curves(i)%y(l), Curves(i)%s(l), Curves(i)%c(l), Curves(i)%J_11(l), Curves(i)%J_12(l), Curves(i)%J_21(l), Curves(i)%J_22(l)
            case (2)
                write(1,"(E15.5, ' ', E15.5)") &
                Curves(i)%x(l), Curves(i)%y(l)
            case DEFAULT
                write(1,"(E15.5, ' ', E15.5)") &
                Curves(i)%x(l), Curves(i)%y(l)
            end select
        end do
    end do
    end subroutine
subroutine build_time_isolines() !строим изолинии по времени
    use mod
    integer(4) :: i, j, nn
    integer(4), allocatable :: index_curves(:), index_curves_temporary(:)
    !complex(8), allocatable :: cord(:)
    real(8) :: t_max, x, y, c
    real(8), allocatable :: t_arr(:)
    nn = 100
    t_max = d0
    allocate(index_curves(size(Curves)))
    do i = 1, size(Curves)
        if (t_max < Curves(i)%t(size(Curves(i)%t))) then
            t_max = Curves(i)%t(size(Curves(i)%t))
        end if 
    end do 
    allocate(t_arr(nn))
    do i = 1, nn
        t_arr(i) = (t_max - d0) * (i - 1) / (nn - 1) + d0
    end do  
    open (3, file='time_isolines.dat')
    write(3,*) 'title = "time_isolines"'
    write(3,*) 'variables = "x", "y", "t", "c"'
    do i = 1, size(Curves)
        index_curves(i) = i         
    end do
    do i = 1, nn
        j = 1
        do 
            if (j > size(index_curves)) exit
            if (t_arr(i) > Curves(index_curves(j))%t(size(Curves(index_curves(j))%t))) then
                allocate(index_curves_temporary(size(index_curves) - 1))
                if (j == 1) then
                    index_curves_temporary(j:) = index_curves(j + 1:)
                elseif (j == size(index_curves)) then
                    index_curves_temporary(:j - 1) = index_curves(:j - 1)
                else
                    index_curves_temporary(:j - 1) = index_curves(:j - 1)
                    index_curves_temporary(j:) = index_curves(j + 1:)
                endif
                deallocate(index_curves)
                allocate(index_curves(size(index_curves_temporary)))
                index_curves = index_curves_temporary
                deallocate(index_curves_temporary)
            end if
            j = j + 1
        end do
        write(3,"('ZONE T=""area',i0,'"", I=', i0, ', F=POINT')") i, size(index_curves)
        do j = 1, size(index_curves)
            call dcsiez(size(Curves(index_curves(j))%t), Curves(index_curves(j))%t,Curves(index_curves(j))%x, 1, t_arr(i), x)
            call dcsiez(size(Curves(index_curves(j))%t), Curves(index_curves(j))%t,Curves(index_curves(j))%y, 1, t_arr(i), y)
            call dcsiez(size(Curves(index_curves(j))%t), Curves(index_curves(j))%t,Curves(index_curves(j))%c, 1, t_arr(i), c)
            !if ((-L1/2 > x > L1/2) .and. (d0 > y > H1) .and. ((x**2 + y**2) > 1)) then
            write(3,"(E15.5, ' ', E15.5, ' ', E15.5, ' ', E15.5)") x, y, t_arr(i), c
            !end if
        end do
    end do 
    deallocate(index_curves, t_arr)
    
    end subroutine 
subroutine build_curve() !поиск кривых
    use mod
    integer(4) :: i, n1, ido
    integer(4), parameter :: n = 4
    real(8) area_quadrilateral, search_for_extreme_particles, alfa
    integer(4), parameter :: mxparm = 50
    real(8) :: param(mxparm), d_s, s, y(n), tol, pg_get_fun_xy
    real(8), allocatable :: Curve_tempr(:,:)
    external fcn_s_1, area_quadrilateral, fcn_s_top_bottom, search_for_extreme_particles
    tol             = 1d-4
    d_s             = 1d-2
    s               = d0
    ido             = 1
    dlt             = d0
    allocate(Curves(num_particle))
    !cord_extreme_particles = search_for_extreme_particles()
    
    !$omp parallel do if (use_parallel_build_cerves == 1) private(i, n1, ido, s, y, Curve_tempr, param)
    do i = 1, num_particle
        allocate(Curve_tempr(N_arr,5)) !Curve_tempr = [x, y, V_x, V_y, s]
        write(*,"('I=',i0)") i
        n1                 = 1
        ido                = 1
        s                  = d0
        Curve_tempr(n1, 1) = -L1/2 
        !if ((H1*((i)/(num_particle-d1))**3>cord_extreme_particles).and.(cord_extreme_particles>H1*((i-d1)/(num_particle-d1))**3)) then
        !    !cord_extreme_particles = search_for_extreme_particles()
        !    Curve_tempr(n1, 2) = cord_extreme_particles - eps
        !    index_extreme_particles = i
        !else
        !    !Curve_tempr(n1, 2) = H1*(i - d1)/(num_particle - d1)
        !    ! измененно !!!!!!!!!!
        !    Curve_tempr(n1, 2) = H1*((i-d1)/(num_particle-d1))**3
        !    ! Curve_tempr(n1, 2) = 0.5d0
        !end if

        y(1)            = -L1/2
        y(2)            = bottom_coordinat + (top_coordinat - bottom_coordinat)*((i-d1)/(num_particle-d1))**2
        call get_uxuy(y(1), y(2), y(3), y(4))
        ! if (((dabs(y(2)) < eps) .or. (dabs(y(2) - H1) < eps)) .and. (-H1 - y(1) < eps) ) then
        !     y(3)            = pg_get_fun_xy(y(1) + eps ,y(2),2,d0,d1,0)
        !     y(4)            = -pg_get_fun_xy(y(1) + eps,y(2) ,2,d1,d0,0)
        ! ! else if(((y(1) + 1d0) < eps) .and. (y2 < eps))then 
        ! !     y(3)            = pg_get_fun_xy(y(1) - eps,y(2),2,d0,d1,0)
        ! !     y(4)            = -pg_get_fun_xy(y(1) - eps,y(2),2,d1,d0,0)
        ! else 
        !     y(3)            = pg_get_fun_xy(y(1),y(2),2,d0,d1,0)
        !     y(4)            = -pg_get_fun_xy(y(1),y(2),2,d1,d0,0)
        ! end if
        !y(5)            = d0
        param           = d0
        param(4)        = N_arr
        param(10)       = 1.0d0
        Curve_tempr(n1, 1) = y(1)
        Curve_tempr(n1, 2) = y(2)
        Curve_tempr(n1, 3) = y(3)
        Curve_tempr(n1, 4) = y(4)
        !Curve_tempr(n1, 5) = y(5)
        Curve_tempr(n1, 5) = s
        do while((dsqrt(y(1)**2+y(2)**2)>d1+dlt+1d-10).and.(-L1/2<=y(1)).and.(y(1)<=L1/2).and.(y(2)<=H1).and.(d0<=y(2)).and.(n1<N_arr))
            !write(*,"('N=',i0)") n1
            if (dabs(Curve_tempr(n1, 2)) < eps)then
                call divprk(ido, n, fcn_s_top_bottom, s, s+d_s, tol, param, y)
                y(2) = d0
            elseif (dabs(Curve_tempr(n1, 2) - H1) < eps) then
                call divprk(ido, n, fcn_s_top_bottom, s, s+d_s, tol, param, y)
                y(2) = H1
            else
                call divprk(ido, n, fcn_s_1, s, s+d_s, tol, param, y)
            end if 
            n1 = n1 + 1
            if (y(3) < 0) then
                write(*,*) 'Error: negative value V_x'
                alfa = datan(y(2)/y(1))
                y(3) = d_s*dcos(alfa)
                y(4) = d_s*dsin(alfa)
            end if 
            Curve_tempr(n1, 1) = y(1)
            Curve_tempr(n1, 2) = y(2)
            Curve_tempr(n1, 3) = y(3)
            Curve_tempr(n1, 4) = y(4)
            !Curve_tempr(n1, 5) = y(5)
            Curve_tempr(n1, 5) = s
        end do
        allocate(Curves(i)%x(n1))
        allocate(Curves(i)%y(n1))
        !allocate(Curves(i)%t(n1))
        allocate(Curves(i)%V_x(n1))
        allocate(Curves(i)%V_y(n1))
        allocate(Curves(i)%s(n1))
        Curves(i)%x(1:n1) = Curve_tempr(1:n1, 1)
        Curves(i)%y(1:n1) = Curve_tempr(1:n1, 2)
        Curves(i)%V_x(1:n1) = Curve_tempr(1:n1, 3)
        Curves(i)%V_y(1:n1) = Curve_tempr(1:n1, 4)
        !Curves(i)%t(1:n1) = Curve_tempr(1:n1, 5)
        Curves(i)%s(1:n1) = Curve_tempr(1:n1, 5)
        Curves(i)%n = n1
        call divprk(3, n, fcn_s_1, s, s+d_s, tol, param, y)
        deallocate(Curve_tempr)
    end do 
    !$OMP END PARALLEL DO
    do i = 1, num_particle
        if (Curves(i)%x(Curves(i)%n)>d0) then
            index_extreme_particles = i - 1
            exit
        end if
        index_extreme_particles = num_particle
    end do
    end subroutine build_curve
subroutine find_derivative() !find the derivative du_i/dx_j of the curve 
    use mod
    integer(4) :: i, j
    real(8) du1ds, du2ds, cosfi, sinfi
    integer(4), parameter :: nn1 = 4
    real(8) ::  a(nn1, nn1), b(nn1), xx(nn1)
    real(8) pg_get_fun_xy
    real(8) :: u_x_1, u_x_2, u_y_1, u_y_2, u_x_0, u_y_0
    do i = 1, num_particle
        allocate(Curves(i)%du1dx1(Curves(i)%n), Curves(i)%du2dx1(Curves(i)%n), Curves(i)%du1dx2(Curves(i)%n), Curves(i)%du2dx2(Curves(i)%n))
        do j = 1, Curves(i)%n
            if ((j > 1) .and.(j < Curves(i)%n)) then
                cosfi = ((Curves(i)%x(j + 1) - Curves(i)%x(j - 1))/(Curves(i)%s(j + 1) - Curves(i)%s(j - 1)))
                sinfi = ((Curves(i)%y(j + 1) - Curves(i)%y(j - 1))/(Curves(i)%s(j + 1) - Curves(i)%s(j - 1)))
                call get_uxuy(Curves(i)%x(j - 1), Curves(i)%y(j - 1), u_x_1, u_y_1)
                call get_uxuy(Curves(i)%x(j + 1), Curves(i)%y(j + 1), u_x_2, u_y_2)
                du1ds =  (u_x_2 - u_x_1)/(Curves(i)%s(j + 1) - Curves(i)%s(j - 1))
                du2ds =  (u_y_2 - u_y_1)/(Curves(i)%s(j + 1) - Curves(i)%s(j - 1))
            elseif (j == 1) then
                cosfi = (3 * Curves(i)%x(j) - 4 * Curves(i)%x(j + 1) + Curves(i)%x(j + 2))/(Curves(i)%s(j + 2) - Curves(i)%s(j))
                sinfi = (3 * Curves(i)%y(j) - 4 * Curves(i)%y(j + 1) + Curves(i)%y(j + 2))/(Curves(i)%s(j + 2) - Curves(i)%s(j))
                call get_uxuy(Curves(i)%x(j + 1), Curves(i)%y(j + 1), u_x_1, u_y_1)
                call get_uxuy(Curves(i)%x(j + 2), Curves(i)%y(j + 2), u_x_2, u_y_2)
                call get_uxuy(Curves(i)%x(j), Curves(i)%y(j), u_x_0, u_y_0)
                du1ds = (3 * u_x_0 - 4 * u_x_1 + u_x_2)/(Curves(i)%s(j + 2) - Curves(i)%s(j))
                du2ds = (3 * u_y_0 - 4 * u_y_1 + u_y_2)/(Curves(i)%s(j + 2) - Curves(i)%s(j))
            elseif (j == Curves(i)%n) then
                cosfi = (3 * Curves(i)%x(j) - 4 * Curves(i)%x(j - 1) + Curves(i)%x(j - 2))/(Curves(i)%s(j) - Curves(i)%s(j - 2))
                sinfi = (3 * Curves(i)%y(j) - 4 * Curves(i)%y(j - 1) + Curves(i)%y(j - 2))/(Curves(i)%s(j) - Curves(i)%s(j - 2))
                call get_uxuy(Curves(i)%x(j - 1), Curves(i)%y(j - 1), u_x_1, u_y_1)
                call get_uxuy(Curves(i)%x(j), Curves(i)%y(j), u_x_2, u_y_2)
                call get_uxuy(Curves(i)%x(j - 2), Curves(i)%y(j - 2), u_x_0, u_y_0)
                du1ds = (3 * u_x_2 - 4 * u_x_1 + u_x_0)/(Curves(i)%s(j) - Curves(i)%s(j - 2))
                du2ds = (3 * u_y_2 - 4 * u_y_1 + u_y_0)/(Curves(i)%s(j) - Curves(i)%s(j - 2))
            end if
            a = reshape((/1d0, d0, d0, 1d0, &
                          d0, -1d0, 1d0, d0, &
                          cosfi, sinfi, d0, d0, &
                          d0, d0, cosfi, sinfi /), (/ 4, 4 /))
            b = [d0, -pg_get_fun_xy(Curves(i)%x(j),Curves(i)%y(j),3,d0,d0,2), du1ds, du2ds]
            call DLSLRG(nn1, a, nn1, b, 2, xx)
            Curves(i)%du1dx1(j) = xx(1)
            Curves(i)%du1dx2(j) = xx(2)
            Curves(i)%du2dx1(j) = xx(3)
            Curves(i)%du2dx2(j) = xx(4)
        end do
    end do
    end subroutine find_derivative
subroutine build_mesh_1 !строим сетку 
    use mod
    integer(4) :: i, l, k
    real(8), allocatable :: tempr_y_c(:), tempr_x_c(:), tempr_c(:), tempr_v_x(:), tempr_v_y(:)
    integer(4) :: index_point_rigth_bound_array_area_1(num_particle), index_point_rigth_bound_array_area_2(num_particle - index_extreme_particles)
    integer(4) :: num_of_partitions_by_x
    real(8) :: x_c, y_c
    integer(4) :: begin_index_z_m, begin_index_trm
    integer(4) :: start_index_2st_area
    N_part_2 = index_extreme_particles
    num_of_partitions_by_x = (N_part_1 + N_part_2 + N_part_3 + N_part_4)
    ! allocate(mesh%x_y_(num_particle * num_of_partitions_by_x, 2),mesh%t(num_particle * num_of_partitions_by_x, 1),mesh%c(num_particle * num_of_partitions_by_x, 1)&
    !         mesh%v(num_particle * num_of_partitions_by_x, 2),mesh%s(num_particle * num_of_partitions_by_x, 2))
    allocate(mesh)
    mesh%n_i = num_of_partitions_by_x
    mesh%n_j = num_particle
    !allocate(mesh%x_y_(num_of_partitions_by_x * num_particle, 2))
    ! область 1 - x меняеться от -L1 до -1, y - от 0 до H1
    ! область 2 - x меняеться от -1 до 0, y - от Y_Last до H1
    ! область 3 - x меняеться от 0 до 1, y - от Y_Last до H1
    ! область 4 - x меняеться от 1 до L1, y - от 0 до H1
    allocate(mesh%z_m((N_part_1 + N_part_2 + N_part_3 + N_part_4) * num_particle))
    allocate(mesh%c((N_part_1 + N_part_2 + N_part_3 + N_part_4) * num_particle))
    allocate(mesh%trm(4,(N_part_1 + N_part_2 + N_part_3 + N_part_4 - 4) * (num_particle - 1)))
    allocate(mesh%v_m((N_part_1 + N_part_2 + N_part_3 + N_part_4) * num_particle))
    !allocate(index_point_rigth_bound_array_area_1(num_particle))
    
    do i = 1, 3        
        if (i == 1) then
            ! область 1
            begin_index_z_m = 1
            begin_index_trm = 1
            allocate(tempr_y_c(N_part_1),tempr_x_c(N_part_1), tempr_c(N_part_1), tempr_v_x(N_part_1), tempr_v_y(N_part_1))
            do l = 1, N_part_1
                x_c = -L1/2 + (L1/2 - d1)*(l-1)/(N_part_1-1)
                tempr_x_c(l) = x_c
            end do
            do k = 1, num_particle
                call dcsiez(Curves(k)%n, Curves(k)%x, Curves(k)%y, N_part_1, tempr_x_c, tempr_y_c)
                call dcsiez(Curves(k)%n, Curves(k)%x, Curves(k)%c, N_part_1, tempr_x_c, tempr_c)
                call dcsiez(Curves(k)%n, Curves(k)%x, Curves(k)%V_x, N_part_1, tempr_x_c, tempr_v_x)
                call dcsiez(Curves(k)%n, Curves(k)%x, Curves(k)%V_y, N_part_1, tempr_x_c, tempr_v_y)
                mesh%z_m(1 + (k-1)*N_part_1: k * N_part_1) = cmplx(tempr_x_c,tempr_y_c)
                mesh%v_m(1 + (k-1)*N_part_1: k * N_part_1) = cmplx(tempr_v_x,tempr_v_y)
                mesh%c(1 + (k-1)*N_part_1: k * N_part_1) = tempr_c
                index_point_rigth_bound_array_area_1(k) = k * N_part_1
            end do
            do l = 1, num_particle - 1
                do k = 1, N_part_1 - 1
                    mesh%trm(1, k + (N_part_1 - 1) * (l-1)) = k + (N_part_1) * (l - 1)     
                    mesh%trm(2, k + (N_part_1 - 1) * (l-1)) = k + 1 + (N_part_1) * (l - 1)
                    mesh%trm(3, k + (N_part_1 - 1) * (l-1)) = k + 1 + (N_part_1) * (l)
                    mesh%trm(4, k + (N_part_1 - 1) * (l-1)) = k + (N_part_1) * (l) 
                end do 
            end do
        elseif (i == 2)then
            ! область 2
            begin_index_z_m = num_particle * N_part_1 + 1
            begin_index_trm = (num_particle - 1) * (N_part_1 - 1) + 1

            if (Curves(1)%y(Curves(1)%n) == d0) then 
                start_index_2st_area = 2
            else 
                start_index_2st_area = 1
            end if
            
            allocate(tempr_y_c(N_part_2 + 1 - start_index_2st_area + 1),&
                    tempr_x_c(N_part_2 + 1 - start_index_2st_area + 1),&
                    tempr_c(N_part_2 + 1 - start_index_2st_area + 1),&
                    tempr_v_x(N_part_2 + 1 - start_index_2st_area + 1),&
                    tempr_v_y(N_part_2 + 1 - start_index_2st_area + 1))
            do l = start_index_2st_area, N_part_2
                tempr_x_c(l - start_index_2st_area + 1) = Curves(l)%x(Curves(l)%n)
                call dcsiez(Curves(l)%n, Curves(l)%x, Curves(l)%y, l - start_index_2st_area + 1, tempr_x_c(1:l - start_index_2st_area + 1), tempr_y_c)
                call dcsiez(Curves(l)%n, Curves(l)%x, Curves(l)%c, l - start_index_2st_area + 1, tempr_x_c(1:l - start_index_2st_area + 1), tempr_c)
                call dcsiez(Curves(l)%n, Curves(l)%x, Curves(l)%V_x, l - start_index_2st_area + 1, tempr_x_c(1:l - start_index_2st_area + 1), tempr_v_x)
                call dcsiez(Curves(l)%n, Curves(l)%x, Curves(l)%V_y, l - start_index_2st_area + 1, tempr_x_c(1:l - start_index_2st_area + 1), tempr_v_y)
                mesh%z_m(begin_index_z_m : begin_index_z_m + l - start_index_2st_area + 1 - 1) = cmplx(tempr_x_c(1 : l - start_index_2st_area + 1),tempr_y_c(1 : l - start_index_2st_area + 1))
                mesh%v_m(begin_index_z_m : begin_index_z_m + l - start_index_2st_area + 1 - 1) = cmplx(tempr_v_x(1 : l - start_index_2st_area + 1),tempr_v_y(1 : l - start_index_2st_area + 1))
                mesh%c(begin_index_z_m : begin_index_z_m + l - start_index_2st_area + 1 - 1) = tempr_c
                begin_index_z_m = begin_index_z_m + l - start_index_2st_area + 1
            end do
            tempr_x_c(N_part_2 + 1 - start_index_2st_area + 1) = d0
            do k = N_part_2 + 1, num_particle
                call dcsiez(Curves(k)%n, Curves(k)%x, Curves(k)%y, N_part_2 + 1 - start_index_2st_area + 1, tempr_x_c, tempr_y_c)
                call dcsiez(Curves(k)%n, Curves(k)%x, Curves(k)%c, N_part_2 + 1 - start_index_2st_area + 1, tempr_x_c, tempr_c)
                call dcsiez(Curves(k)%n, Curves(k)%x, Curves(k)%V_x, N_part_2 + 1 - start_index_2st_area + 1, tempr_x_c, tempr_v_x)
                call dcsiez(Curves(k)%n, Curves(k)%x, Curves(k)%V_y, N_part_2 + 1 - start_index_2st_area + 1, tempr_x_c, tempr_v_y)
                mesh%z_m(begin_index_z_m : begin_index_z_m + N_part_2 - start_index_2st_area + 1) = cmplx(tempr_x_c,tempr_y_c)
                mesh%v_m(begin_index_z_m : begin_index_z_m + N_part_2 - start_index_2st_area + 1) = cmplx(tempr_v_x,tempr_v_y)
                mesh%c(begin_index_z_m: begin_index_z_m + N_part_2 - start_index_2st_area + 1) = tempr_c
                begin_index_z_m = begin_index_z_m + N_part_2 + 1 - start_index_2st_area + 1
                index_point_rigth_bound_array_area_2(k - N_part_2) = begin_index_z_m - 1
            end do 
            begin_index_z_m = num_particle * N_part_1 + 1
            do l = 1, N_part_2
                if (l == num_particle) then 
                    exit
                end if
                do k = 1, l + 1 - start_index_2st_area + 1
                    if (k == 1) then
                        mesh%trm(1, begin_index_trm + k - 1) = index_point_rigth_bound_array_area_1(l)
                        mesh%trm(4, begin_index_trm + k - 1) = index_point_rigth_bound_array_area_1(l + 1)
                    else
                        mesh%trm(1, begin_index_trm + k - 1) = begin_index_z_m + k - 2
                        mesh%trm(4, begin_index_trm + k - 1) = begin_index_z_m + k - 2 + l - start_index_2st_area + 1
                    end if
                    
                    if (k == l + 1 - start_index_2st_area + 1) then
                        mesh%trm(2, begin_index_trm + k - 1) = begin_index_z_m + k - 1 + l - start_index_2st_area + 1
                        mesh%trm(3, begin_index_trm + k - 1) = begin_index_z_m + k - 1 + l - start_index_2st_area + 1
                    else
                        mesh%trm(2, begin_index_trm + k - 1) = begin_index_z_m + k - 1
                        mesh%trm(3, begin_index_trm + k - 1) = begin_index_z_m + k - 1 + l - start_index_2st_area + 1
                    end if
                end do
                begin_index_trm = begin_index_trm + l + 1 - start_index_2st_area + 1
                begin_index_z_m = begin_index_z_m + l - start_index_2st_area + 1
            end do
            do l = N_part_2 + 1, num_particle - 1
                do k = 1, N_part_2 + 1 - start_index_2st_area + 1
                    if (k == 1) then
                        mesh%trm(1, begin_index_trm + k - 1) = index_point_rigth_bound_array_area_1(l)
                        mesh%trm(4, begin_index_trm + k - 1) = index_point_rigth_bound_array_area_1(l + 1)
                    else 
                        mesh%trm(1, begin_index_trm + k - 1) = begin_index_z_m + k - 2
                        mesh%trm(4, begin_index_trm + k - 1) = begin_index_z_m + k - 2 + N_part_2 + 1 - start_index_2st_area + 1
                    end if 
                    mesh%trm(2, begin_index_trm + k - 1) = begin_index_z_m + k - 1
                    mesh%trm(3, begin_index_trm + k - 1) = begin_index_z_m + k - 1 + N_part_2 + 1 - start_index_2st_area + 1
                end do 
                begin_index_trm = begin_index_trm + N_part_2 + 1 - start_index_2st_area + 1
                begin_index_z_m = begin_index_z_m + N_part_2 + 1 - start_index_2st_area + 1
                index_point_rigth_bound_array_area_2(l - N_part_2) = begin_index_z_m - 1
            end do
            begin_index_z_m = begin_index_z_m + N_part_2 + 1 - start_index_2st_area + 1
            index_point_rigth_bound_array_area_2(num_particle - index_extreme_particles) = begin_index_z_m - 1
        elseif (i == 3) then
            allocate(tempr_y_c(N_part_3 - 1),tempr_x_c(N_part_3 - 1), tempr_c(N_part_3 - 1), tempr_v_x(N_part_3 - 1), tempr_v_y(N_part_3 - 1))
            do l = 2, N_part_3
                x_c = (L1/2)*(l-1)/(N_part_3-1)
                tempr_x_c(l - 1) = x_c
            end do
            do l = 1, num_particle - index_extreme_particles
                call dcsiez(Curves(index_extreme_particles + l)%n, Curves(index_extreme_particles + l)%x, Curves(index_extreme_particles + l)%y, N_part_3 - 1, tempr_x_c, tempr_y_c)
                call dcsiez(Curves(index_extreme_particles + l)%n, Curves(index_extreme_particles + l)%x, Curves(index_extreme_particles + l)%c, N_part_3 - 1, tempr_x_c, tempr_c)
                call dcsiez(Curves(index_extreme_particles + l)%n, Curves(index_extreme_particles + l)%x, Curves(index_extreme_particles + l)%V_x, N_part_3 - 1, tempr_x_c, tempr_v_x)
                call dcsiez(Curves(index_extreme_particles + l)%n, Curves(index_extreme_particles + l)%x, Curves(index_extreme_particles + l)%V_y, N_part_3 - 1, tempr_x_c, tempr_v_y)
                mesh%z_m(begin_index_z_m: begin_index_z_m + N_part_3 - 2) = cmplx(tempr_x_c,tempr_y_c)
                mesh%V_m(begin_index_z_m: begin_index_z_m + N_part_3 - 2) = cmplx(tempr_v_x,tempr_v_y)
                mesh%c(begin_index_z_m: begin_index_z_m + N_part_3 - 2) = tempr_c
                if (l /= num_particle - index_extreme_particles) then
                    do k = 1, N_part_3 - 1
                        if (k == 1) then 
                            mesh%trm(1, begin_index_trm + k - 1) = index_point_rigth_bound_array_area_2(l)
                            mesh%trm(4, begin_index_trm + k - 1) = index_point_rigth_bound_array_area_2(l + 1)
                        else 
                            mesh%trm(1, begin_index_trm + k - 1) = begin_index_z_m + k - 2
                            mesh%trm(4, begin_index_trm + k - 1) = begin_index_z_m + k - 2 + N_part_3 - 1
                        end if 
                        mesh%trm(2, begin_index_trm + k - 1) = begin_index_z_m + k - 1
                        mesh%trm(3, begin_index_trm + k - 1) = begin_index_z_m + k - 1 + N_part_3 - 1
                    end do
                    begin_index_trm = begin_index_trm + N_part_3 - 1
                end if 
                begin_index_z_m = begin_index_z_m + N_part_3 - 1
            end do
        end if
        mesh%ntr = begin_index_trm - 1
        mesh%n = begin_index_z_m - 1
        ! end do
        deallocate(tempr_y_c,tempr_x_c,tempr_c, tempr_v_x, tempr_v_y)
        
    end do 
    end subroutine build_mesh_1
subroutine draw_mesh(par) !вывод сетки
    use mod
    !par=
    !   1-вывод сетки с концертацией
    !   2-вывод сетки
    !   3-вывод сетки с силой f в узлах
    integer(4) i,k,par
    character(200) formatstr
    OPEN (1,FILE='tr_mesh.dat')
    write(1,*) 'TITLE = "Triangle Mesh"'
    select case (par)
    case (1)
        write(1,*) 'VARIABLES = "X", "Y", "C"'
        formatstr="('ZONE T=""area_e"" N=', i0, ', E=', i0, ', F=FEPOINT, ET=QUADRILATERAL')"
        write(1,trim(formatstr)) mesh%n, mesh%ntr
        do i = 1, mesh%n
            write(1,"(F9.5, ' ', F9.5, ' ', F9.5)") dreal(mesh%z_m(i)), dimag(mesh%z_m(i)), mesh%c(i)
        end do
    case(2)
        write(1,*) 'VARIABLES = "X", "Y"'
        formatstr="('ZONE T=""area_e"" N=', i0, ', E=', i0, ', F=FEPOINT, ET=QUADRILATERAL')"
        write(1,trim(formatstr)) mesh%n, mesh%ntr
        do i = 1, mesh%n
            write(1,"(F9.5, ' ', F9.5)") dreal(mesh%z_m(i)), dimag(mesh%z_m(i))
        end do
    case(3)
        write(1,*) 'VARIABLES = "X", "Y", "F_x", "F_y"'
        formatstr="('ZONE T=""area_e"" N=', i0, ', E=', i0, ', F=FEPOINT, ET=QUADRILATERAL')"
        write(1,trim(formatstr)) mesh%n, mesh%ntr
        do i = 1, mesh%n
            write(1,"(F9.5, ' ', F9.5, ' ', F9.5, ' ', F9.5)") dreal(mesh%z_m(i)), dimag(mesh%z_m(i)), dreal(mesh%f_m(i)), dimag(mesh%f_m(i))
        end do
    end select
    do i = 1, mesh%ntr
        WRITE(1,"(4(' ',i0))") (mesh%trm(k,i),k=1,mesh%npe)
    end do 
    close(1) 
    end subroutine
subroutine draw_mesh_with_cell(par)
    use mod
    !par=
    !   1-вывод сетки с F
    integer(4) i,k,par
    character(200) formatstr
    OPEN (1,FILE='mesh_data.txt')
    write(1,*) 'NODES CELLS NODES_PER_CELL'
    write(1,"(I, ' ', I, ' ', I)") mesh%n, mesh%ntr, mesh%npe

    select case (par)
    case (1)
        write(1,*) 'NODE_COORDINATES'
        
        do i = 1, mesh%n
            write(1,"(F9.5, ' ', F9.5)") dreal(mesh%z_m(i)), dimag(mesh%z_m(i))
        end do
        
        write(1,*) 'CONNECTIVITY'
        
        do i = 1, mesh%ntr
            write(1,"(I, ' ', I, ' ', I, ' ', I)") mesh%trm(1, i), mesh%trm(2, i), mesh%trm(3, i), mesh%trm(4, i)
        end do
        
        write(1,*) 'CELL_VALUES'
        
        do i = 1, mesh%ntr
            write(1,"(F9.5)") mesh%F_trm(i)
        end do
    end select
    close(1) 
    end subroutine
subroutine filling_in_the_cells()
    use mod
    integer(4) i, k
    real(8) u_x, u_y, fi1, fi2, tetta, sum_integral, l, area_quadrilateral, S_area
    external area_quadrilateral
    dlt             = 0.01d0
    allocate(mesh%f_m(mesh%n), mesh%F_trm(mesh%ntr))
    do i=1,mesh%n
        call get_uxuy(REAL(mesh%z_m(i)), AIMAG(mesh%z_m(i)), u_x, u_y)
        mesh%f_m(i) = 3 * pi * mu * cmplx(u_x - REAL(mesh%v_m(i)), u_y - AIMAG(mesh%v_m(i))) * dlt * mesh%c(i)
    end do
    do i=1, mesh%ntr
        sum_integral = d0
        do k=1, 4
            if (mesh%trm(k, i) == mesh%trm(MODULO(k, 4) + 1, i)) exit
            fi1 = datan2(REAL(mesh%z_m(mesh%trm(k, i))),AIMAG(mesh%z_m(mesh%trm(k, i))))
            fi2 = datan2(REAL(mesh%z_m(mesh%trm(MODULO(k, 4) + 1, i))),AIMAG(mesh%z_m(mesh%trm(MODULO(k, 4) + 1, i))))
            tetta = datan2(REAL(mesh%z_m(mesh%trm(MODULO(k, 4) + 1, i)) - mesh%z_m(mesh%trm(k , i))), AIMAG(mesh%z_m(mesh%trm(MODULO(k, 4) + 1, i)) - mesh%z_m(mesh%trm(k , i))))
            l = abs(mesh%z_m(mesh%trm(MODULO(k, 4) + 1, i)) - mesh%z_m(mesh%trm(k, i)))
            sum_integral = sum_integral + (mesh%f_m(mesh%trm(k, i))*dcos(fi1 - tetta) + mesh%f_m(mesh%trm(MODULO(k, 4) + 1, i))*dcos(fi2 - tetta))*l/2
        end do
        S_area = area_quadrilateral(mesh%z_m(mesh%trm(1, i)), mesh%z_m(mesh%trm(2, i)),mesh%z_m(mesh%trm(3, i)),mesh%z_m(mesh%trm(4, i)))
        if (S_area < 1d-5) then 
            S_area = 1d-5
        end if
        mesh%F_trm(i) = sum_integral / S_area
    end do
    end subroutine
subroutine find_concentration() !поиск концентрации
    use mod
    integer(4) :: i, j
    real(8) V_0, area_quadrilateral 
    !real(8), allocatable :: value_(:), xvec(:), xdata(:), fdata(:)
    type(Curve) :: Curve_tempr_top, Curve_tempr_bottom
    complex(8) p1, p2, p3, p4
    external area_quadrilateral
    V_0 = d1
    do i = 1, num_particle
        current_Curve=>Curves(i)
        if ((i == 1).or.(i == index_extreme_particles + 1)) then
            Curve_tempr_top%n = current_Curve%n
            allocate(Curve_tempr_top%x(current_Curve%n), Curve_tempr_top%y(current_Curve%n))
            call dcsiez(Curves(i+1)%n, Curves(i+1)%t, Curves(i+1)%x, Curve_tempr_top%n, current_Curve%t, Curve_tempr_top%x)!построение сплайна 
            call dcsiez(Curves(i+1)%n, Curves(i+1)%t, Curves(i+1)%y, Curve_tempr_top%n, current_Curve%t, Curve_tempr_top%y)
            do j = 1, current_Curve%n
                if (j == 1) then
                    p1 = cmplx(current_Curve%x(j),current_Curve%y(j))
                    p2 = cmplx(current_Curve%x(j+1),current_Curve%y(j+1))
                    p3 = cmplx(Curve_tempr_top%x(j+1),Curve_tempr_top%y(j+1))
                    p4 = cmplx(Curve_tempr_top%x(j),Curve_tempr_top%y(j))
                    V_0 = 4*area_quadrilateral(p1,p2,p3,p4)
                    current_Curve%c(j)=d1
                elseif ( j == current_Curve%n ) then
                    p1 = cmplx(current_Curve%x(j-1),current_Curve%y(j-1))
                    p2 = cmplx(current_Curve%x(j),current_Curve%y(j))
                    p3 = cmplx(Curve_tempr_top%x(j),Curve_tempr_top%y(j))
                    p4 = cmplx(Curve_tempr_top%x(j-1),Curve_tempr_top%y(j-1))
                    current_Curve%c(j) = V_0/(4*area_quadrilateral(p1,p2,p3,p4))
                else
                    p1 = cmplx(current_Curve%x(j-1),current_Curve%y(j-1))
                    p2 = cmplx(current_Curve%x(j+1),current_Curve%y(j+1))
                    p3 = cmplx(Curve_tempr_top%x(j+1),Curve_tempr_top%y(j+1))
                    p4 = cmplx(Curve_tempr_top%x(j-1),Curve_tempr_top%y(j-1))
                    current_Curve%c(j) = V_0/(2*area_quadrilateral(p1,p2,p3,p4))
                end if
            end do
            call finalize(Curve_tempr_top)
        ! else if (i == index_extreme_particles): 
        else if ((i == num_particle) .or. (i == index_extreme_particles)) then 
            Curve_tempr_bottom%n = current_Curve%n
            allocate(Curve_tempr_bottom%x(current_Curve%n), Curve_tempr_bottom%y(current_Curve%n))
            call dcsiez(Curves(i-1)%n, Curves(i-1)%t, Curves(i-1)%x, Curve_tempr_bottom%n, current_Curve%t, Curve_tempr_bottom%x)
            call dcsiez(Curves(i-1)%n, Curves(i-1)%t, Curves(i-1)%y, Curve_tempr_bottom%n, current_Curve%t, Curve_tempr_bottom%y)
            do j = 1, current_Curve%n
                if (j == 1) then
                    p1 = cmplx(Curve_tempr_bottom%x(j),Curve_tempr_bottom%y(j))
                    p2 = cmplx(Curve_tempr_bottom%x(j+1),Curve_tempr_bottom%y(j+1))
                    p3 = cmplx(Curves(i)%x(j+1),Curves(i)%y(j+1))
                    p4 = cmplx(Curves(i)%x(j),Curves(i)%y(j))
                    V_0 = 4*area_quadrilateral(p1,p2,p3,p4)
                    Curves(i)%c(j)=d1
                elseif ( j == current_Curve%n ) then
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
subroutine fcn_s_1(n, s, y, yprime) !интегрирование по длине дуги s
    use mod
    integer(4) :: n
    real(8) s, y(n), yprime(n), u_x, u_y, V
    V           = dsqrt(y(3)**2+y(4)**2)
    call get_uxuy(y(1), y(2), u_x, u_y)
    yprime(1)   = y(3)/V
    yprime(2)   = y(4)/V
    yprime(3)   = (u_x-y(3))/(st*V)
    yprime(4)   = (u_y-y(4))/(st*V)
    end subroutine fcn_s_1
subroutine fcn_s_top_bottom(n, s, y, yprime) !интегрирование по длине дуги s для нижней и верхней границы
    use mod
    integer(4) :: n
    real(8) s, y(n), yprime(n), u_x, u_y, V, r_eps
    !real(8), parameter :: r_eps = 1d-3
    
    if (y(3) < d0) then
        ! alfa = datan2(y(2), y(1))
        y(3)    = 1d-4
    end if
    
    V           = y(3)
    r_eps       = 0.5d0
    call get_uxuy(y(1), y(2), u_x, u_y)
    if ((V<u_x).and.(y(1)>-d1-r_eps))then 
       V       = u_x
       y(3)    = u_x
    end if 

    yprime(1)   = d1
    yprime(2)   = d0
    yprime(3)   = (u_x-V)/(st*V)
    yprime(4)   = d0
    end subroutine fcn_s_top_bottom
subroutine find_Jacobian()
    use mod 
    integer(4) :: i, ido, j
    integer(4), parameter :: n = 8
    integer(4), parameter :: mxparm = 50
    real(8) :: y(n), s, d_s, tol, param(mxparm)
    external fcn_t_Jacobian
    write(*,*) 'find Jacobian'
    !!$omp parallel do private(i, t, y, param, current_Curve)
    do i = 1, num_particle
        write(*,"('I=',i0)") i
        allocate(Curves(i)%J_11(Curves(i)%n))
        allocate(Curves(i)%J_12(Curves(i)%n))
        allocate(Curves(i)%J_21(Curves(i)%n))
        allocate(Curves(i)%J_22(Curves(i)%n))
        allocate(Curves(i)%J(Curves(i)%n))
        current_Curve => Curves(i)
        ido         = 1
        j           = 1
        tol         = 1d-3
        s           = d0
        y           = d0 ! y = [J11, J12, J21, J22, w11, w12, w21, w22]
        y(1)        = d1
        y(4)        = d1
        param       = d0
        !param(4)    = current_Curve%n+10
        param(4)    = N_arr
        param(10)   = 1.0d0
        current_Curve%J_11(j)   = y(1)
        current_Curve%J_12(j)   = y(2)
        current_Curve%J_21(j)   = y(3)
        current_Curve%J_22(j)   = y(4)
        current_Curve%J(j)      = dabs(y(1)*y(4) - y(2)*y(3))
        do while (j < current_Curve%n)
            !write(*,"('N=',i0)") j
            d_s = current_Curve%s(j + 1) - current_Curve%s(j)
            call divprk(ido, n, fcn_t_Jacobian, s, s + d_s, tol, param, y)
            j = j + 1
            current_Curve%J_11(j)   = y(1)
            current_Curve%J_12(j)   = y(2)
            current_Curve%J_21(j)   = y(3)
            current_Curve%J_22(j)   = y(4)
            current_Curve%J(j)      = dabs(y(1)*y(4) - y(2)*y(3))
        end do
        call divprk(3, n, fcn_t_Jacobian, current_Curve%s(j - 1), current_Curve%s(j), tol, param, y)
    end do
    !!$OMP END PARALLEL DO
    end subroutine find_Jacobian
subroutine fcn_t_Jacobian(n, s, y, yprime)
    use mod 
    integer(4) :: n 
    real(8) s, y(n), yprime(n), fcn_derivative_for_u_ij
    real(8) :: V, V_x, V_y
    external fcn_derivative_for_u_ij
    call dcsiez(current_Curve%n, current_Curve%s, current_Curve%V_x, 1, s, V_x)
    call dcsiez(current_Curve%n, current_Curve%s, current_Curve%V_y, 1, s, V_y)
    V = dsqrt(V_x**2+V_y**2)
    ! y = [J11, J12, J21, J22, w11, w12, w21, w22]
    ! yprime = [d_J11, d_J12, d_J21, d_J22, d_w11, d_w12, d_w21, d_w22]
    yprime(1) = y(5)/V
    yprime(2) = y(6)/V
    yprime(3) = y(7)/V
    yprime(4) = y(8)/V
    yprime(5) = (y(1)*fcn_derivative_for_u_ij(1, 1, s) + y(3)*fcn_derivative_for_u_ij(1, 2, s) - y(5))/st/V
    yprime(6) = (y(2)*fcn_derivative_for_u_ij(1, 1, s) + y(4)*fcn_derivative_for_u_ij(1, 2, s) - y(6))/st/V
    yprime(7) = (y(1)*fcn_derivative_for_u_ij(2, 1, s) + y(3)*fcn_derivative_for_u_ij(2, 2, s) - y(7))/st/V
    yprime(8) = (y(2)*fcn_derivative_for_u_ij(2, 1, s) + y(4)*fcn_derivative_for_u_ij(2, 2, s) - y(8))/st/V
    end subroutine fcn_t_Jacobian
function fcn_derivative_for_u_ij(i, j, s)
    use mod
    integer(4) :: i, j
    real(8) :: s, fcn_derivative_for_u_ij
    if ((s >= current_Curve%s(1)) .and. (s<= current_Curve%s(current_Curve%n))) then 
        if (i == 1) then
            if (j == 1) then
                call dcsiez(current_Curve%n, current_Curve%s, current_Curve%du1dx1, 1, s, fcn_derivative_for_u_ij)
            else
                call dcsiez(current_Curve%n, current_Curve%s, current_Curve%du1dx2, 1, s, fcn_derivative_for_u_ij)
            end if
        else
            if (j == 1) then
                call dcsiez(current_Curve%n, current_Curve%s, current_Curve%du2dx1, 1, s, fcn_derivative_for_u_ij)
            else
                call dcsiez(current_Curve%n, current_Curve%s, current_Curve%du2dx2, 1, s, fcn_derivative_for_u_ij)
            end if
        end if 
    else
        write(*,*) 'ERROR: fcn_derivative_for_u_ij'
    end if 
    end function fcn_derivative_for_u_ij
subroutine find_concentration_by_Jacobian()
    use mod 
    integer(4) :: i, k
    do i = 1, size(Curves)
        allocate(Curves(i)%c(Curves(i)%n))
        do k = 1, Curves(i)%n
            Curves(i)%c(k) = 1/Curves(i)%J(k)
            !if (Curves(i)%c(k) > 1) then
            !    Curves(i)%c(k) = 2d0 - 1d0/Curves(i)%c(k)
            !end if 
        end do 
    end do
    end subroutine find_concentration_by_Jacobian
subroutine checking_derivatives() ! проверка производных du_i/dx_j при помощи апроксимации
    use mod
    integer(4) :: i, j
    real(8) :: h, du1dx1, du1dx2, du2dx1, du2dx2
    real(8) :: pg_get_fun_xy
    h = 0.00001d0
    do j = 1, num_particle
        current_Curve=>Curves(j)
        allocate(current_Curve%discrepancy_du1dx1(current_Curve%n))
        allocate(current_Curve%discrepancy_du1dx2(current_Curve%n))
        allocate(current_Curve%discrepancy_du2dx1(current_Curve%n))
        allocate(current_Curve%discrepancy_du2dx2(current_Curve%n))
        allocate(current_Curve%aprox_du1dx1(current_Curve%n))
        allocate(current_Curve%aprox_du2dx1(current_Curve%n))
        allocate(current_Curve%aprox_du1dx2(current_Curve%n))
        allocate(current_Curve%aprox_du2dx2(current_Curve%n))
        do i = 1, current_Curve%n
            if (i == 1) then
                du1dx1 = (pg_get_fun_xy(current_Curve%x(i) + h,current_Curve%y(i),2,d0,d1,0) - &
                    4*pg_get_fun_xy(current_Curve%x(i) + h/2,current_Curve%y(i),2,d0,d1,0) + &
                    3*pg_get_fun_xy(current_Curve%x(i),current_Curve%y(i),2,d0,d1,0))/(h)
                du2dx1 = (pg_get_fun_xy(current_Curve%x(i) + h,current_Curve%y(i),2,d1,d0,0) - &
                    4*pg_get_fun_xy(current_Curve%x(i) + h/2,current_Curve%y(i),2,d1,d0,0) + &
                    3*pg_get_fun_xy(current_Curve%x(i),current_Curve%y(i),2,d1,d0,0))/(h)
            elseif (i == current_Curve%n) then
                du1dx1 = (3*pg_get_fun_xy(current_Curve%x(i),current_Curve%y(i),2,d0,d1,0) - &
                    4*pg_get_fun_xy(current_Curve%x(i) - h/2,current_Curve%y(i),2,d0,d1,0) + &
                    pg_get_fun_xy(current_Curve%x(i) - h,current_Curve%y(i),2,d0,d1,0))/(h)
                du2dx1 = (3*pg_get_fun_xy(current_Curve%x(i),current_Curve%y(i),2,d1,d0,0) - &
                    4*pg_get_fun_xy(current_Curve%x(i) - h/2,current_Curve%y(i),2,d1,d0,0) + &
                    pg_get_fun_xy(current_Curve%x(i) - h,current_Curve%y(i),2,d1,d0,0))/(h)
            else 
                du1dx1 = (pg_get_fun_xy(current_Curve%x(i) + h/2,current_Curve%y(i),2,d0,d1,0) - &
                    pg_get_fun_xy(current_Curve%x(i) - h/2,current_Curve%y(i),2,d0,d1,0))/(h)
                du2dx1 = (-pg_get_fun_xy(current_Curve%x(i) + h/2,current_Curve%y(i),2,d1,d0,0) + &
                    pg_get_fun_xy(current_Curve%x(i) - h/2,current_Curve%y(i),2,d1,d0,0))/(h)
            end if
            if (current_Curve%x(i)**2 + (current_Curve%y(i) + h/2)**2 < d1) then 
                du1dx2 = (pg_get_fun_xy(current_Curve%x(i),current_Curve%y(i) + h/2,2,d0,d1,0) - &
                    pg_get_fun_xy(current_Curve%x(i),current_Curve%y(i) - h/2,2,d0,d1,0))/h
                du2dx2 = (-pg_get_fun_xy(current_Curve%x(i),current_Curve%y(i) + h/2,2,d1,d0,0) + &
                    pg_get_fun_xy(current_Curve%x(i),current_Curve%y(i) - h/2,2,d1,d0,0))/h
            else 
                du1dx2 = (pg_get_fun_xy(current_Curve%x(i),current_Curve%y(i) + h/2,2,d0,d1,0) - &
                    pg_get_fun_xy(current_Curve%x(i),current_Curve%y(i),2,d0,d1,0))/(0.5d0*h)
                du2dx2 = (-pg_get_fun_xy(current_Curve%x(i),current_Curve%y(i) + h/2,2,d1,d0,0) + &
                    pg_get_fun_xy(current_Curve%x(i),current_Curve%y(i),2,d1,d0,0))/(0.5d0*h)
            end if 
            current_Curve%discrepancy_du1dx1(i) = dabs(current_Curve%du1dx1(i) - du1dx1)
            current_Curve%discrepancy_du2dx1(i) = dabs(current_Curve%du2dx1(i) - du2dx1)
            current_Curve%discrepancy_du1dx2(i) = dabs(current_Curve%du1dx2(i) - du1dx2)
            current_Curve%discrepancy_du2dx2(i) = dabs(current_Curve%du2dx2(i) - du2dx2)
            current_Curve%aprox_du1dx1(i) = du1dx1
            current_Curve%aprox_du1dx2(i) = du1dx2
            current_Curve%aprox_du2dx1(i) = du2dx1
            current_Curve%aprox_du2dx2(i) = du2dx2
        end do
    end do
    end subroutine checking_derivatives  

subroutine get_uxuy(x,y,ux,uy)
    use mod 
    real(8) rr,tt,ux,uy,x,y,ff1(2),vr,vt
    complex(8) z1,z2,z3,z4,z5,z6,z
    logical solved
    !заплатка на вычисление скорости вблизи угла на цилиндре
    rr=d1+d_zapl
    solved=.false.
    z1=dcmplx(-d1,d0)
    z2=dcmplx(d1,d0)
    z3=dcmplx(-L1/2,d0)
    z4=dcmplx(-L1/2,H1)
    z5=dcmplx(L1/2,d0)
    z6=dcmplx(L1/2,H1)
    z=dcmplx(x,y)
    if (cdabs(z-z1)<d_zapl .or. cdabs(z-z2)<d_zapl) then
        call rr_tt(x,y,rr,tt)
        vr=(cc1_zapl*(rr-d1/rr)**2+dd1_zapl*d5*(-d1+d1/rr**2+d2*dlog(rr)))*dcos(tt)
        vt=(cc1_zapl*(d2+d1/rr**2-3*rr**2)-dd1_zapl*d5*(d1-d1/rr**2+d2*dlog(rr)))*dsin(tt)
        if (zaplat==3) then
            vr=vr+(cc2_zapl*(d1/rr**4-4*rr**2+3*rr**4)-dd2_zapl*(d2/rr**4-3/rr**2+rr**2))*dcos(3*tt)
            vt=vt+(cc2_zapl*(d1/rr**4+4*rr**2-5*rr**4)+dd2_zapl*(-d2/rr**4+d1/rr**2+rr**2))*dsin(3*tt)
        endif
        call uvrt_to_xy(tt,vr,vt,ux,uy)
    elseif (cdabs(z-z3)<ds_pg .or. cdabs(z-z4)<ds_pg .or. cdabs(z-z5)<ds_pg .or. cdabs(z-z6)<ds_pg) then 
        ux = d1
        uy = d0
    else
        call pg_get_fun_xy_gradient(x,y,5,0,ff1)
        ux=ff1(2)
        uy=-ff1(1)
    endif
    end

subroutine init_zaplat
    use mod
    real(8) rr,pg_get_fun_xy,ux,uy,tt,x,y,ff1(2)
    real(8) m2(4,4),b2(4),res2(4)
    d_zapl=3*ds_pg
    rr=d1+d_zapl
    tt=pi
    !vr
    m2(1,1)=(rr-d1/rr)**2*dcos(tt)
    m2(1,2)=d5*(-d1+d1/rr**2+d2*dlog(rr))*dcos(tt)
    m2(1,3)=(d1/rr**4-4*rr**2+3*rr**4)*dcos(3*tt)
    m2(1,4)=-(d2/rr**4-3/rr**2+rr**2)*dcos(3*tt)
    b2(1)=-pg_get_fun_xy(-rr,d0,2,d0,d1,0)
    !!dvr/dr
    x=-d1-d_zapl*dcos(pi/4)
    y=d_zapl*dsin(pi/4)
    call rr_tt(x,y,rr,tt)
    !om
    m2(2,1)=-8*rr*dsin(tt)
    m2(2,2)=-2/rr*dsin(tt)
    m2(2,3)=-16*rr**3*dsin(3*tt)
    m2(2,4)=8/rr**3*dsin(3*tt)
    b2(2)=-pg_get_fun_xy(x,y,3,d0,d0,0)
    !vr
    m2(3,1)=(rr-d1/rr)**2*dcos(tt)
    m2(3,2)=d5*(-d1+d1/rr**2+d2*dlog(rr))*dcos(tt)
    m2(3,3)=(d1/rr**4-4*rr**2+3*rr**4)*dcos(3*tt)
    m2(3,4)=-(d2/rr**4-3/rr**2+rr**2)*dcos(3*tt)
    !vt
    m2(4,1)=(d2+d1/rr**2-3*rr**2)*dsin(tt)
    m2(4,2)=-d5*(d1-d1/rr**2+d2*dlog(rr))*dsin(tt)
    m2(4,3)=(d1/rr**4+4*rr**2-5*rr**4)*dsin(3*tt)
    m2(4,4)=(-d2/rr**4+d1/rr**2+rr**2)*dsin(3*tt)
    call pg_get_fun_xy_gradient(x,y,5,0,ff1)
    ux=ff1(2)
    uy=-ff1(1)
    b2(3)=ux*dcos(tt)+uy*dsin(tt) !vr
    b2(4)=-ux*dsin(tt)+uy*dcos(tt) !vt
    call DLSLRG(4,m2,4,b2,1,res2)
    cc1_zapl=res2(1)
    dd1_zapl=res2(2)
    cc2_zapl=res2(3)
    dd2_zapl=res2(4)
    end

subroutine uvrt_to_xy(tt,vr,vtt,u,v)
    !перевод скоростей из полярной в декартову
    use mod
    real(8) tt,vr,vtt,u,v
    u=vr*dcos(tt)-vtt*dsin(tt)
    v=vr*dsin(tt)+vtt*dcos(tt)
    end
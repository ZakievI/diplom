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
        write(1,*) 'variables = "x", "y", "t", "consetr", "J_11", "J_12", "J_21", "J_22"'
    case (2)
        write(1,*) 'variables = "x", "y"'
    case DEFAULT
        write(1,*) 'variables = "x", "y", "t"'
    end select
    do i = 1, size(Curves)
        write(1,"('ZONE T=""area',i0,'"", I=', i0, ', F=POINT')") i, Curves(i)%n
        do l = 1, Curves(i)%n
            select case (par)
            case (1)
                write(1,"(E15.5, ' ', E15.5, ' ', E15.5, ' ', E15.5, ' ', E15.5, ' ', E15.5, ' ', E15.5, ' ', E15.5)") &
                Curves(i)%x(l), Curves(i)%y(l), Curves(i)%t(l), Curves(i)%c(l), Curves(i)%J_11(l), Curves(i)%J_12(l), Curves(i)%J_21(l), Curves(i)%J_22(l)
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
    integer(4), parameter :: n = 5
    real(8) area_quadrilateral, search_for_extreme_particles
    integer(4), parameter :: mxparm = 50
    real(8) :: param(mxparm), d_s, s, y(n), dlt, tol, pg_get_fun_xy
    real(8), allocatable :: Curve_tempr(:,:)
    external fcn_s_t, area_quadrilateral, fcn_s_t_top, fcn_s_t_bottom, search_for_extreme_particles
    tol             = 0.0001d0
    d_s             = 0.01d0
    s               = d0
    ido             = 1
    dlt             = d0
    allocate(Curves(num_particle))
    !cord_extreme_particles = search_for_extreme_particles()
    
    !$omp parallel do if (use_parallel_build_cerves == 1) private(i, n1, ido, s, y, Curve_tempr, param)
    do i = 1, num_particle
        allocate(Curve_tempr(N_arr,6)) !Curve_tempr = [x, y, V_x, V_y, t, s]
        write(*,"('I=',i0)") i
        n1                 = 1
        ido                = 1
        s                  = d0
        Curve_tempr(n1, 1) = -L1/2 
        ! if ((H1*((i)/(num_particle-d1))**3>cord_extreme_particles).and.(cord_extreme_particles>H1*((i-d1)/(num_particle-d1))**3)) then
        !     cord_extreme_particles = search_for_extreme_particles()
        !     Curve_tempr(n1, 2) = cord_extreme_particles - eps
        !     index_extreme_particles = i
        ! else
        !     !Curve_tempr(n1, 2) = H1*(i - d1)/(num_particle - d1)
        !     ! измененно !!!!!!!!!!
        !     Curve_tempr(n1, 2) = H1*((i-d1)/(num_particle-d1))**3
        !     ! Curve_tempr(n1, 2) = 0.5d0
        ! end if

        
        y(1)            = -L1/2
        y(2)            = bottom_coordinat + (top_coordinat - bottom_coordinat)*((i-d1)/(num_particle-d1))**1
        if (dabs(y(2)) < eps) then
            y(3)            = pg_get_fun_xy(y(1),y(2) + eps,2,d0,d1,0)
            y(4)            = -pg_get_fun_xy(y(1),y(2) + eps,2,d1,d0,0)
        else if(dabs(y(2) - H1) < eps)then 
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
        Curve_tempr(n1, 1) = y(1)
        Curve_tempr(n1, 2) = y(2)
        Curve_tempr(n1, 3) = y(3)
        Curve_tempr(n1, 4) = y(4)
        Curve_tempr(n1, 5) = y(5)
        Curve_tempr(n1, 6) = s
        do while((sqrt(y(1)**2+y(2)**2)>d1+dlt).and.(-L1/2<=y(1)).and.(y(1)<=L1/2).and.(y(2)<=H1).and.(d0<=y(2)).and.(n1<N_arr))
            ! write(*,"('N=',i0)") n1
            if (dabs(Curve_tempr(n1, 2)) < eps) then
                call divprk(ido, n, fcn_s_t_bottom, s, s+d_s, tol, param, y)
                y(2) = d0
            else if(dabs(Curve_tempr(n1, 2) - H1) < eps)then 
                call divprk(ido, n, fcn_s_t_top, s, s+d_s, tol, param, y)
                y(2) = H1
            else
                call divprk(ido, n, fcn_s_t, s, s+d_s, tol, param, y)
            end if 
            n1 = n1 + 1
            Curve_tempr(n1, 1) = y(1)
            Curve_tempr(n1, 2) = y(2)
            Curve_tempr(n1, 3) = y(3)
            Curve_tempr(n1, 4) = y(4)
            Curve_tempr(n1, 5) = y(5)
            Curve_tempr(n1, 6) = s
        end do
        allocate(Curves(i)%x(n1))
        allocate(Curves(i)%y(n1))
        allocate(Curves(i)%t(n1))
        allocate(Curves(i)%V_x(n1))
        allocate(Curves(i)%V_y(n1))
        allocate(Curves(i)%s(n1))
        Curves(i)%x(1:n1) = Curve_tempr(1:n1, 1)
        Curves(i)%y(1:n1) = Curve_tempr(1:n1, 2)
        Curves(i)%V_x(1:n1) = Curve_tempr(1:n1, 3)
        Curves(i)%V_y(1:n1) = Curve_tempr(1:n1, 4)
        Curves(i)%t(1:n1) = Curve_tempr(1:n1, 5)
        Curves(i)%s(1:n1) = Curve_tempr(1:n1, 6)
        Curves(i)%n = n1
        call divprk(3, n, fcn_s_t, s, s+d_s, tol, param, y)
        deallocate(Curve_tempr)
    end do 
    !$OMP END PARALLEL DO
    end subroutine build_curve
subroutine find_derivative() !find the derivative du_i/dx_j of the curve 
    use mod
    integer(4) :: i, j
    real(8) du1ds, du2ds, cosfi, sinfi
    integer(4), parameter :: nn1 = 4
    real(8) ::  a(nn1, nn1), b(nn1), xx(nn1)
    real(8) pg_get_fun_xy
    do i = 1, num_particle
        allocate(Curves(i)%du1dx1(Curves(i)%n), Curves(i)%du2dx1(Curves(i)%n), Curves(i)%du1dx2(Curves(i)%n), Curves(i)%du2dx2(Curves(i)%n))
        do j = 1, Curves(i)%n
            if ((j > 1) .and.(j < Curves(i)%n)) then
                cosfi = ((Curves(i)%x(j + 1) - Curves(i)%x(j - 1))/(Curves(i)%s(j + 1) - Curves(i)%s(j - 1)))
                sinfi = ((Curves(i)%y(j + 1) - Curves(i)%y(j - 1))/(Curves(i)%s(j + 1) - Curves(i)%s(j - 1)))
                du1ds =  (pg_get_fun_xy(Curves(i)%x(j + 1),Curves(i)%y(j + 1),2,d0,d1,0) - pg_get_fun_xy(Curves(i)%x(j - 1),Curves(i)%y(j - 1),2,d0,d1,0))/(Curves(i)%s(j + 1) - Curves(i)%s(j - 1))
                du2ds =  (-pg_get_fun_xy(Curves(i)%x(j + 1),Curves(i)%y(j + 1),2,d1,d0,0) + pg_get_fun_xy(Curves(i)%x(j - 1),Curves(i)%y(j - 1),2,d1,d0,0))/(Curves(i)%s(j + 1) - Curves(i)%s(j - 1))
            elseif (j == 1) then
                cosfi = (3 * Curves(i)%x(j) - 4 * Curves(i)%x(j + 1) + Curves(i)%x(j + 2))/(Curves(i)%s(j + 2) - Curves(i)%s(j))
                sinfi = (3 * Curves(i)%y(j) - 4 * Curves(i)%y(j + 1) + Curves(i)%y(j + 2))/(Curves(i)%s(j + 2) - Curves(i)%s(j))
                du1ds = (3 * pg_get_fun_xy(Curves(i)%x(j),Curves(i)%y(j),2,d0,d1,0) - 4 * pg_get_fun_xy(Curves(i)%x(j + 1),Curves(i)%y(j + 1),2,d0,d1,0)&
                        + pg_get_fun_xy(Curves(i)%x(j + 2),Curves(i)%y(j + 2),2,d0,d1,0))/(Curves(i)%s(j + 2) - Curves(i)%s(j))
                du2ds = (-3 * pg_get_fun_xy(Curves(i)%x(j),Curves(i)%y(j),2,d1,d0,0) + 4 * pg_get_fun_xy(Curves(i)%x(j + 1),Curves(i)%y(j + 1),2,d1,d0,0)&
                        - pg_get_fun_xy(Curves(i)%x(j + 2),Curves(i)%y(j + 2),2,d1,d0,0))/(Curves(i)%s(j + 2) - Curves(i)%s(j))
            elseif (j == Curves(i)%n) then
                cosfi = (3 * Curves(i)%x(j) - 4 * Curves(i)%x(j - 1) + Curves(i)%x(j - 2))/(Curves(i)%s(j) - Curves(i)%s(j - 2))
                sinfi = (3 * Curves(i)%y(j) - 4 * Curves(i)%y(j - 1) + Curves(i)%y(j - 2))/(Curves(i)%s(j) - Curves(i)%s(j - 2))
                du1ds = (3 * pg_get_fun_xy(Curves(i)%x(j),Curves(i)%y(j),2,d0,d1,0) - 4 * pg_get_fun_xy(Curves(i)%x(j - 1),Curves(i)%y(j - 1),2,d0,d1,0)&
                        + pg_get_fun_xy(Curves(i)%x(j - 2),Curves(i)%y(j - 2),2,d0,d1,0))/(Curves(i)%s(j) - Curves(i)%s(j - 2))
                du2ds = (-3 * pg_get_fun_xy(Curves(i)%x(j),Curves(i)%y(j),2,d1,d0,0) + 4 * pg_get_fun_xy(Curves(i)%x(j - 1),Curves(i)%y(j - 1),2,d1,d0,0)&
                        - pg_get_fun_xy(Curves(i)%x(j - 2),Curves(i)%y(j - 2),2,d1,d0,0))/(Curves(i)%s(j) - Curves(i)%s(j - 2))
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
subroutine fcn_s_t(n, t, y, yprime) !интегрирование по длине дуги s с поиском времени
    use mod
    integer(4) :: n
    real(8) t, y(n), yprime(n), u_x, u_y, pg_get_fun_xy, V
    if (y(3) < eps) then
        y(3)           = eps
    end if
    
    if (y(4) < eps) then 
        y(4)           = eps
    end if 
    
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
    if (y(3) < eps) then
        y(3)           = eps
    end if
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
    if (y(3) < eps) then
        y(3)           = eps
    end if
    V           = abs(y(3))
    u_x         = pg_get_fun_xy(y(1),y(2)+eps,2,d0,d1,0)
    yprime(1)   = y(3)/V
    yprime(2)   = d0
    yprime(3)   = (u_x-y(3))/(st*V)
    yprime(4)   = d0
    yprime(5)   = 1/V
    end subroutine fcn_s_t_bottom
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
        tol         = 0.0001d0
        s           = d0
        y           = d0 ! y = [J11, J12, J21, J22, w11, w12, w21, w22]
        y(1)        = d1
        y(4)        = d1
        param       = d0
        !param(4)    = current_Curve%n
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
        ! current_Curve%J_11(j)   = y(1)
        ! current_Curve%J_12(j)   = y(2)
        ! current_Curve%J_21(j)   = y(3)
        ! current_Curve%J_22(j)   = y(4)
        ! current_Curve%J(j)      = abs(y(1)*y(4) - y(2)*y(3))
        call divprk(3, n, fcn_t_Jacobian, current_Curve%t(j - 1), current_Curve%t(j), tol, param, y)
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
    if (V < eps) then
        V           = eps
    end if
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
    integer(4) index_t
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
            if (Curves(i)%c(k) > 1) then
                Curves(i)%c(k) = 2d0 - 1d0/Curves(i)%c(k)
            end if 
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
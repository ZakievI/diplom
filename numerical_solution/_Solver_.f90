subroutine solve
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
subroutine build_bound()
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
subroutine fcn(n, t, y, yprime)
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
subroutine fcn_s(n, t, y, yprime)
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
function function_impact_testing(n, y, dlt)
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
subroutine impact_test(n, y , y_out, num, dlt)
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
subroutine coordinate_first_particle(n, y, arr_bound, num, dlt, size_arr)
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
subroutine build_curve()
    use mod
    type :: Curve
        integer(4) :: n
        real(8) :: x(N_arr), y(N_arr), t(N_arr)
    end type
    integer(4) :: num_particle, i, n1, ido
    integer(4), parameter :: n      = 5
    !integer(4), parameter :: n      = 4
    integer(4), parameter :: mxparm = 50
    real(8) :: param(mxparm), d_s, s, y(n), dlt, tol, pg_get_fun_xy
    type(Curve), allocatable :: Curves(:)
    external fcn_s_t
    num_particle    = 5
    tol             = 0.001d0
    d_s             = 0.05d0
    s               = d0
    st              = 1d0
    ido             = 1
    dlt             = d0
    allocate(Curves(num_particle))
    n1 = 1
    Curves(1)%y(n1) = eps
    Curves(1)%x(n1) = -L1/2
    Curves(1)%t(n1) = d0
    y(1)            = -L1/2
    y(2)            = eps 
    y(3)            = pg_get_fun_xy(y(1),y(2),2,d0,d1,0)
    ! y(4)            = -pg_get_fun_xy(y(1),y(2),2,d1,d0,0)
    y(4)            = d0
    y(5)            = d0
    param           = d0
    param(4)        = N_arr
    param(10)       = 1.0d0
    do while((sqrt(y(1)**2+y(2)**2)>d1+dlt).and.(-L1/2<=y(1)).and.(y(1)<=L1/2).and.(n1<=N_arr))
        n1 = n1 + 1
        call divprk(ido, n, fcn_s_t, s, s+d_s, tol, param, y)
        Curves(1)%x(n1) = y(1)
        Curves(1)%y(n1) = d0
        Curves(1)%t(n1) = y(5)
    end do
    Curves(1)%n = n1
    !call divprk(3, n, fcn_s_t, s, s+d_s, tol, param, y)
    do i = 2, num_particle-1
        n1 = 1
        Curves(i)%y(n1) = H1*(i - d1)/(num_particle - d1)
        Curves(i)%x(n1) = -L1/2
        Curves(i)%t(n1) = d0
        y(1)            = -L1/2
        y(2)            = H1*(i - d1)/(num_particle - d1) 
        y(3)            = pg_get_fun_xy(y(1),y(2),2,d0,d1,0)
        y(4)            = -pg_get_fun_xy(y(1),y(2),2,d1,d0,0)
        y(5)            = Curves(i)%t(n1)
        param           = d0
        param(4)        = N_arr
        param(10)       = 1.0d0
        do while((sqrt(y(1)**2+y(2)**2)>d1+dlt).and.(-L1/2<=y(1)).and.(y(1)<=L1/2).and.(y(2)<=H1).and.(d0<=y(2)).and.(n1<=N_arr))
            n1 = n1 + 1
            call divprk(ido, n, fcn_s_t, s, s+d_s, tol, param, y)
            Curves(i)%x(n1) = y(1)
            Curves(i)%y(n1) = y(2)
            Curves(i)%t(n1) = y(5)
        end do
        Curves(i)%n = n1
    end do 
    n1 = 1
    Curves(num_particle)%y(n1) = H1
    Curves(num_particle)%x(n1) = -L1/2
    Curves(num_particle)%t(n1) = d0
    y(1)            = -L1/2
    y(2)            = H1 - eps
    y(3)            = pg_get_fun_xy(y(1),y(2),2,d0,d1,0)
    ! y(4)            = -pg_get_fun_xy(y(1),y(2),2,d1,d0,0)
    y(4)            = d0
    y(5)            = d0
    param           = d0
    param(4)        = N_arr
    param(10)       = 1.0d0
    do while((sqrt(y(1)**2+y(2)**2)>d1+dlt).and.(-L1/2<=y(1)).and.(y(1)<=L1/2).and.(n1<=N_arr))
        n1 = n1 + 1
        call divprk(ido, n, fcn_s_t, s, s+d_s, tol, param, y)
        Curves(num_particle)%x(n1) = y(1)
        Curves(num_particle)%y(n1) = H1
        Curves(num_particle)%t(n1) = y(5)
    end do
    Curves(num_particle)%n = n1
    call divprk(3, n, fcn_s_t, s, s+d_s, tol, param, y)
end subroutine build_curve
subroutine fcn_s_t(n, t, y, yprime)
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
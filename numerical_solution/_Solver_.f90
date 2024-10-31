subroutine solve
    use mod
    implicit none
    integer(4), parameter::n = 4, number = 100
    integer(4)::ido, k1, i, num, f_t(3)
    integer(4) function_impact_testing
    real(8), allocatable :: y_out(:,:) 
    real(8):: y(n), dlt, stocs, T_num, y_t(3)
    real(8):: pg_get_fun_xy, y_0
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
    dlt             = 0.1d0
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
        !print *, y_0, " ", T_num*H1/number
    end do 
    deallocate(y_out)
end subroutine solve
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
    s           = d0                !Начальные условия при интегрировании по дуге
    tol         = 0.001d0           !Допустимая ошибка
    param       = d0                !Действуем по умолчанию
    param(4)    = num               !Максимальное кол-во итераций
    param(10)   = 1.0d0             !Контролируем абсолютную ошибику
    d_s         = 0.05d0            !Шаг по дуге
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
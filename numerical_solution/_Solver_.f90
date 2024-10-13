subroutine solve
    use mod
    implicit none
    integer(4), parameter::mxparm = 100, n = 4, number = 10
    integer(4)::ido, nout, k1, i, j, num
    real(8), allocatable :: y_out(:,:) 
    real(8)::param(mxparm), t, tend, tol, y(n), istep , s, dlt
    real(8):: y_old(n), sred,dt, pg_get_fun_xy
    external fcn, fcn_s
    open (1, FILE='traektorie.dat')
    write(1,*) 'TITLE = "traektorie"'
    write(1,*) 'VARIABLES = "X", "Y"'
    !write(1,"('ZONE T=""area"", I=', i0, ', J=', i0, ', F=POINT')") nr+1,ng+1
    print "(4x, 'istep' , 5x, 'time', 9x, 'x', 11x, 'y')"
    num = 5000
    dlt         = 0.001d0
    allocate(y_out(num,2))

        t           = d0                !Начальные условия при интегрировании по времени
        y(1)        = -L1/2
        y(2)        = H1*(i)/(2*(number+1))
        y(3)        = pg_get_fun_xy(y(1),y(2),2,d0,d1,0)
        y(4)        = -pg_get_fun_xy(y(1),y(2),2,d1,d0,0)
        s           = d0                !Начальные условия при интегрировании по дуге
        tol         = 0.001d0         !Допустимая ошибка
        param       = d0                !Действуем по умолчанию
        param(4)    = num
        param(10)   = 1.0d0             !Контролируем абсолютную ошибику
        !write(2,*) "(4x, 'istep' , 5x, 'time', 9x, 'x', 11x, 'y', 13x, 'v_x', 15x, 'v_y')"
        ido         = 1
        dt          = 0.1d0
        ds          = 0.1d0
        k1          = 1
        y_out(k1,1) =y(1)
        y_out(k1,2) =y(2)
        k1 = 2
        do while ((sqrt(y(1)**2+y(2)**2)>d1+dlt).and.(-L1/2<=y(1)).and.(y(1)<L1/2).and.(y(2)<H1).and.(d0<y(2)).and.(k1<=num))
            !call divprk(ido, n, fcn, t, t+dt, tol, param, y)
            call divprk(ido, n, fcn_s, s, s+ds, tol, param, y)
            !print '(i6, 6f12.3)', k1, t, y
            y_out(k1,1)=y(1)
            y_out(k1,2)=y(2)
            k1 = k1 + 1
        end do
        print '(i0)', i
        write(1,"('ZONE T=""area',i0,'"", I=', i0, ', F=POINT')") i, k1-1
        do j = 1, k1-1
            write(1,"(F9.5, ' ', F9.5)") y_out(j,1),y_out(j,2)
            end do 
        !call divprk(3, n, fcn, t, t+dt, tol, param, y)
        call divprk(3, n, fcn_s, s, s+ds, tol, param, y)
    end do
    deallocate(y_out)
end subroutine solve
subroutine fcn(n, t, y, yprime)
    use mod
    integer(4) :: n
    real(8) t, y(n), yprime(n), u_x, u_y, pg_get_fun_xy
    ! st=1, 10, 100
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
    ! st=1, 10, 100, 0.02d0
    real(8), parameter :: g = 9.81d0, v0 = 10.0d0
    V           = dsqrt(y(3)**2+y(4)**2)
    u_x         = pg_get_fun_xy(y(1),y(2),2,d0,d1,0)
    u_y         = -pg_get_fun_xy(y(1),y(2),2,d1,d0,0)
    yprime(1)   = y(3)/V
    yprime(2)   = y(4)/V
    yprime(3)   = (u_x-y(3))/(st*V)
    yprime(4)   = (u_y-y(4))/(st*V)
end subroutine fcn_s
subroutine testing()
    use mod
    integer(4), parameter::mxparm = 100, n = 4 
    integer(4)::ido, nout
    real(8)::param(mxparm), t, tend, tol, y(n), istep, alfa, y_analit_max, x_analit_max, y_max, x_max, x_zero
    real(8), parameter :: g = 9.81, v0 = 10.0
    real(8):: y_old(n), sred,dt
    external fcn1
    alfa        = pi/4
    t           = 0.0           !Начальные условия
    y(1)        = 0.0
    y(2)        = 0.0
    y(3)        = v0*dcos(alfa)
    y(4)        = v0*dsin(alfa)
    tol         = 0.00001       !Допустимая ошибка
    param       = 0.0           !Действуем по умолчанию
    param(10)   = 1.0           !Контролируем абсолютную ошибику
    y_analit_max= v0**2*dsin(alfa)*dsin(alfa)/(2*g)
    x_analit_max= 2*v0**2*dcos(alfa)*dsin(alfa)/(g)
    !Вывод заголовока таблицы результатов
    open (2, FILE='test.dat')
    print "(4x, 'istep' , 5x, 'time', 9x, 'x', 11x, 'y')"
    write(2,*) "(4x, 'istep' , 5x, 'time', 9x, 'x', 11x, 'y', 13x, 'v_x', 15x, 'v_y')"
    ido = 1
    dt=0.1d0
    do while (ido.ne.3)
        call divprk(ido, n, fcn, t, t+dt, tol, param, y)
        print '(i6, 5f12.3)', istep, t, y
        WRITE(2,*) istep, t, y 
        if ((y_old(2)>0).and.(y(2)<0)) then
            x_zero = sred(y_old(2),y(2),y_old(1),y(1),d0)
            print *, x_zero, v0**2*dsin(2*alfa)/g
            ido =3
        end if 
        y_old = y
        !if (istep == 2) ido = 3 !Освобождаем память
    end do    
    call divprk(ido, n, fcn1, t, t+dt, tol, param, y)
end subroutine testing

subroutine fcn1(n, t, y, yprime)
    use mod
    integer(4) :: n
    real(8) t, y(n), yprime(n)
    real(8), parameter :: g = 9.81d0, v0 = 10.0d0
    yprime(1) = y(3) 
    yprime(2) = y(4)
    yprime(3) = d0
    yprime(4) = -g
end subroutine fcn1
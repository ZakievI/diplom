
subroutine testing()
    integer(4), parameter::mxparm = 50, n = 2 
    integer(4)::ido, nout
    real(4)::param(mxparm), t, tend, tol, y(n), istep
    external fcn
    t           = 0.0           !Начальные условия
    y(1)        = 0.0
    y(2)        = 0.0
    tol         = 0.00001       !Допустимая ошибка
    param       = 0.0           !Действуем по усолчанию
    param(10)   = 1.0           !Контролируем абсолютную ошибику
    !Вывод заголовока таблицы результатов
    print "(4x, 'istep' , 5x, 'time', 9x, 'x', 11x, 'y')"
    ido = 1
    istep = 0
    do while (istep<2)
        istep = istep + 0.1
        tend = istep
        call ivprk(ido, n, fcn, t, tend, tol, param, y)
        print '(i6, 3f12.3)', istep, t, y
        if (istep == 2) ido = 3 !Освобождаем память
    end do    
end subroutine testing

subroutine fcn(n, t, y, yprime)
    integer(4) :: n
    real(4) :: t, y(n), yprime(n)
    real(8) :: g = 9.81, v0 = 10
    yprime(1) = v0/sqrt(2.0)
    yprime(2) = v0/sqrt(2.0) - g*t
end subroutine fcn
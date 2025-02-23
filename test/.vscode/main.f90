! program main
!     real(4) x1,x2,x3,y1,y2,y3,x4,y4,h
!     print *,'enter x1 and y1'
!     read*, x1,y1
!     print *,'enter x2 and y2'
!     read*, x2,y2
!     print *,'enter x3 and y3'
!     read*, x3,y3
!     print *,'A(',x1,' ',y1,')','B(',x2,' ',y2,')','C(',x3,' ',y3,')'
    
!     if ((y3==y2 .and. x3==x2) .or. (x1==x2 .and.y1==y2) .or.(x1==x3 .and. y1==y3)) then
!         print*,'EROOR!!!, this is LINE!!!!!!'
!     else if(x1==x2 .and. y1==y2 .and. y2==y3 .and.x2==x3) then
!         h=0
!         x4=x1
!         y4=y1
!     else if(y3==y2) then
!         h=ABS(y1-y2)
!         x4=x1
!         y4=y2
!     else if(y1==y3) then
!         h=ABS(y1-y2)
!         x4=x1
!         y4=y2
!     else if(x3==x2) then
!         h=ABS(x1-x2)
!         x4=x2
!         y4=y1
!     else 
!         h=ABS((y3-y2)*(x1-x2)-(x3-x2)*(y1-y2))/(SQRT((y3-y2)**2+(x3-x2)**2))
!         y4=(x1-x2+((y3-y2)/(x3-x2))*y1+((x3-x2)/(y3-y2))*y2)/((((y3-y2)/(x3-x2))+((x3-x2)/(y3-y2))))
!         x4=((x3-x2)*(y4-y2))/(y3-y2)+x2
!     end if
!     print *, 'h=',h
!     print *, 'D=(',x4,' ',y4,')'
! end program main
! program main
!     use module1
!     use module2
!     implicit none

!     call hello_from_module1()
!     call hello_from_module2()
! end program main
program main
    use imsl
    implicit none
    integer :: n = 10
    real(8) :: a(n), b(n), c(n)
    integer :: i

    ! Initialize arrays
    do i = 1, n
        a(i) = dble(i)
        b(i) = dble(i)
    end do

    ! Call IMSL function
    call vdAdd(n, a, b, c)

    ! Print result
    do i = 1, n
        print *, c(i)
    end do
end program main

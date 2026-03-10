!итерационный процесс для решения уравнения Навье-Стокса
  
subroutine pg_navie_stoks_iterate
!dec$ attributes dllexport:: pg_navie_stoks_iterate
!итерационный процесс для решения уравнения Навье-Стокса
use pgmod
external closing_navie_stoks_empty
real(8) init_val(1)
call navie_stoks_iterate(closing_navie_stoks_empty,init_val,init_val,.false.)
end

subroutine pg_navie_stoks_iterate1(f_closing)
!dec$ attributes dllexport:: pg_navie_stoks_iterate1
!итерационный процесс для решения уравнения Навье-Стокса
!c процедурой замыкания
use pgmod
external f_closing
real(8) init_val(1)
call navie_stoks_iterate(f_closing,init_val,init_val,.false.)
end

subroutine pg_navie_stoks_iterate2(init_val1,init_val2)
!dec$ attributes dllexport:: pg_navie_stoks_iterate2
!итерационный процесс для решения уравнения Навье-Стокса
!с начальными приближениями
use pgmod
external closing_navie_stoks_empty
real(8) init_val1(gsarea%a%ntr), init_val2(gsarea%a%ntr)
call navie_stoks_iterate(closing_navie_stoks_empty,init_val1,init_val2,.true.)
end

subroutine pg_navie_stoks_iterate3(f_closing,init_val1,init_val2)
!dec$ attributes dllexport:: pg_navie_stoks_iterate3
!итерационный процесс для решения уравнения Навье-Стокса
!c процедурой замыкания и начальными приближениями
use pgmod
external f_closing
real(8) init_val1(gsarea%a%ntr),init_val2(gsarea%a%ntr)
call navie_stoks_iterate(f_closing,init_val1,init_val2,.true.)
end

subroutine navie_stoks_iterate(f_closing,init_val1,init_val2,have_init_val)
use pgmod
external f_closing
integer(4) i,imax
real(8) err,Re
logical Re_iter
real(8) time,rtc,time1
real(8) init_val1(gsarea%a%ntr), init_val2(gsarea%a%ntr)
logical have_init_val
character(50) str
Re_iter=.false.
imax=1000
Re=gsarea%const%Re
!if (Re_iter.and.gsarea%const%Re>d1) Re=d1
if (Re_iter.and.gsarea%const%Re>10.0d0) Re=10.0d0
time=rtc()
select case (gsarea%type_eq(1))
case (21)
  call pg_init_need_cash_for_func(2,1,2,.true.)
  call pg_init_need_cash_for_func(2,2,2,.false.)
  call pg_init_need_cash_for_func(4,1,2,.false.)
  call pg_init_need_cash_for_func(4,2,2,.false.)
case (24)
  call pg_init_need_cash_for_func(4,1,2,.true.)
  call pg_init_need_cash_for_func(4,2,2,.false.)
case (25)
  call pg_init_need_cash_for_func(2,1,2,.true.)
  call pg_init_need_cash_for_func(2,2,2,.false.)
end select
call pg_initCash_area
!call init_ba_cash(1,1)
!call init_aa_cash(1)
do i=1,imax
  call init_meshval_ns(Re,i==1,err,init_val1,init_val2,have_init_val)
  call gs_print('-----------------------')
  write(str,"('i=',i0,'  Re=',F10.2,'  err=',E13.5)") i,Re,err
  call gs_print(trim(str))
  call gs_print('-----------------------')
  if (i>1) then
    if (err<1.0d-6) exit
    if (Re_iter.and.Re<gsarea%const%Re.and.err<1.0d-1) then
      !Re=Re*2.0d0
      Re=Re+2.0d0
      if (Re>gsarea%const%Re) Re=gsarea%const%Re
    endif
  endif
  call pg_get_matrix
  call f_closing
  call pg_solve
enddo
time1=rtc()
call gs_print('Navie-Stoks iter time=  '//ftoa(time1-time))
end

subroutine init_meshval_ns(Re,is_init,err,init_val1,init_val2,have_init_val)
use pgmod
real(8) pg_get_fun_1_ar,pg_get_fun2_1_ar,Re
real(8) val(gsarea%a%ntr)
integer(4) i,k,nk,j,type_eq
logical is_init
real(8) err,erre,errt
integer(4) ierr
real(8) init_val1(gsarea%a%ntr), init_val2(gsarea%a%ntr)
logical have_init_val
type_eq=gsarea%type_eq(1)
err=d0
nk=1
if (type_eq==24.or.type_eq==25) nk=2
do k=1,nk
  select case (type_eq)
  case (21)
    j=1
  case (24)
    j=3+k
  case (25)
    j=5+k
  end select
  if (Re==d0) then
    val=d0
  elseif (is_init) then
    if (have_init_val) then
      if (k==1) then
        val=init_val1
      else
        val=init_val2
      endif
    else
      val=d0
    endif
  else
    select case (type_eq)
    case (21)
      do i=1,gsarea%a%ntr
        val(i)=Re*(pg_get_fun_1_ar(i,2)*pg_get_fun2_1_ar(i,1)-pg_get_fun_1_ar(i,1)*pg_get_fun2_1_ar(i,2))
      enddo
    case (24)
      if (k==1) then
        do i=1,gsarea%a%ntr
          val(i)=-Re*pg_get_fun2_1_ar(i,2)
        enddo
      else
        do i=1,gsarea%a%ntr
          val(i)=Re*pg_get_fun2_1_ar(i,1)
        enddo
      endif
    case (25)
      if (k==1) then
        do i=1,gsarea%a%ntr
          val(i)=Re*pg_get_fun_1_ar(i,2)
        enddo
      else
        do i=1,gsarea%a%ntr
          val(i)=-Re*pg_get_fun_1_ar(i,1)
        enddo
      endif
    end select
  endif
  call pg_init_area_gu_test(val,j,errt,ierr,erre)
  if (err<errt) err=errt
enddo
end

subroutine closing_navie_stoks_empty
end
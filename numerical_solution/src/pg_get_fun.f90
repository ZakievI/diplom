function pg_get_fun(x,y)
!dec$ attributes dllexport:: pg_get_fun
!вычисление главной функции в текущей области
use pgmod
real(8) x,y
real(8) pg_get_fun,get_fun_
pg_get_fun=get_fun_(0,x,y,0,0)
end

function pg_get_fun_m(j,knd)
!dec$ attributes dllexport:: pg_get_fun_m
!вычисление главной функции в текущей области
!в контрольной точке (центр панели)
use pgmod
integer(4) j    !номер панели
integer(4) knd  !номер границы
real(8) pg_get_fun_m,get_fun_
pg_get_fun_m=get_fun_(1,d0,d0,knd,j)
end

function pg_get_fun_ar(j)
!dec$ attributes dllexport:: pg_get_fun_ar
!вычисление главной функции в текущей области
!в контрольной точке (центр ячейки)
use pgmod
integer(4) j
real(8) pg_get_fun_ar,get_fun_
pg_get_fun_ar=get_fun_(2,d0,d0,0,j)
end

function get_fun_(mode,x,y,i_0,i_1)
!вычисление главной функции в текущей области
use pgmod
integer(4) type_eq
integer(4) mode,i_0,i_1
real(8) x,y
real(8) get_fun_,get_fun,get_fun2,get_funp,get_funh,get_funh2,get_funbri,get_funoss,get_fun2p
call bind_AreaConst(gsarea%i)
type_eq=gsarea%type_eq(1)
select case (type_eq)
case (1)
  get_fun_=get_fun(mode,x,y,i_0,i_1,0,gsarea%i)
case (3)
  get_fun_=get_fun2(mode,x,y,i_0,i_1,gsarea%i)
case (4:6,8,11,16:19,26)
  get_fun_=get_funp(mode,x,y,i_0,i_1,0,gsarea%i)
case (9)
  get_fun_=get_funh(mode,x,y,i_0,i_1,0,gsarea%i)
case (10)
  get_fun_=get_funh2(mode,x,y,i_0,i_1,0,gsarea%i)
case (15)
  get_fun_=get_funbri(mode,x,y,i_0,i_1,gsarea%i)
case (20)
  get_fun_=get_funoss(mode,x,y,i_0,i_1,gsarea%i)
case (21,22,24,25)
  get_fun_=get_fun2p(mode,x,y,i_0,i_1,gsarea%i)
end select
end

function pg_get_fun_1(x,y,np)
!dec$ attributes dllexport:: pg_get_fun_1
!вычисление производной главной функции в текущей области
!np - тип производной    1 - dpsi/dx, 2 - dpsi/dy
use pgmod
integer(4) np
real(8) x,y
real(8) pg_get_fun_1,get_fun_1_
pg_get_fun_1=get_fun_1_(0,x,y,0,0,np)
end

function pg_get_fun_1_m(j,knd,np)
!dec$ attributes dllexport:: pg_get_fun_1_m
!вычисление производной главной функции в текущей области
!в контрольной точке (центр панели)
!np - тип производной    1 - dpsi/dx, 2 - dpsi/dy
use pgmod
integer(4) np
integer(4) j    !номер панели
integer(4) knd  !номер границы
real(8) pg_get_fun_1_m,get_fun_1_
pg_get_fun_1_m=get_fun_1_(1,d0,d0,knd,j,np)
end

function pg_get_fun_1_ar(j,np)
!dec$ attributes dllexport:: pg_get_fun_1_ar
!вычисление производной главной функции в текущей области
!в контрольной точке (центр ячейки)
!np - тип производной    1 - dpsi/dx, 2 - dpsi/dy
use pgmod
integer(4) np
integer(4) j
real(8) pg_get_fun_1_ar,get_fun_1_
pg_get_fun_1_ar=get_fun_1_(2,d0,d0,0,j,np)
end

function get_fun_1_(mode,x,y,i_0,i_1,np)
!вычисление производной главной функции в текущей области
!np - тип производной    1 - dpsi/dx, 2 - dpsi/dy
use pgmod
integer(4) mode,i_0,i_1
real(8) x,y
integer(4) type_eq,np
real(8) get_fun_1_,get_fun_1,get_fun2_1,get_funp_1,get_funbri_1,get_funoss_1,get_fun2p_1,get_funh2_1 !,get_funh_1
call bind_AreaConst(gsarea%i)
type_eq=gsarea%type_eq(1)
select case (type_eq)
case (1)
  get_fun_1_=get_fun_1(mode,x,y,i_0,i_1,0,np,gsarea%i)
case (3)
  get_fun_1_=get_fun2_1(mode,x,y,i_0,i_1,np,gsarea%i)
case (4:6,8,11,16:19,26)
  get_fun_1_=get_funp_1(mode,x,y,i_0,i_1,0,np,gsarea%i)
!case (9)
!  get_fun_1_=get_funh_1(mode,x,y,i_0,i_1,0,np,gsarea%i)
case (10)
  get_fun_1_=get_funh2_1(mode,x,y,i_0,i_1,0,np,gsarea%i)
case (15)
  get_fun_1_=get_funbri_1(mode,x,y,i_0,i_1,np,gsarea%i)
case (20)
  get_fun_1_=get_funoss_1(mode,x,y,i_0,i_1,np,gsarea%i)
case (21,22,24,25)
  get_fun_1_=get_fun2p_1(mode,x,y,i_0,i_1,np,gsarea%i)
case default
  call gs_print_stop("Error get_fun_1_")
end select
end

subroutine pg_get_fun_1_dual(x,y,val,val2)
!dec$ attributes dllexport:: pg_get_fun_1_dual
!вычисление градиента главной функции в текущей области с использованием одновременного вычисления производных
use pgmod
real(8) x,y,val,val2
integer(4) type_eq
call bind_AreaConst(gsarea%i)
type_eq=gsarea%type_eq(1)
select case (type_eq)
case (1)
  call get_fun_1_dual(x,y,0,gsarea%i,val,val2)
case (3)
  call get_fun2_1_dual(x,y,gsarea%i,val,val2)
case (4:6,8,11,16:19,26)
  call get_funp_1_dual(x,y,0,gsarea%i,val,val2)
!case (9)
!  call get_funh_1_dual(x,y,0,gsarea%i,val,val2)
case (10)
  call get_funh2_1_dual(x,y,0,gsarea%i,val,val2)
case (15)
  call get_funbri_1_dual(x,y,gsarea%i,val,val2)
!case (20)
!  get_fun_1_=get_funoss_1(mode,x,y,i_0,i_1,np,gsarea%i)
case (21,22,24,25)
  call get_fun2p_1_dual(x,y,gsarea%i,val,val2)
case default
  call gs_print_stop("Error pg_get_fun_1_dual")
end select
end

function pg_get_fun2(x,y)
!dec$ attributes dllexport:: pg_get_fun2
!вычисление второй функции в текущей области
use pgmod
real(8) x,y
real(8) pg_get_fun2,get_fun2_
pg_get_fun2=get_fun2_(0,x,y,0,0)
end

function pg_get_fun2_m(j,knd)
!dec$ attributes dllexport:: pg_get_fun2_m
!вычисление второй функции в текущей области
!в контрольной точке (центр панели)
use pgmod
integer(4) j    !номер панели
integer(4) knd  !номер границы
real(8) pg_get_fun2_m,get_fun2_
pg_get_fun2_m=get_fun2_(1,d0,d0,knd,j)
end

function pg_get_fun2_ar(j)
!dec$ attributes dllexport:: pg_get_fun2_ar
!вычисление второй функции в текущей области
!в контрольной точке (центр ячейки)
use pgmod
integer(4) j
real(8) pg_get_fun2_ar,get_fun2_
pg_get_fun2_ar=get_fun2_(2,d0,d0,0,j)
end

function get_fun2_(mode,x,y,i_0,i_1)
!вычисление второй функции в текущей области
use pgmod
integer(4) type_eq
integer(4) mode,i_0,i_1
real(8) x,y
real(8) get_fun2_,get_fun,get_funp,get_funh,get_funh2
call bind_AreaConst(gsarea%i)
get_fun2_=d0
if (gsarea%nu==1) return
type_eq=gsarea%type_eq(2)
select case (type_eq)
case (2)
  get_fun2_=get_fun(mode,x,y,i_0,i_1,2,gsarea%i)
case (7,14,23)
  get_fun2_=get_funp(mode,x,y,i_0,i_1,2,gsarea%i)
case (12)
  get_fun2_=get_funh(mode,x,y,i_0,i_1,2,gsarea%i)
case (13)
  get_fun2_=get_funh2(mode,x,y,i_0,i_1,2,gsarea%i)
end select
end

function pg_get_fun2_1(x,y,np)
!dec$ attributes dllexport:: pg_get_fun2_1
!вычисление производной второй функции в текущей области
!np - тип производной    1 - dpsi/dx, 2 - dpsi/dy
use pgmod
integer(4) np
real(8) x,y
real(8) pg_get_fun2_1,get_fun2_1_
pg_get_fun2_1=get_fun2_1_(0,x,y,0,0,np)
end

function pg_get_fun2_1_m(j,knd,np)
!dec$ attributes dllexport:: pg_get_fun2_1_m
!вычисление производной второй функции в текущей области
!в контрольной точке (центр панели)
!np - тип производной    1 - dpsi/dx, 2 - dpsi/dy
use pgmod
integer(4) np
integer(4) j    !номер панели
integer(4) knd  !номер границы
real(8) pg_get_fun2_1_m,get_fun2_1_
pg_get_fun2_1_m=get_fun2_1_(1,d0,d0,knd,j,np)
end

function pg_get_fun2_1_ar(j,np)
!dec$ attributes dllexport:: pg_get_fun2_1_ar
!вычисление производной второй функции в текущей области
!в контрольной точке (центр ячейки)
!np - тип производной    1 - dpsi/dx, 2 - dpsi/dy
use pgmod
integer(4) np
integer(4) j
real(8) pg_get_fun2_1_ar,get_fun2_1_
pg_get_fun2_1_ar=get_fun2_1_(2,d0,d0,0,j,np)
end

function get_fun2_1_(mode,x,y,i_0,i_1,np)
!вычисление производной второй функции в текущей области
!np - тип производной    1 - dpsi/dx, 2 - dpsi/dy
use pgmod
integer(4) type_eq,np
integer(4) mode,i_0,i_1
real(8) x,y
real(8) get_fun2_1_,get_fun_1,get_funp_1,get_funh2_1 !,get_funh_1
call bind_AreaConst(gsarea%i)
get_fun2_1_=d0
if (gsarea%nu==1) return
type_eq=gsarea%type_eq(2)
select case (type_eq)
case (2)
  get_fun2_1_=get_fun_1(mode,x,y,i_0,i_1,2,np,gsarea%i)
case (7,14,23)
  get_fun2_1_=get_funp_1(mode,x,y,i_0,i_1,2,np,gsarea%i)
!case (12)
!  get_fun2_1_=get_funh_1(mode,x,y,i_0,i_1,2,np,gsarea%i)
case (13)
  get_fun2_1_=get_funh2_1(mode,x,y,i_0,i_1,2,np,gsarea%i)
case default
  call gs_print_stop("Error get_fun2_1_")
end select
end

subroutine pg_get_fun2_1_dual(x,y,val,val2)
!dec$ attributes dllexport:: pg_get_fun2_1_dual
!вычисление градиента второй функции в текущей области с использованием одновременного вычисления производных
use pgmod
integer(4) type_eq
real(8) x,y
real(8) val,val2
call bind_AreaConst(gsarea%i)
val=d0
val2=d0
if (gsarea%nu==1) return
type_eq=gsarea%type_eq(2)
select case (type_eq)
case (2)
  call get_fun_1_dual(x,y,2,gsarea%i,val,val2)
case (7,14,23)
  call get_funp_1_dual(x,y,2,gsarea%i,val,val2)
!case (12)
!  call get_funh_1_dual(x,y,2,gsarea%i,val,val2)
case (13)
  call get_funh2_1_dual(x,y,2,gsarea%i,val,val2)
case default
  call gs_print_stop("Error pg_get_fun2_1_dual")
end select
end

function get_fun(mode,x0,y0,i_0,i_1,k,ia)
!решение уравнения лапласа в области
use pgmod
integer(4) mode !0 - (x0,y0), 1 - (i_0-bnd, i_1-panel), 2 - (i_1-ntr)
integer(4) i_0,i_1
real(8) x0,y0,s,val,ff
real(8) calcint,get_fun,cp_calcint,calcintm,calcint_ar
!real(8) calcint2
integer(4) j,k,j0,ia,i,k1
type(TArea), pointer :: a
type(TBound), pointer :: b
type(TBoundline2), pointer :: bl2
type(TBounLineFuncApprox), pointer :: ga
integer(4) nf(3),nf2(3),nf1(2)
real(8) kk(2),fff(3)
data nf /1,22,23/
data nf2 /0,20,21/
data nf1 /1,0/
data kk /d1,-d1/
complex(8) ztc
a=>gs%a(ia)
s=d0
!$omp parallel do if (gs_use_parallel_get_fun==1) reduction(+ : s) private (j0,j,i,b,bl2,ztc,ga,val,ff,k1,fff) shared (a,x0,y0,k,ia,nf,nf2,gs,mode,i_0,i_1,kk,nf1)
do j0=1,a%nb
  b=>a%bnd(j0)
  if (b%skip_getfun) cycle
  select case (b%boundLineType)
  case (1)
    !$omp parallel do if (gs_use_parallel_get_fun==1) reduction(+ : s) private (j,val,ff,k1) shared (b,k,x0,y0,ia,j0,gs,mode,i_0,i_1,kk,nf1)
    do j=1,b%npanel
      do k1=1,2
        val=b%psiom(j,k+k1)
        if (val.ne.d0) then
          select case (mode)
          case (1)
            ff=calcintm(j,i_1,nf1(k1),ia,j0,i_0)
          case (2)
            ff=calcint_ar(j,i_1,nf1(k1),ia,j0)
          case default
            ff=calcint(j,x0,y0,nf1(k1),ia,j0)
          end select
          s=s+kk(k1)*val*ff
        endif
	    enddo
    enddo
    !omp end parallel do nowait
  case (2)
    if (mode.ne.0) call gs_print_stop("Error get_fun")
    do j=1,b%nline
	    bl2=>b%line2(j)
	    ztc=bl2%zc-dcmplx(x0,y0)
      if (gs_use_dual_integral_solve) then
        ga=>bl2%ga(k+1)
        val=ga%valkk(1)
        if (val.ne.d0) s=s+val*cp_calcint(ztc,bl2%r,1,bl2%g0,bl2%dir,0)
        call cp_calcint_nf2223_2(ztc,bl2%r,0,fff(2),fff(3))
	      do i=2,ga%n
	        val=ga%valkk(i)
		      if (val.ne.d0) s=s+val*fff(i)
	      enddo	
	      ga=>bl2%ga(k+2)
        val=ga%valkk(1)
        if (val.ne.d0) s=s-val*cp_calcint(ztc,bl2%r,0,bl2%g0,bl2%dir,0)
        call cp_calcint_nf2021_2(ztc,bl2%r,0,fff(2),fff(3))
	      do i=2,ga%n
	        val=ga%valkk(i)
		      if (val.ne.d0) s=s-val*fff(i)
	      enddo	
      else
	      ga=>bl2%ga(k+1)
	      do i=1,ga%n
	        val=ga%valkk(i)
		      if (val.ne.d0) s=s+val*cp_calcint(ztc,bl2%r,nf(i),bl2%g0,bl2%dir,0)
	      enddo	
	      ga=>bl2%ga(k+2)
	      do i=1,ga%n
	        val=ga%valkk(i)
		      if (val.ne.d0) s=s-val*cp_calcint(ztc,bl2%r,nf2(i),bl2%g0,bl2%dir,0)
        enddo	
      endif
	  enddo
  end select
enddo
!omp end parallel do nowait
get_fun=s/pi2
end

function get_fun_1(mode,x0,y0,i_0,i_1,k,np,ia)
!вычисление производной функции - решения уравнения лапласа в области
!np - тип производной    1 - dpsi/dx, 2 - dpsi/dy
use pgmod
integer(4) mode !0 - (x0,y0), 1 - (i_0-bnd, i_1-panel), 2 - (i_0-ntr)
integer(4) i_0,i_1,np1
real(8) x0,y0,s,val,ff
real(8) calcint,get_fun_1,cp_calcint,calcintm,calcint_ar
integer(4) j,k,np,j0,ia,i,k1
type(TArea), pointer :: a
type(TBound), pointer :: b
type(TBoundline2), pointer :: bl2
type(TBounLineFuncApprox), pointer :: ga
integer(4) nf(3),nf2(3),nf1(2),nf1_(2)
real(8) kk(2)
data nf /6,32,34/
data nf2 /4,28,30/
data nf1 /6,4/
data kk /d1,-d1/
complex(8) ztc,z0
np1=np-1
nf1_=nf1+np1
a=>gs%a(ia)
s=d0
z0=dcmplx(x0,y0)
!$omp parallel do if (gs_use_parallel_get_fun==1) reduction(+ : s) private (j0,j,b,val, i,bl2,ga,ztc,k1,ff) shared (a,x0,y0,z0,k,ia,np1,gs, nf,nf2,mode,i_0,i_1,kk,nf1_)
do j0=1,a%nb
  b=>a%bnd(j0)
  if (b%skip_getfun) cycle
  select case (b%boundLineType)
  case (1)
    !$omp parallel do if (gs_use_parallel_get_fun==1) reduction(+ : s) private (j,val,k1,ff) shared (j0,b,x0,y0,k,ia,np1,gs,mode,i_0,i_1,kk,nf1_)
    do j=1,b%npanel
      do k1=1,2
        val=b%psiom(j,k+k1)
        if (val.ne.d0) then
          select case (mode)
          case (1)
            ff=calcintm(j,i_1,nf1_(k1),ia,j0,i_0)
          case (2)
            ff=calcint_ar(j,i_1,nf1_(k1),ia,j0)
          case default
            ff=calcint(j,x0,y0,nf1_(k1),ia,j0)
          end select
          s=s+kk(k1)*val*ff
        endif
      enddo
    enddo
    !omp end parallel do nowait
  case (2)
    if (mode.ne.0) call gs_print_stop("Error get_fun_1")
    do j=1,b%nline
	    bl2=>b%line2(j)
	    ztc=bl2%zc-z0
	    ga=>bl2%ga(k+1)
	    do i=1,ga%n
	      val=ga%valkk(i)
		    if (val.ne.d0) s=s+val*cp_calcint(ztc,bl2%r,nf(i)+np1,bl2%g0,bl2%dir,0)
	    enddo	
	    ga=>bl2%ga(k+2)
	    do i=1,ga%n
	      val=ga%valkk(i)
		    if (val.ne.d0) s=s-val*cp_calcint(ztc,bl2%r,nf2(i)+np1,bl2%g0,bl2%dir,0)
	    enddo	
	  enddo
  end select
enddo
!omp end parallel do nowait
get_fun_1=s/pi2
end

subroutine get_fun_1_dual(x0,y0,k,ia,s,s2)
!вычисление производной функции - решения уравнения лапласа в области
use pgmod
!integer(4) mode=0!!!
real(8) x0,y0,s,s2,val,ff,ff2
real(8) fff(4,3) !для nf = 6,32,34;7,33,35;4,28,30;5,29,31
integer(4) j,k,j0,ia,k1,i
type(TArea), pointer :: a
type(TBound), pointer :: b
type(TBoundline2), pointer :: bl2
type(TBounLineFuncApprox), pointer :: ga
integer(4) nf1(2)
real(8) kk(2)
data nf1 /6,4/
data kk /d1,-d1/
complex(8) ztc,z0
a=>gs%a(ia)
s=d0
s2=d0
z0=dcmplx(x0,y0)
do j0=1,a%nb
  b=>a%bnd(j0)
  if (b%skip_getfun) cycle
  select case (b%boundLineType)
  case (1)
    !$omp parallel do if (gs_use_parallel_get_fun==1) reduction(+ : s,s2) private (j,val,k1,ff,ff2) shared (j0,b,x0,y0,k,ia,gs,kk,nf1)
    do j=1,b%npanel
      do k1=1,2
        val=b%psiom(j,k+k1)
        if (val.ne.d0) then
          call calcinta_2(j,x0,y0,nf1(k1),ia,j0,ff,ff2)
          s=s+val*kk(k1)*ff
          s2=s2+val*kk(k1)*ff2
        endif
      enddo
    enddo
    !omp end parallel do nowait
  end select
enddo
!$omp parallel do if (gs_use_parallel_get_fun==1) reduction(+ : s,s2) private (j0,j,b,val,i,bl2,ga,ztc,fff) shared (a,z0,k,gs)
do j0=1,a%nb
  b=>a%bnd(j0)
  if (b%skip_getfun) cycle
  select case (b%boundLineType)
  case (2)
    do j=1,b%nline
	    bl2=>b%line2(j)
	    ztc=bl2%zc-z0
	    ga=>bl2%ga(k+1)
      fff(1:2,1)=d0
      call cp_calcint_nf45_2(ztc,bl2%r,0,fff(3,1),fff(4,1))
      call cp_calcint_nf2831_2(ztc,bl2%r,0,fff(3,2),fff(4,2),fff(3,3),fff(4,3))
      call cp_calcint_nf3235_2(ztc,bl2%r,0,fff(1,2),fff(2,2),fff(1,3),fff(2,3))
	    do i=1,ga%n
	      val=ga%valkk(i)
		    if (val.ne.d0) then
          s=s+val*fff(1,i)
          s2=s2+val*fff(2,i)
        endif
	    enddo	
	    ga=>bl2%ga(k+2)
	    do i=1,ga%n
	      val=ga%valkk(i)
		    if (val.ne.d0) then
          s=s-val*fff(3,i)
          s2=s2-val*fff(4,i)
        endif
	    enddo	
	  enddo
  end select
enddo
!omp end parallel do nowait
s=s/pi2
s2=s2/pi2
end

function get_fun2(mode,x0,y0,i_0,i_1,ia)
!решение бигармонического уравнения в области
use pgmod
integer(4) mode !0 - (x0,y0), 1 - (i_0-bnd, i_1-panel), 2 - (i_1-ntr)
integer(4) i_0,i_1
real(8) x0,y0,s,val,ff
real(8) calcint,get_fun2,cp_calcint,calcintm,calcint_ar
integer(4) j,j0,ia,i,k
type(TArea), pointer :: a
type(TBound), pointer :: b
type(TBoundline2), pointer :: bl2
type(TBounLineFuncApprox), pointer :: ga
integer(4) nf0(3),nf1(3),nf2(3),nf3(3),nf(4)
real(8) kk(4),fff(3)
data nf0 /0,20,21/  !psi
data nf1 /1,22,23/  !psi'
data nf2 /2,24,25/  !eta
data nf3 /3,26,27/  !eta'
data nf /1,0,3,2/
data kk /d1,-d1,d1,-d1/
complex(8) ztc
a=>gs%a(ia)
s=d0
!$omp parallel do if (gs_use_parallel_get_fun==1) reduction(+ : s) private (j0,j,i,b,bl2,ztc,ga,val,k,ff,fff) shared (a,x0,y0,ia,nf0,nf1,nf2,nf3,gs,mode,i_0,i_1,kk,nf)
do j0=1,a%nb
  b=>a%bnd(j0)
  if (b%skip_getfun) cycle
  select case (b%boundLineType)
  case (1)
    !$omp parallel do if (gs_use_parallel_get_fun==1) reduction(+ : s) private (j,val,k,ff) shared (b,x0,y0,ia,j0,gs,mode,i_0,i_1,kk,nf)
    do j=1,b%npanel
      do k=1,4
        val=b%psiom(j,k)
        if (val.ne.d0) then
          select case (mode)
          case (1)
            ff=calcintm(j,i_1,nf(k),ia,j0,i_0)
          case (2)
            ff=calcint_ar(j,i_1,nf(k),ia,j0)
          case default
            ff=calcint(j,x0,y0,nf(k),ia,j0)
          end select
          s=s+kk(k)*val*ff
        endif
      enddo
    enddo
    !omp end parallel do nowait
  case (2)
    if (mode.ne.0) call gs_print_stop("Error get_fun2")
    do j=1,b%nline
	    bl2=>b%line2(j)
	    ztc=bl2%zc-dcmplx(x0,y0)
      if (gs_use_dual_integral_solve) then
	      ga=>bl2%ga(1)
        val=ga%valkk(1)
        if (val.ne.d0) s=s+val*cp_calcint(ztc,bl2%r,1,bl2%g0,bl2%dir,0)
        call cp_calcint_nf2223_2(ztc,bl2%r,0,fff(2),fff(3))
	      do i=2,ga%n
	        val=ga%valkk(i)
		      if (val.ne.d0) s=s+val*fff(i)
	      enddo	
	      ga=>bl2%ga(2)
        val=ga%valkk(1)
        if (val.ne.d0) s=s-val*cp_calcint(ztc,bl2%r,0,bl2%g0,bl2%dir,0)
        call cp_calcint_nf2021_2(ztc,bl2%r,0,fff(2),fff(3))
	      do i=2,ga%n
	        val=ga%valkk(i)
		      if (val.ne.d0) s=s-val*fff(i)
	      enddo	
	      ga=>bl2%ga(3)
        val=ga%valkk(1)
        if (val.ne.d0) s=s+val*cp_calcint(ztc,bl2%r,3,bl2%g0,bl2%dir,0)
        call cp_calcint_nf2627_2(ztc,bl2%r,0,fff(2),fff(3))
	      do i=2,ga%n
	        val=ga%valkk(i)
		      if (val.ne.d0) s=s+val*fff(i)
	      enddo	
	      ga=>bl2%ga(4)
        val=ga%valkk(1)
        if (val.ne.d0) s=s-val*cp_calcint(ztc,bl2%r,2,bl2%g0,bl2%dir,0)
        call cp_calcint_nf2425_2(ztc,bl2%r,0,fff(2),fff(3))
	      do i=2,ga%n
	        val=ga%valkk(i)
		      if (val.ne.d0) s=s-val*fff(i)
        enddo	
      else
        ga=>bl2%ga(1)
	      do i=1,ga%n
	        val=ga%valkk(i)
		      if (val.ne.d0) s=s+val*cp_calcint(ztc,bl2%r,nf1(i),bl2%g0,bl2%dir,0)
	      enddo	
	      ga=>bl2%ga(2)
	      do i=1,ga%n
	        val=ga%valkk(i)
		      if (val.ne.d0) s=s-val*cp_calcint(ztc,bl2%r,nf0(i),bl2%g0,bl2%dir,0)
	      enddo	
	      ga=>bl2%ga(3)
	      do i=1,ga%n
	        val=ga%valkk(i)
		      if (val.ne.d0) s=s+val*cp_calcint(ztc,bl2%r,nf3(i),bl2%g0,bl2%dir,0)
	      enddo	
	      ga=>bl2%ga(4)
	      do i=1,ga%n
	        val=ga%valkk(i)
		      if (val.ne.d0) s=s-val*cp_calcint(ztc,bl2%r,nf2(i),bl2%g0,bl2%dir,0)
        enddo	
      endif
	  enddo
  end select
enddo
!omp end parallel do nowait
get_fun2=s/pi2
end

function get_fun2_1(mode,x0,y0,i_0,i_1,np,ia)
!вычисление производной функции - решения бигармонического уравнения в области
!np - тип производной    1 - dpsi/dx, 2 - dpsi/dy
use pgmod
integer(4) mode !0 - (x0,y0), 1 - (i_0-bnd, i_1-panel), 2 - (i_1-ntr)
integer(4) i_0,i_1,np1
real(8) x0,y0,s,val,ff
real(8) calcint,get_fun2_1,cp_calcint,calcintm,calcint_ar
integer(4) j,np,j0,ia,i,k
type(TArea), pointer :: a
type(TBound), pointer :: b
type(TBoundline2), pointer :: bl2
type(TBounLineFuncApprox), pointer :: ga
integer(4) nf0(3),nf1(3),nf2(3),nf3(3),nf(4),nf_(4)
real(8) kk(4)
data nf0 /4,28,30/  !psi'
data nf1 /6,32,34/  !psi
data nf2 /8,36,38/  !eta'
data nf3 /10,40,42/  !eta
data nf /6,4,10,8/
data kk /d1,-d1,d1,-d1/
complex(8) ztc
np1=np-1
nf_=nf+np1
a=>gs%a(ia)
s=d0
!$omp parallel do if (gs_use_parallel_get_fun==1) reduction(+ : s) private (j0,j,b,val, i,bl2,ga,ztc,k,ff) shared (a,x0,y0,ia,np,gs, nf0,nf1,nf2,nf3,mode,i_0,i_1,kk,nf_)
do j0=1,a%nb
  b=>a%bnd(j0)
  if (b%skip_getfun) cycle
  select case (b%boundLineType)
  case (1)
    !$omp parallel do if (gs_use_parallel_get_fun==1) reduction(+ : s) private (j,val,k,ff) shared (j0,b,x0,y0,ia,np,gs,mode,i_0,i_1,kk,nf_)
    do j=1,b%npanel
      do k=1,4
        val=b%psiom(j,k)
        if (val.ne.d0) then
          select case (mode)
          case (1)
            ff=calcintm(j,i_1,nf_(k),ia,j0,i_0)
          case (2)
            ff=calcint_ar(j,i_1,nf_(k),ia,j0)
          case default
            ff=calcint(j,x0,y0,nf_(k),ia,j0)
          end select
          s=s+kk(k)*val*ff
        endif
      enddo
    enddo
    !omp end parallel do nowait
  case (2)
    if (mode.ne.0) call gs_print_stop("Error get_fun2_1")
    do j=1,b%nline
	    bl2=>b%line2(j)
	    ztc=bl2%zc-dcmplx(x0,y0)
	    ga=>bl2%ga(1)
	    do i=1,ga%n
	      val=ga%valkk(i)
		    if (val.ne.d0) s=s+val*cp_calcint(ztc,bl2%r,nf1(i)+np1,bl2%g0,bl2%dir,0)
	    enddo	
	    ga=>bl2%ga(2)
	    do i=1,ga%n
	      val=ga%valkk(i)
		    if (val.ne.d0) s=s-val*cp_calcint(ztc,bl2%r,nf0(i)+np1,bl2%g0,bl2%dir,0)
	    enddo	
	    ga=>bl2%ga(3)
	    do i=1,ga%n
	      val=ga%valkk(i)
		    if (val.ne.d0) s=s+val*cp_calcint(ztc,bl2%r,nf3(i)+np1,bl2%g0,bl2%dir,0)
	    enddo	
	    ga=>bl2%ga(4)
	    do i=1,ga%n
	      val=ga%valkk(i)
		    if (val.ne.d0) s=s-val*cp_calcint(ztc,bl2%r,nf2(i)+np1,bl2%g0,bl2%dir,0)
	    enddo	
	  enddo
  end select
enddo
!omp end parallel do nowait
get_fun2_1=s/pi2
end

subroutine get_fun2_1_dual(x0,y0,ia,s,s2)
!вычисление производной функции - решения бигармонического уравнения в области
use pgmod
!integer(4) mode=0!!!
real(8) x0,y0,s,val,ff,ff2,s2
integer(4) j,j0,ia,k,i
type(TArea), pointer :: a
type(TBound), pointer :: b
type(TBoundline2), pointer :: bl2
type(TBounLineFuncApprox), pointer :: ga
integer(4) nf(4)
real(8) kk(4)
real(8) fff(8,3) !для nf = 6,32,34;7,33,35;4,28,30;5,29,31;10,40,42;11,41,43;8,36,38;9,37,39
data nf /6,4,10,8/
data kk /d1,-d1,d1,-d1/
complex(8) ztc,z0
a=>gs%a(ia)
s=d0
s2=0
z0=dcmplx(x0,y0)
do j0=1,a%nb
  b=>a%bnd(j0)
  if (b%skip_getfun) cycle
  select case (b%boundLineType)
  case (1)
    !$omp parallel do if (gs_use_parallel_get_fun==1) reduction(+ : s,s2) private (j,val,k,ff,ff2) shared (j0,b,x0,y0,ia,gs,kk,nf)
    do j=1,b%npanel
      do k=1,4
        val=b%psiom(j,k)
        if (val.ne.d0) then
          call calcinta_2(j,x0,y0,nf(k),ia,j0,ff,ff2)
          s=s+kk(k)*val*ff
          s2=s2+kk(k)*val*ff2
        endif
      enddo
    enddo
    !omp end parallel do nowait
  end select
enddo
!$omp parallel do if (gs_use_parallel_get_fun==1) reduction(+ : s,s2) private (j0,j,b,val, i,bl2,ga,ztc,fff) shared (a,z0,gs)
do j0=1,a%nb
  b=>a%bnd(j0)
  if (b%skip_getfun) cycle
  select case (b%boundLineType)
  case (2)
    do j=1,b%nline
	    bl2=>b%line2(j)
	    ztc=bl2%zc-z0
      !для nf = 6,32,34;7,33,35;4,28,30;5,29,31; 10,40,42;11,41,43;8,36,38;9,37,39
      fff(1:2,1)=d0
      call cp_calcint_nf45_2(ztc,bl2%r,0,fff(3,1),fff(4,1))
      call cp_calcint_nf2831_2(ztc,bl2%r,0,fff(3,2),fff(4,2),fff(3,3),fff(4,3))
      call cp_calcint_nf3235_2(ztc,bl2%r,0,fff(1,2),fff(2,2),fff(1,3),fff(2,3))
      call cp_calcint_nf89_2(ztc,bl2%r,0,fff(7,1),fff(8,1))
      call cp_calcint_nf1011_2(ztc,bl2%r,0,fff(5,1),fff(6,1))
      call cp_calcint_nf3639_2(ztc,bl2%r,0,fff(7,2),fff(8,2),fff(7,3),fff(8,3))
      call cp_calcint_nf4043_2(ztc,bl2%r,0,fff(5,2),fff(6,2),fff(5,3),fff(6,3))
      ga=>bl2%ga(1)
	    do i=1,ga%n
	      val=ga%valkk(i)
		    if (val.ne.d0) then
          s=s+val*fff(1,i)
          s2=s2+val*fff(2,i)
        endif
	    enddo	
	    ga=>bl2%ga(2)
	    do i=1,ga%n
	      val=ga%valkk(i)
		    if (val.ne.d0) then
          s=s-val*fff(3,i)
          s2=s2-val*fff(4,i)
        endif
	    enddo	
	    ga=>bl2%ga(3)
	    do i=1,ga%n
	      val=ga%valkk(i)
		    if (val.ne.d0) then
          s=s+val*fff(5,i)
          s2=s2+val*fff(6,i)
        endif
	    enddo	
	    ga=>bl2%ga(4)
	    do i=1,ga%n
	      val=ga%valkk(i)
		    if (val.ne.d0) then
          s=s-val*fff(7,i)
          s2=s2-val*fff(8,i)
        endif
	    enddo	
	  enddo
  end select
enddo
!omp end parallel do nowait
s=s/pi2
s2=s2/pi2
end

function get_funbri(mode,x0,y0,i_0,i_1,ia)
!решение уравнения Бринкмана в области через единую функцию Грина
use pgmod
integer(4) mode !0 - (x0,y0), 1 - (i_0-bnd, i_1-panel), 2 - (i_1-ntr)
integer(4) i_0,i_1
real(8) x0,y0,s,val,k2,ff
real(8) calcint,get_funbri,calcintn,calcintm,calcint_ar,calcintn2,calcintn_ar
integer(4) j,j0,ia,k
type(TArea), pointer :: a
type(TBound), pointer :: b
integer(4) nf(4)
real(8) kk(4)
data nf /1,0,15,14/
data kk /d1,-d1,d1,-d1/
a=>gs%a(ia)
s=d0
k2=d1/a%const%k_helm**2
!$omp parallel do if (gs_use_parallel_get_fun==1) reduction(+ : s) private (j0,j,b,val,k,ff) shared (x0,y0,a,ia,k2,gs,mode,i_0,i_1,nf,kk)
do j0=1,a%nb
  b=>a%bnd(j0)
  if (b%skip_getfun) cycle
  !$omp parallel do if (gs_use_parallel_get_fun==1) reduction(+ : s) private (j,val,k,ff) shared (j0,b,x0,y0,a,ia,k2,gs,mode,i_0,i_1,nf,kk)
  do j=1,b%npanel
    do k=1,2
      val=b%psiom(j,k)-b%psiom(j,k+2)*k2
      if (val.ne.d0) then
        select case (mode)
        case (1)
          ff=calcintm(j,i_1,nf(k),ia,j0,i_0)
        case (2)
          ff=calcint_ar(j,i_1,nf(k),ia,j0)
        case default
          ff=calcint(j,x0,y0,nf(k),ia,j0)
        end select
        s=s+kk(k)*val*ff
      endif
    enddo
    do k=3,4
      val=b%psiom(j,k)*k2
      if (val.ne.d0) then
        select case (mode)
        case (1)
          ff=calcintn2(j,i_1,nf(k),ia,j0,i_0)
        case (2)
          ff=calcintn_ar(j,i_1,nf(k),ia,j0)
        case default
          ff=calcintn(j,x0,y0,nf(k),ia,j0)
        end select
        s=s+kk(k)*val*ff
      endif
    enddo
  enddo
  !omp end parallel do nowait
enddo
!omp end parallel do nowait
get_funbri=s/pi2
end

function get_funbri_1(mode,x0,y0,i_0,i_1,np,ia)
!вычисление производной функции - решения уравнения Бринкмана в области через единую функцию Грина
!np - тип производной    1 - dpsi/dx, 2 - dpsi/dy
use pgmod
integer(4) mode !0 - (x0,y0), 1 - (i_0-bnd, i_1-panel), 2 - (i_1-ntr)
integer(4) i_0,i_1
real(8) x0,y0,s,val,k2,ff
real(8) calcint,get_funbri_1,calcintn,calcintm,calcint_ar,calcintn2,calcintn_ar
integer(4) j,j0,ia,np,k
type(TArea), pointer :: a
type(TBound), pointer :: b
integer(4) nf(4),nf_(4)
real(8) kk(4)
data nf_ /5,3,17,15/
data kk /d1,-d1,d1,-d1/
nf=nf_+np
a=>gs%a(ia)
s=d0
k2=d1/a%const%k_helm**2
!$omp parallel do if (gs_use_parallel_get_fun==1) reduction(+ : s) private (j0,j,b,val,k,ff) shared (x0,y0,a,ia,k2,nf,gs,mode,i_0,i_1,kk)
do j0=1,a%nb
  b=>a%bnd(j0)
  if (b%skip_getfun) cycle
  !$omp parallel do if (gs_use_parallel_get_fun==1) reduction(+ : s) private (j,val,k,ff) shared (j0,b,x0,y0,a,ia,k2,nf,gs,mode,i_0,i_1,kk)
  do j=1,b%npanel
    do k=1,2
      val=b%psiom(j,k)-b%psiom(j,k+2)*k2
      if (val.ne.d0) then
        select case (mode)
        case (1)
          ff=calcintm(j,i_1,nf(k),ia,j0,i_0)
        case (2)
          ff=calcint_ar(j,i_1,nf(k),ia,j0)
        case default
          ff=calcint(j,x0,y0,nf(k),ia,j0)
        end select
        s=s+kk(k)*val*ff
      endif
    enddo
    do k=3,4
      val=b%psiom(j,k)*k2
      if (val.ne.d0) then
        select case (mode)
        case (1)
          ff=calcintn2(j,i_1,nf(k),ia,j0,i_0)
        case (2)
          ff=calcintn_ar(j,i_1,nf(k),ia,j0)
        case default
          ff=calcintn(j,x0,y0,nf(k),ia,j0)
        end select
        s=s+kk(k)*val*ff
      endif
    enddo
  enddo
  !omp end parallel do nowait
enddo
!omp end parallel do nowait
get_funbri_1=s/pi2
end

subroutine get_funbri_1_dual(x0,y0,ia,s,s2)
!вычисление производной функции - решения уравнения Бринкмана в области через единую функцию Грина
use pgmod
!integer(4) mode=0
real(8) x0,y0,s,val,k2,ff,s2,ff2
integer(4) j,j0,ia,k
type(TArea), pointer :: a
type(TBound), pointer :: b
integer(4) nf(4)
real(8) kk(4)
data nf/6,4,18,16/
data kk /d1,-d1,d1,-d1/
a=>gs%a(ia)
s=d0
s2=d0
k2=d1/a%const%k_helm**2
!!$omp parallel do if (gs_use_parallel_get_fun==1) reduction(+ : s,s2) private (j0,j,b,val,k,ff,ff2) shared (x0,y0,a,ia,k2,nf,kk)
!$omp parallel do if (gs_use_parallel_get_fun==1) private (j0,j,b,val,k,ff,ff2) shared (x0,y0,a,ia,k2,nf,kk)
do j0=1,a%nb
  b=>a%bnd(j0)
  if (b%skip_getfun) cycle
  !$omp parallel do if (gs_use_parallel_get_fun==1) reduction(+ : s,s2) private (j,val,k,ff,ff2) shared (j0,b,x0,y0,a,ia,k2,nf,kk)
  do j=1,b%npanel
    do k=1,2
      val=b%psiom(j,k)-b%psiom(j,k+2)*k2
      if (val.ne.d0) then
        call calcinta_2(j,x0,y0,nf(k),ia,j0,ff,ff2)  
        s=s+kk(k)*val*ff
        s2=s2+kk(k)*val*ff2
      endif
    enddo
    do k=3,4
      val=b%psiom(j,k)*k2
      if (val.ne.d0) then
        call calcintn_dual(j,x0,y0,nf(k),ia,j0,ff,ff2)
        s=s+kk(k)*val*ff
        s2=s2+kk(k)*val*ff2
      endif
    enddo
  enddo
  !omp end parallel do nowait
enddo
!omp end parallel do nowait
s=s/pi2
s2=s2/pi2
end

function get_funp(mode,x0,y0,i_0,i_1,k,ia)
!решение уравнения Пуассона, Гельмгольца или переноса в области
use pgmod
integer(4) mode !0 - (x0,y0), 1 - (i_0-bnd, i_1-panel), 2 - (i_1-ntr)
integer(4) i_0,i_1
real(8) x0,y0,s,val,ff
real(8) get_funp,get_fun,intd_z,intd_m,intd_ar
integer(4) j,ia,k,kp
complex(8) z0
type(areatype), pointer :: a
type(TAreaValue), pointer :: av
a=>gs%a(ia)%a
s=d0
kp=k/2
z0=dcmplx(x0,y0)
av=>a%arval_funp(1+kp)
!$omp parallel do if (gs_use_parallel_get_fun==1) reduction(+ : s) private (j,val,ff) shared (a,av,z0,ia,gs,mode,i_0,i_1)
do j=1,a%ntr
  val=av%v(j)
  if (val.ne.d0) then
    select case (mode)
    case (1)
      ff=intd_m(ia,j,i_1,i_0,0)
    case (2)
      ff=intd_ar(ia,j,i_1,0)
    case default
      ff=intd_z(z0,ia,j,0) !intd_an(z0,zmz(1,j,ia),zmz(2,j,ia),zmz(3,j,ia),0)
    end select
    s=s+val*ff
  endif
enddo
!omp end parallel do nowait
get_funp=s/pi2+get_fun(mode,x0,y0,i_0,i_1,k,ia)
end

function get_funp_1(mode,x0,y0,i_0,i_1,k,np,ia)
!вычисление производной функции - решения уравнения Пуассона, Гельмгольца или переноса в области
!np - тип производной    1 - dpsi/dx, 2 - dpsi/dy
use pgmod
integer(4) mode !0 - (x0,y0), 1 - (i_0-bnd, i_1-panel), 2 - (i_1-ntr)
integer(4) i_0,i_1
real(8) x0,y0,s,val,ff
real(8) get_funp_1,get_fun_1,intd_z,intd_m,intd_ar
integer(4) j,ia,k,kp,np,nf
complex(8) z0
type(areatype), pointer :: a
type(TAreaValue), pointer :: av
a=>gs%a(ia)%a
s=d0
kp=k/2
nf=np+3
z0=dcmplx(x0,y0)
av=>a%arval_funp(1+kp)
!$omp parallel do if (gs_use_parallel_get_fun==1) reduction(+ : s) private (j,val,ff) shared (a,av,z0,ia,nf,gs,mode,i_0,i_1)
do j=1,a%ntr
  val=av%v(j)
  if (val.ne.d0) then
    select case (mode)
    case (1)
      ff=intd_m(ia,j,i_1,i_0,nf)
    case (2)
      ff=intd_ar(ia,j,i_1,nf)
    case default
      ff=intd_z(z0,ia,j,nf) !intd_an(z0,zmz(1,j,ia),zmz(2,j,ia),zmz(3,j,ia),np+3)
    end select
    s=s+val*ff
  endif
enddo
!omp end parallel do nowait
get_funp_1=s/pi2+get_fun_1(mode,x0,y0,i_0,i_1,k,np,ia)
end

subroutine get_funp_1_dual(x0,y0,k,ia,s,s2)
!вычисление производной функции - решения уравнения Пуассона, Гельмгольца или переноса в области
use pgmod
!integer(4) mode=0!!!
real(8) x0,y0,s,s2,val,ff,ff2
integer(4) j,ia,k,kp,nf
complex(8) z0
type(areatype), pointer :: a
type(TAreaValue), pointer :: av
a=>gs%a(ia)%a
s=d0
s2=d0
kp=k/2
nf=4
z0=dcmplx(x0,y0)
av=>a%arval_funp(1+kp)
!$omp parallel do if (gs_use_parallel_get_fun==1) reduction(+ : s,s2) private (j,val,ff,ff2) shared (a,av,z0,ia,nf)
do j=1,a%ntr
  val=av%v(j)
  if (val.ne.d0) then
    call intd_z_2(z0,ia,j,nf,ff,ff2)
    s=s+val*ff
    s2=s2+val*ff2
  endif
enddo
!omp end parallel do nowait
call get_fun_1_dual(x0,y0,k,ia,ff,ff2)
s=s/pi2+ff
s2=s2/pi2+ff2
end

function get_fun2p(mode,x0,y0,i_0,i_1,ia)
!решение неоднородного бигармонического уравнения
use pgmod
integer(4) mode !0 - (x0,y0), 1 - (i_0-bnd, i_1-panel), 2 - (i_1-ntr)
integer(4) i_0,i_1
real(8) x0,y0,s,val,ff
real(8) get_fun2p,get_fun2,intd_z,intd_m,intd_ar
integer(4) j,ia
complex(8) z0
type(areatype), pointer :: a
type(TAreaValue), pointer :: av
a=>gs%a(ia)%a
s=d0
z0=dcmplx(x0,y0)
av=>a%arval_funp(1)
!$omp parallel do if (gs_use_parallel_get_fun==1) reduction(+ : s) private (j,val,ff) shared (a,z0,ia,gs,mode,i_0,i_1,av)
do j=1,a%ntr
  val=av%v(j)
  if (val.ne.d0) then
    select case (mode)
    case (1)
      ff=intd_m(ia,j,i_1,i_0,2)
    case (2)
      ff=intd_ar(ia,j,i_1,2)
    case default
      ff=intd_z(z0,ia,j,2) !intd_an(z0,zmz(1,j,ia),zmz(2,j,ia),zmz(3,j,ia),2)
    end select
    s=s+val*ff
  endif
enddo
!omp end parallel do nowait
get_fun2p=s/pi2+get_fun2(mode,x0,y0,0,0,ia)
end

function get_fun2p_1(mode,x0,y0,i_0,i_1,np,ia)
!вычисление производной функции - решения неоднородного бигармонического уравнения
!np - тип производной    1 - dpsi/dx, 2 - dpsi/dy
use pgmod
integer(4) mode !0 - (x0,y0), 1 - (i_0-bnd, i_1-panel), 2 - (i_1-ntr)
integer(4) i_0,i_1
real(8) x0,y0,s,val,ff
real(8) get_fun2p_1,get_fun2_1,intd_z,intd_m,intd_ar
integer(4) j,ia,np,nf
complex(8) z0
type(areatype), pointer :: a
type(TAreaValue), pointer :: av
a=>gs%a(ia)%a
s=d0
z0=dcmplx(x0,y0)
nf=np+7
av=>a%arval_funp(1)
!$omp parallel do if (gs_use_parallel_get_fun==1) reduction(+ : s) private (j,val,ff) shared (a,z0,ia,nf,gs,mode,i_0,i_1,av)
do j=1,a%ntr
  val=av%v(j)
  if (val.ne.d0) then
    select case (mode)
    case (1)
      ff=intd_m(ia,j,i_1,i_0,nf)
    case (2)
      ff=intd_ar(ia,j,i_1,nf)
    case default
      ff=intd_z(z0,ia,j,nf) !intd_an(z0,zmz(1,j,ia),zmz(2,j,ia),zmz(3,j,ia),np+7)
    end select
    s=s+val*ff
  endif
enddo
!omp end parallel do nowait
get_fun2p_1=s/pi2+get_fun2_1(mode,x0,y0,i_0,i_1,np,ia)
end

subroutine get_fun2p_1_dual(x0,y0,ia,s,s2)
!вычисление производной функции - решения неоднородного бигармонического уравнения
use pgmod
!integer(4) mode=0
real(8) x0,y0,s,val,ff,ff2,s2
integer(4) j,ia,nf
complex(8) z0
type(areatype), pointer :: a
type(TAreaValue), pointer :: av
a=>gs%a(ia)%a
s=d0
s2=d0
z0=dcmplx(x0,y0)
nf=8
av=>a%arval_funp(1)
!$omp parallel do if (gs_use_parallel_get_fun==1) reduction(+ : s,s2) private (j,val,ff,ff2) shared (a,z0,ia,nf,av)
do j=1,a%ntr
  val=av%v(j)
  if (val.ne.d0) then
    call intd_z_2(z0,ia,j,nf,ff,ff2)
    s=s+val*ff
    s2=s2+val*ff2
  endif
enddo
!omp end parallel do nowait
call get_fun2_1_dual(x0,y0,ia,ff,ff2)
s=s/pi2+ff
s2=s2/pi2+ff2
end

function get_funoss(mode,x0,y0,i_0,i_1,ia)
!решение осесимметричного уравнения в области через вычисление интеграла с учетом 1/y
use pgmod
integer(4) mode !0 - (x0,y0), 1 - (i_0-bnd, i_1-panel), 2 - (i_1-ntr)
integer(4) i_0,i_1
real(8) x0,y0,s,val,ff
real(8) get_funoss,get_fun,intd_oss_z,intd_oss_m,intd_oss_ar
integer(4) j,ia
complex(8) z0
type(areatype), pointer :: a
type(TAreaValue), pointer :: av
a=>gs%a(ia)%a
s=d0
z0=dcmplx(x0,y0)
av=>a%arval_funp(1)
!$omp parallel do if (gs_use_parallel_get_fun==1) reduction(+ : s) private (j,val,ff) shared (a,z0,ia,gs,mode,i_0,i_1,av)
do j=1,a%ntr
  val=av%v(j)
  if (val.ne.d0) then
    select case (mode)
    case (1)
      ff=intd_oss_m(ia,j,i_1,i_0,0)
    case (2)
      ff=intd_oss_ar(ia,j,i_1,0)
    case default
      ff=intd_oss_z(z0,ia,j,0) !intd_oss(z0,zmz(1,j,ia),zmz(2,j,ia),zmz(3,j,ia),0)
    end select
    s=s+val*ff
  endif
enddo
!omp end parallel do nowait
get_funoss=s/pi2+get_fun(mode,x0,y0,i_0,i_1,0,ia)
end

function get_funoss_1(mode,x0,y0,i_0,i_1,np,ia)
!решение осесимметричного уравнения в области через вычисление интеграла с учетом 1/y
!np - тип производной    1 - dpsi/dx, 2 - dpsi/dy
use pgmod
integer(4) mode !0 - (x0,y0), 1 - (i_0-bnd, i_1-panel), 2 - (i_1-ntr)
integer(4) i_0,i_1
real(8) x0,y0,s,val,ff
real(8) get_funoss_1,get_fun_1,intd_oss_z,intd_oss_m,intd_oss_ar
integer(4) j,ia,np,nf
complex(8) z0
type(areatype), pointer :: a
type(TAreaValue), pointer :: av
a=>gs%a(ia)%a
s=d0
z0=dcmplx(x0,y0)
nf=np+3
av=>a%arval_funp(1)
!$omp parallel do if (gs_use_parallel_get_fun==1) reduction(+ : s) private (j,val,ff) shared (a,z0,ia,nf,gs,mode,i_0,i_1,av)
do j=1,a%ntr
  val=av%v(j)
  if (val.ne.d0) then
    select case (mode)
    case (1)
      ff=intd_oss_m(ia,j,i_1,i_0,nf)
    case (2)
      ff=intd_oss_ar(ia,j,i_1,nf)
    case default
      ff=intd_oss_z(z0,ia,j,nf) !intd_oss(z0,zmz(1,j,ia),zmz(2,j,ia),zmz(3,j,ia),np+3)
    end select
    s=s+val*ff
  endif
enddo
!omp end parallel do nowait
get_funoss_1=s/pi2+get_fun_1(mode,x0,y0,i_0,i_1,0,np,ia)
end

function get_funh(mode,x0,y0,i_0,i_1,k,ia)
!решение уравнения Гельмгольца через функции Бесселя в области \Delta\psi+k_helm^2\psi=0
use pgmod
integer(4) mode !0 - (x0,y0), 1 - (i_0-bnd, i_1-panel), 2 - (i_1-ntr)
integer(4) i_0,i_1
real(8) x0,y0,s,val,ff
real(8) calcintn,get_funh,calcintn2,calcintn_ar
integer(4) j,k,ia,j0,k1
type(TBound), pointer :: b
type(TArea), pointer :: a
integer(4) nf(2)
real(8) kk(2)
data nf /13,12/
data kk /d1,-d1/
a=>gs%a(ia)
s=d0
!$omp parallel do if (gs_use_parallel_get_fun==1) reduction(+ : s) private (j,j0,b,val,k1,ff) shared (a,x0,y0,ia,k,gs,mode,i_0,i_1,nf,kk)
do j0=1,a%nb
  b=>a%bnd(j0)
  if (b%skip_getfun) cycle
  !$omp parallel do if (gs_use_parallel_get_fun==1) reduction(+ : s) private (j,val,k1,ff) shared (a,j0,b,x0,y0,ia,k,gs,mode,i_0,i_1,nf,kk)
  do j=1,b%npanel
    do k1=1,2
      val=b%psiom(j,k+k1)
      if (val.ne.d0) then
        select case (mode)
        case (1)
          ff=calcintn2(j,i_1,nf(k1),ia,j0,i_0)
        case (2)
          ff=calcintn_ar(j,i_1,nf(k1),ia,j0)
        case default
          ff=calcintn(j,x0,y0,nf(k1),ia,j0)
        end select
        s=s+kk(k1)*val*ff
      endif
	  enddo
  enddo
  !omp end parallel do nowait
enddo
!omp end parallel do nowait
get_funh=s/pi2
end

function get_funh2(mode,x0,y0,i_0,i_1,k,ia)
!решение уравнения Гельмгольца через функции Бесселя в области \Delta\psi-k_helm^2\psi=0
use pgmod
integer(4) mode !0 - (x0,y0), 1 - (i_0-bnd, i_1-panel), 2 - (i_1-ntr)
integer(4) i_0,i_1
real(8) x0,y0,s,val,ff
real(8) calcintn,get_funh2,calcintn2,calcintn_ar
integer(4) j,k,ia,j0,k1
type(TBound), pointer :: b
type(TArea), pointer :: a
integer(4) nf(2)
real(8) kk(2)
data nf /15,14/
data kk /d1,-d1/
a=>gs%a(ia)
s=d0
!$omp parallel do if (gs_use_parallel_get_fun==1) reduction(+ : s) private (j,j0,b,val,k1,ff) shared (a,x0,y0,ia,k,mode,i_0,i_1,nf,kk)
do j0=1,a%nb
  b=>a%bnd(j0)
  if (b%skip_getfun) cycle
  !$omp parallel do if (gs_use_parallel_get_fun==1) reduction(+ : s) private (j,val,k1,ff) shared (a,j0,b,x0,y0,ia,k,mode,i_0,i_1,nf,kk)
  do j=1,b%npanel
    do k1=1,2
      val=b%psiom(j,k+k1)
      if (val.ne.d0) then
        select case (mode)
        case (1)
          ff=calcintn2(j,i_1,nf(k1),ia,j0,i_0)
        case (2)
          ff=calcintn_ar(j,i_1,nf(k1),ia,j0)
        case default
          ff=calcintn(j,x0,y0,nf(k1),ia,j0)
        end select
        s=s+kk(k1)*val*ff
      endif
	  enddo
  enddo
  !omp end parallel do nowait
enddo
!omp end parallel do nowait
get_funh2=s/pi2
end

function get_funh2_1(mode,x0,y0,i_0,i_1,k,np,ia)
!!вычисление производной функции - решение уравнения Гельмгольца через функции Бесселя в области \Delta\psi-k_helm^2\psi=0
!np - тип производной    1 - dpsi/dx, 2 - dpsi/dy
use pgmod
integer(4) mode !0 - (x0,y0), 1 - (i_0-bnd, i_1-panel), 2 - (i_1-ntr)
integer(4) i_0,i_1,np
real(8) x0,y0,s,val,ff
real(8) calcintn,get_funh2_1,calcintn2,calcintn_ar
integer(4) j,k,ia,j0,k1
type(TBound), pointer :: b
type(TArea), pointer :: a
integer(4) nf(2),nf_(2)
real(8) kk(2)
data nf_ /17,15/
data kk /d1,-d1/
nf=nf_+np
a=>gs%a(ia)
s=d0
!$omp parallel do if (gs_use_parallel_get_fun==1) reduction(+ : s) private (j,j0,b,val,k1,ff) shared (a,x0,y0,ia,k,mode,i_0,i_1,nf,kk)
do j0=1,a%nb
  b=>a%bnd(j0)
  if (b%skip_getfun) cycle
  !$omp parallel do if (gs_use_parallel_get_fun==1) reduction(+ : s) private (j,val,k1,ff) shared (a,j0,b,x0,y0,ia,k,mode,i_0,i_1,nf,kk)
  do j=1,b%npanel
    do k1=1,2
      val=b%psiom(j,k+k1)
      if (val.ne.d0) then
        select case (mode)
        case (1)
          ff=calcintn2(j,i_1,nf(k1),ia,j0,i_0)
        case (2)
          ff=calcintn_ar(j,i_1,nf(k1),ia,j0)
        case default
          ff=calcintn(j,x0,y0,nf(k1),ia,j0)
        end select
        s=s+kk(k1)*val*ff
      endif
	  enddo
  enddo
  !omp end parallel do nowait
enddo
!omp end parallel do nowait
get_funh2_1=s/pi2
end

subroutine get_funh2_1_dual(x0,y0,k,ia,s,s2)
!вычисление производной функции - решения уравнения лапласа в области
use pgmod
!integer(4) mode=0!!!
real(8) x0,y0,s,s2,val,ff,ff2
integer(4) j,k,j0,ia,k1
type(TArea), pointer :: a
type(TBound), pointer :: b
integer(4) nf(2)
real(8) kk(2)
data nf /18,16/
data kk /d1,-d1/
a=>gs%a(ia)
s=d0
s2=d0
do j0=1,a%nb
  b=>a%bnd(j0)
  if (b%skip_getfun) cycle
  select case (b%boundLineType)
  case (1)
    !$omp parallel do if (gs_use_parallel_get_fun==1) reduction(+ : s,s2) private (j,val,k1,ff,ff2) shared (j0,b,x0,y0,k,ia,kk,nf)
    do j=1,b%npanel
      do k1=1,2
        val=b%psiom(j,k+k1)
        if (val.ne.d0) then
          call calcintn_dual(j,x0,y0,nf(k1),ia,j0,ff,ff2)
          s=s+val*kk(k1)*ff
          s2=s2+val*kk(k1)*ff2
        endif
      enddo
    enddo
    !omp end parallel do nowait
  end select
enddo
s=s/pi2
s2=s2/pi2
end

subroutine solve_postprocessor
use pgmod
integer(4) type_eq,i0,type_eq2,i,j,k,k1
type(TArea), pointer :: a
type(TBound), pointer :: b
type(TBoundline2), pointer :: bl2
type(TBounLineFuncApprox), pointer :: ga
real(8) sigm,kap
do i0=1,gs%na
  a=>gs%a(i0)
  type_eq=a%type_eq(1)
  select case (type_eq)
  case (4:6,8,16:22,24:26)
    call get_arval_funp(type_eq,i0)
    select case (type_eq)
	  case (8)
	    type_eq2=a%type_eq(2)
	    if (type_eq2==7.or.type_eq2==14) call get_arval_funp(type_eq2,i0)
    case (21)
	    type_eq2=a%type_eq(2)
	    if (type_eq2==7) call get_arval_funp(type_eq2,i0)
    case (22,24,25)
	    type_eq2=a%type_eq(2)
	    if (type_eq2==23) call get_arval_funp(type_eq2,i0)
	  end select
  case (11)
    type_eq2=a%type_eq(2)
    select case (type_eq2)
    case (2,12,13)
	    call get_arval_funp2(type_eq,type_eq2,i0)
    end select
  end select
  do i=1,a%nb
    b=>a%bnd(i)
	  if (b%boundLineType==2) then
	    do j=1,b%nline
	      bl2=>b%line2(j)
	  	  do k=1,a%umax
	  	    ga=>bl2%ga(k)
	  	    if (ga%typea==1) then
	  	      do k1=1,ga%n
	  	  	    ga%valkk(k1)=ga%val(k1)*ga%kk(k1)
	  	  	  enddo
	  	      ga%cc=dsqrt(ga%valkk(2)**2+ga%valkk(3)**2)
            if (ga%cc>d0) then
              ga%delta=datan2(ga%valkk(3),ga%valkk(2))-pi
              if (ga%delta<-pi) ga%delta=ga%delta+pi2
            endif
	  	    endif
        enddo
        !!!эти формулы работают более менее только в области, где внешняя (осредненная) завихренность близка к 0
        sigm=bl2%r*bl2%ga(4)%cc/bl2%ga(3)%cc
        bl2%m1=(sigm-d1)/(sigm+d1)
        kap=bl2%m1-bl2%m1**2-d5*dlog(bl2%m1)-0.75d0
        bl2%u=bl2%ga(3)%cc*bl2%r*kap*d5/(d1-bl2%m1)
        bl2%tet=(bl2%ga(3)%delta+bl2%ga(4)%delta)*d5
	    enddo
	  endif
  enddo
enddo
end

subroutine get_arval_funp(ku,knd)
use pgmod
integer(4) ku,knd
call get_arval_funp2(ku,0,knd)
end

subroutine do_get_arval_funp2(knd,kp,iav,ipsiar,mode)
use pgmod
integer(4) knd,kp,iav,ipsiar,j
integer(4) mode !1 "+", -1 "-"
type(TAreaValue), pointer :: av,avfp
type(areatype), pointer :: a
real(8), pointer :: psiar(:)
real(8) val
a=>gs%a(knd)%a
avfp=>a%arval_funp(kp)
av=>a%areaval(iav)
psiar=>a%psiarea(:,ipsiar)
if (mode>0) then
  !$omp parallel do if (gs_use_parallel_get_fun==1) private (j,val) shared (a,av,avfp,psiar)
  do j=1,a%ntr 
    val=av%v(j)
    if (val/=d0) avfp%v(j)=avfp%v(j)+val*psiar(j)
  enddo
  !omp end parallel do nowait
else
  !$omp parallel do if (gs_use_parallel_get_fun==1) private (j,val) shared (a,av,avfp,psiar)
  do j=1,a%ntr 
    val=av%v(j)
    if (val/=d0) avfp%v(j)=avfp%v(j)-val*psiar(j)
  enddo
  !omp end parallel do nowait
endif
end

subroutine get_arval_funp2(ku1,ku2,knd)
use pgmod
integer(4) ku1,ku2,j,kp,knd,i
real(8) get_fun,get_funh,get_funh2,val,fff
type(areatype), pointer :: a
type(TArea), pointer :: ar
type(TAreaValue), pointer :: av,avfp
ar=>gs%a(knd)
a=>ar%a
kp=1
if ((ku1.eq.7).or.(ku1.eq.14).or.(ku1.eq.23)) kp=2 
avfp=>a%arval_funp(kp)
av=>null()
select case (ku1) !неоднородные уравнения
case (4,6,7,17,18,21)
  av=>a%areaval(1)
case (22,24:26)
  if (ar%eq_var(1)) av=>a%areaval(1)
end select
if (associated(av)) then
  avfp%v=av%v
else
  avfp%v=d0
endif
select case (ku1)
case (5,6,14)
  call do_get_arval_funp2(knd,kp,2,1,-1)
case (8)
  call do_get_arval_funp2(knd,kp,3,1,1)
case (11)
  av=>a%areaval(3)
  !$omp parallel do if (gs_use_parallel_get_fun==1) private (j,val,fff) shared (a,av,avfp,ku2,knd)
  do j=1,a%ntr 
    val=av%v(j)
	  if (val.ne.d0) then 
      select case (ku2)
      case (2)
        fff=get_fun(2,d0,d0,0,j,2,knd)
      case (12)
	      fff=get_funh(2,d0,d0,0,j,2,knd)
	    case (13)
	      fff=get_funh2(2,d0,d0,0,j,2,knd)
	    end select
	    avfp%v(j)=val*fff
    endif
  enddo
  !omp end parallel do nowait
case (16:20)
  call do_get_arval_funp2(knd,kp,4,1,1)
  if ((ku1.eq.16).or.(ku1.eq.17)) call do_get_arval_funp2(knd,kp,5,2,1)
case (22,24:26) 
  do i=2,ar%n_eq_var
    if (ar%eq_var(i)) call do_get_arval_funp2(knd,kp,i,ar%var_ind(i),1)
  enddo
case (23) 
  avfp%v=a%arval_funp(1)%v
endselect
end

subroutine pg_get_psioml
!dec$ attributes dllexport:: pg_get_psioml
!нахождение функций в узлах области, построение сплайнов
use pgmod
integer(4) i0,j0
do i0=1,gs%na
  do j0=1,gs%a(i0)%nb
    if (gs%a(i0)%bnd(j0)%boundLineType.ne.1) cycle
    !call get_psioml_OLD(i0,j0,gs%a(i0)%bnd(j0)%npanel+1)
    call get_psioml(i0,j0)
  enddo
enddo
end

subroutine get_psioml(ia,knd)
use pgmod
integer(4) ia,knd
type(TArea), pointer :: a
type(TBound), pointer :: b
integer(4) di,i,nn,j,i1,k,i2,shift_inn,j1,j2,i_,ii1,ii2,n
real(8) epst,ds
complex(8) ett,et
integer(4), allocatable :: inn(:)
real(8), allocatable :: cc_(:,:),bb_(:),xx(:),yy(:),psioml(:,:),psioml2(:,:)
real(8), allocatable,target :: psiom_(:,:),sc_(:),s_(:)
real(8), pointer :: psiom(:,:),sc(:),s(:)
logical, allocatable :: is_nn1(:)
di=gs%const%di !сколько точек скраю не берем для производных, т.к. там содержатся большие прогрешности
a=>gs%a(ia)
b=>a%bnd(knd)
allocate(inn(b%npanel+1))
!ищем точки излома
if (gs%const%test_izlom) then
  epst=gs%const%angle_izlom*pi/180 !угол допустимого поворота, который не считается изломом
  epst=d2*dsin(epst*d5)
  ett=b%ett(b%npanel)
  nn=0
  do i=1,b%npanel
    et=b%ett(i)
    if (cdabs(et-ett)>epst) then
      nn=nn+1
      inn(nn)=i
    endif
    ett=et
  enddo
else
  inn(1)=1
  nn=1
endif
!строим сплайны
if (allocated(b%bb_psioml)) deallocate(b%bb_psioml)
if (allocated(b%cc_psioml)) deallocate(b%cc_psioml)
allocate(b%bb_psioml(b%npanel+1,a%umax))
allocate(b%cc_psioml(4,b%npanel+1,a%umax))
if (nn==0) then
  allocate(xx(b%npanel+1),yy(b%npanel+1))
  xx(1:b%npanel)=b%sc(1:b%npanel)
  ds=b%s(b%npanel+1)-b%s(1)
  xx(b%npanel+1)=b%sc(1)+ds
  do j=1,a%umax
    yy(1:b%npanel)=b%psiom(1:b%npanel,j)
    yy(b%npanel+1)=b%psiom(1,j)
    call dcsper(b%npanel+1,xx,yy,b%bb_psioml(:,j),b%cc_psioml(:,:,j))
  enddo
else
  if (inn(1)/=1) then
    allocate(sc_(b%npanel),s_(b%npanel+1),psiom_(b%npanel,a%umax))
    i1=inn(1)
    i2=b%npanel-i1+1
    ds=b%s(b%npanel+1)-b%s(1)
    forall (j=i1:b%npanel) sc_(j-i1+1)=b%sc(j)
    forall (j=1:i1-1) sc_(j+i2)=b%sc(j)+ds
    forall (j=i1:b%npanel) s_(j-i1+1)=b%s(j)
    forall (j=1:i1) s_(j+i2)=b%s(j)+ds
    forall (j=i1:b%npanel) psiom_(j-i1+1,:)=b%psiom(j,:)
    forall (j=1:i1-1) psiom_(j+i2,:)=b%psiom(j,:)
    psiom=>psiom_
    sc=>sc_
    s=>s_
    shift_inn=inn(1)-1
    inn(1:nn)=inn(1:nn)-shift_inn
  else
    psiom=>b%psiom
    sc=>b%sc
    s=>b%s
    shift_inn=0
  endif
  allocate(psioml(b%npanel,a%umax))
  allocate(psioml2(nn,a%umax))
  allocate(bb_(b%npanel+1),cc_(4,b%npanel+1))
  allocate(is_nn1(nn))
  do i=1,nn
    !индекс первой панели гладкого участка
    ii1=inn(i) 
    !индекс последней панели гладкого участка
    i_=i+1
    if (i_>nn) then
      ii2=b%npanel
    else
      ii2=inn(i_)-1
    endif
    !вычисляем значения функций на концах панелей
    is_nn1(i)=ii1==ii2
    if (ii1==ii2) then
      !для длинных панелей, занимающих целый участок, где полагаются все функции постоянными
      psioml(ii1,:)=psiom(ii1,:)
      psioml2(i,:)=psiom(ii1,:)
    else
      do j2=0,1
        i1=ii1
	      i2=ii2
	      if (j2==1) then
	        i1=ii1+di
	        i2=ii2-di
        endif
	      n=i2-i1+1
        do j1=1,a%umax,2
          j=j1+j2
          call dcsakm(n,sc(i1:i2),psiom(i1:i2,j),bb_,cc_)
          do k=ii1,ii2
	          psioml(k,j)=dcsval(s(k),n-1,bb_,cc_)
          enddo
          psioml2(i,j)=dcsval(s(ii2+1),n-1,bb_,cc_)
        enddo
      enddo
    endif
  enddo
  !предполагаем что функции непрерывны, а производные - нет
  if (.true.) then
    do i=1,nn
      !индекс первой точки гладкого участка
      ii1=inn(i) 
      !индекс предыдущего гладкого участка
      i_=i-1
      if (i_==0) i_=nn
      do j=1,a%umax,2
        if (is_nn1(i)==is_nn1(i_)) then
          psioml(ii1,j)=(psioml(ii1,j)+psioml2(i_,j))*d5
          psioml2(i_,j)=psioml(ii1,j)
        elseif (is_nn1(i)) then
          psioml2(i_,j)=psioml(ii1,j)
        else
          psioml(ii1,j)=psioml2(i_,j)
        endif
      enddo
    enddo
  endif
  !вычисляем коэффициенты сплайнов
  allocate(yy(b%npanel+1))
  do i=1,nn
    !индекс первой панели гладкого участка
    i1=inn(i) 
    !индекс последней панели гладкого участка
    i_=i+1
    if (i_>nn) then
      i2=b%npanel
    else
      i2=inn(i_)-1
    endif
    n=i2-i1+2
    do j=1,a%umax
      forall (k=i1:i2) yy(k-i1+1)=psioml(k,j)
      yy(n)=psioml2(i,j)
      if (is_nn1(i)) then
        call spline_const(n,s(i1:i2+1),yy,bb_,cc_)
      else
        call dcsakm(n,s(i1:i2+1),yy,bb_,cc_)
      endif
      forall (k=i1:i2+1) b%bb_psioml(k,j)=bb_(k-i1+1)
      forall (k=i1:i2+1) b%cc_psioml(:,k,j)=cc_(:,k-i1+1)
    enddo
  enddo
endif
deallocate(inn)
if (allocated(s_)) deallocate(s_)
if (allocated(sc_)) deallocate(sc_)
if (allocated(psiom_)) deallocate(psiom_)
if (allocated(psioml)) deallocate(psioml)
if (allocated(psioml2)) deallocate(psioml2)
if (allocated(xx)) deallocate(xx)
if (allocated(yy)) deallocate(yy)
if (allocated(bb_)) deallocate(bb_)
if (allocated(cc_)) deallocate(cc_)
if (allocated(is_nn1)) deallocate(is_nn1)
end

subroutine get_psioml_old(ia,knd,nmax)
use pgmod
integer(4) nmax
integer(4) i,bndlt,i1,i2,nn,j,k,jmax,ir(nmax),nir,knd,di,ia,i11,i22
logical is_nn1(nmax)
real(8) cc_(4,nmax),bb_(nmax),yy(nmax)
complex(8) et,ett
real(8) epst
!real(8) qqq(nmax,4)
type(TArea), pointer :: a
type(TBound), pointer :: b
real(8), allocatable :: psioml(:,:,:) !(npanel,umax,2) значения искомых функций на концах панелей (3-й индекс =1,2 - начало и конец панели) 
integer(4), allocatable :: bndl(:)
!real(8) fff(10001),sss(10001),pg_f_psioml
di=gs%const%di !сколько точек скраю не берем для производных, т.к. там содержатся большие прогрешности
a=>gs%a(ia)
b=>a%bnd(knd)
if (allocated(b%bb_psioml)) deallocate(b%bb_psioml)
if (allocated(b%cc_psioml)) deallocate(b%cc_psioml)
allocate(psioml(b%npanel,a%umax,2))
allocate(b%bb_psioml(b%npanel+1,a%umax))
allocate(b%cc_psioml(4,b%npanel+1,a%umax))
allocate(bndl(b%npanel))
k=0
do i=1,b%nline
  do j=1,b%line(i)%npanel
    k=k+1
	  bndl(k)=i
  enddo
enddo
psioml=d0
b%bb_psioml=d0
b%cc_psioml=d0
bndlt=bndl(1)
ett=b%ett(1)

epst=gs%const%angle_izlom*pi/180 !угол допустимого поворота, который не считается изломом
epst=d2*dsin(epst*d5)

i11=1
nir=1
ir(1)=1
do i=1,b%npanel
  et=b%ett(i)
  !if ((bndlt.ne.bndl(i)).or.(gs%const%test_izlom.and.cdabs(et-ett)>epst)) then
  !переход на новый участок может не быть изломом
  if (gs%const%test_izlom.and.cdabs(et-ett)>epst) then
    i22=i-1
  elseif (i==b%npanel) then
    i22=i
  else
    i22=0
  endif
  if (i22>0) then
    do j=1,a%umax
      if (i11==i22) then
        !для длинных панелей, занимающих целый участок, где полагаются все функции постоянными
        psioml(i11,j,:)=b%psiom(i11,j)
      else
	      i1=i11
	      i2=i22
	      if (mod(j,2)==0) then
	        i1=i11+di
	        i2=i22-di
        endif
	      nn=i2-i1+1
        call dcsakm(nn,b%sc(i1:i2),b%psiom(i1:i2,j),bb_,cc_)
	      do k=i11,i22
	        psioml(k,j,1)=dcsval(b%s(k),nn-1,bb_,cc_)
	      enddo
	      do k=i11,i22-1
	        psioml(k,j,2)=psioml(k+1,j,1)
	      enddo
	      psioml(i22,j,2)=dcsval(b%s(i22+1),nn-1,bb_,cc_)
      endif
    enddo
    is_nn1(nir)=i11==i22
	  if (i22<b%npanel) then
  	  i11=i22+1 
      bndlt=bndl(i11)
	    nir=nir+1
      ir(nir)=i11
    endif
  endif
  ett=et
enddo
!предполагаем что функции непрерывны, а производные - нет
if (.true.) then
  jmax=1
  if (a%umax>2) jmax=3
  do i=1,jmax,2
    do j=1,nir
	    i1=ir(j)
	    i2=ir(j)-1
      k=j-1
	    if (i2==0) then
        i2=b%npanel
        k=nir
      endif
      if (is_nn1(j)==is_nn1(k)) then
	      psioml(i1,i,1)=(psioml(i1,i,1)+psioml(i2,i,2))*d5
	      psioml(i2,i,2)=psioml(i1,i,1)
      elseif (is_nn1(j)) then
	      psioml(i2,i,2)=psioml(i1,i,1)
      else
        psioml(i1,i,1)=psioml(i2,i,2)
      endif
	  enddo
  enddo
endif
!qqq(:,:)=psioml(:,:,1,i0)
!qqq(:,:)=psioml(:,:,2,i0)
do i=1,a%umax
  do j=1,nir
	  i1=ir(j)
	  if (j==nir) then
	    i2=b%npanel+1
	  else
	    i2=ir(j+1)
    endif
    !if (mod(i,2)==0) then
	  !  i1=i1+di
	  !  i2=i2-di
    !endif
	  yy=d0
	  yy(i1:i2-1)=psioml(i1:i2-1,i,1)
	  yy(i2)=psioml(i2-1,i,2)
    nn=i2-i1+1
    if (is_nn1(j)) then
      call spline_const(nn,b%s(i1:i2),yy(i1:i2),bb_,cc_)
    else
	    !call dcsakm(i2-i1+1,b%s(i1:i2),yy(i1:i2),b%bb_psioml(i1:i2,i),b%cc_psioml(:,i1:i2,i)) непонятная ошибка при таком вызове
      call dcsakm(nn,b%s(i1:i2),yy(i1:i2),bb_,cc_)
    endif
    b%bb_psioml(i1:i2,i)=bb_(1:nn)
    b%cc_psioml(:,i1:i2,i)=cc_(:,1:nn)
  enddo
  !do j=1,10001
  !  sss(j)=(j-d1)/10000*s1(npanel(i0)+1,i0)
  !  fff(j)=pg_f_psioml(sss(j),i,i0)
  !enddo
  !nir=nir
enddo 
deallocate(psioml)
end

recursive function pg_f_psioml(ibnd,s,nf,nder)
!dec$ attributes dllexport:: pg_f_psioml
!вычисление функции и производной на границе по дуговой абсциссе
use pgmod 
real(8) s  !дуговая абсцисса
integer(4) nf !номер функции
integer(4) nder !порядок производной
integer(4) ibnd !номер текущей границы
real(8) pg_f_psioml,g,ds,s0,s1
type(TBound), pointer :: b
type(TBoundline2), pointer :: bl2
type(TBounLineFuncApprox), pointer :: ga
b=>gsarea%bnd(ibnd)
pg_f_psioml=d0
if (gsarea%nu==1.and.nf>2) return
if (b%boundLineType==2) then
  bl2=>b%line2(1)  !!!тут только одна линия
  ga=>bl2%ga(nf)
  if (ga%typea==1) then
    g=s/bl2%r
    if (nder==0) then
      pg_f_psioml=ga%valkk(1)+ga%valkk(2)*dsin(g)+ga%valkk(3)*dcos(g)
	  elseif (nder==1) then
      pg_f_psioml=(ga%valkk(2)*dcos(g)-ga%valkk(3)*dsin(g))/bl2%r
	  endif
  endif 
else
  s0=b%bb_psioml(1,nf)
  s1=b%bb_psioml(b%npanel+1,nf)
  if (s<s0) then
    ds=s1-s0
    pg_f_psioml=pg_f_psioml(ibnd,s+ds,nf,nder)
  elseif (s>b%bb_psioml(b%npanel+1,nf)) then
    ds=s1-s0
    pg_f_psioml=pg_f_psioml(ibnd,s-ds,nf,nder)
  else
    pg_f_psioml=dcsder(nder,s,b%npanel,b%bb_psioml(:,nf),b%cc_psioml(:,:,nf))
  endif
endif
end

recursive function pg_int_psioml(ibnd,s0,s1,nf)
!dec$ attributes dllexport:: pg_int_psioml
!интеграл от функции и производной на границе по дуговой абсциссе
use pgmod 
real(8) s0,s1  !дуговая абсцисса начала и конца
integer(4) nf !номер функции
integer(4) ibnd !номер текущей границы
real(8), parameter :: eps=1.0d-8
real(8) pg_int_psioml,g0,g1,s00,s11,ds
type(TBound), pointer :: b
type(TBoundline2), pointer :: bl2
type(TBounLineFuncApprox), pointer :: ga
b=>gsarea%bnd(ibnd)
pg_int_psioml=d0
if (gsarea%nu==1.and.nf>2) return
if (b%boundLineType==2) then
  bl2=>b%line2(1)  !!!тут только одна линия
  ga=>bl2%ga(nf)
  if (ga%typea==1) then
    g0=s0/bl2%r
	  g1=s1/bl2%r
	  !!!не тестировал
	  pg_int_psioml=(ga%valkk(1)*(g1-g0)-ga%valkk(2)*(dcos(g1)-dcos(g0))+ga%valkk(3)*(dsin(g1)-dsin(g0)))*bl2%r
  endif
else
  s00=b%bb_psioml(1,nf)
  s11=b%bb_psioml(b%npanel+1,nf)
  if (s0<s00-eps) then
    ds=s11-s00
    pg_int_psioml=pg_int_psioml(ibnd,s00,s1,nf)+pg_int_psioml(ibnd,s0+ds,s11,nf)
  elseif (s1>s11+eps) then
    ds=s11-s00
    pg_int_psioml=pg_int_psioml(ibnd,s0,s11,nf)+pg_int_psioml(ibnd,s00,s1-ds,nf)
  else
    pg_int_psioml=dcsitg(s0,s1,b%npanel,b%bb_psioml(:,nf),b%cc_psioml(:,:,nf))
  endif
endif
end

function pg_f_psiomlxy(xx,yy,nf,nder)
!dec$ attributes dllexport:: pg_f_psiomlxy
!вычисление функции и производной на границе по координатам (x,y)
use pgmod 
integer(4) nf,nder,ibnd
real(8) pg_f_psioml,s,pg_f_psiomlxy,xx,yy
integer(4) in_bound,pg_get_bound_s
in_bound=pg_get_bound_s(xx,yy,ibnd,s)
if (in_bound==1) then
  pg_f_psiomlxy=pg_f_psioml(ibnd,s,nf,nder)
else
  pg_f_psiomlxy=d0
endif
end

function pg_get_bound_s_ind(xx,yy,ibnd,j,smin)
!dec$ attributes dllexport:: pg_get_bound_s_ind
!вычисление индекса участка границы по координатам (x,y)
!!!устаревшая, для совместимости со старыми программами
use pgmod 
integer(4) j  !возвращаемый номер панели, на которой лежит точка (x,y)
real(8) smin  !найденная локальная дуговая абсцисса, 0 - начало j-й панели, 1 - конец j-й панели
integer(4) ibnd
real(8) xx,yy
logical pg_get_bound_s_ind,pg_get_bound_s_ind_bt
pg_get_bound_s_ind=pg_get_bound_s_ind_bt(xx,yy,ibnd,j,smin,.false.)
end

function pg_get_bound_s_ind_bt(xx,yy,ibnd,j,smin,only_bt1)
!dec$ attributes dllexport:: pg_get_bound_s_ind_bt
!вычисление индекса участка границы по координатам (x,y)
!с условием на тип границы
!!!устаревшая, для совместимости со старыми программами
use pgmod 
integer(4) j  !возвращаемый номер панели, на которой лежит точка (x,y)
real(8) smin  !найденная локальная дуговая абсцисса, 0 - начало j-й панели, 1 - конец j-й панели
integer(4) i,ibnd,k
logical only_bt1 !искать только для boundLineType==1
logical pg_get_bound_s_ind_bt
real(8) s,xx,yy,minDist,dist,DistanceToLine !,px,py,DistanceToSegment
type(TBound), pointer :: b
type(TBoundline2), pointer :: bl2
complex(8) z
minDist=1.0d10
j=0
do k=1,gsarea%nb
  b=>gsarea%bnd(k)
  if (b%boundLineType==2) then
    if (only_bt1) cycle
    do i=1,b%nline
      bl2=>b%line2(i)  
	    z=dcmplx(xx,yy)-bl2%zc
      dist=dabs(cdabs(z)-bl2%r)
	    if (dist<minDist) then
        minDist=dist
	      j=i
	      smin=(zarg(z)-bl2%g0)*bl2%dir/pi2
		    if (smin<d0) smin=smin+d1
	      ibnd=k
      endif
	  enddo
  else
    dist=DistanceToLine(xx,yy,b%x,b%y,b%npanel,i,s)
    if (dist<minDist) then
      minDist=dist
	    j=i
	    smin=s
	    ibnd=k
    endif
  endif
enddo
!if (minDist>1.0d-2) then
!  call gs_print_stop("!!!!!point not in bound")
!endif
b=>gsarea%bnd(ibnd)
if (b%boundLineType==2) then
  bl2=>b%line2(j)  
  pg_get_bound_s_ind_bt=minDist/bl2%r<1.0d-4
else
  !pg_get_bound_s_ind_bt=minDist/b%l(j)<1.0d-2
  pg_get_bound_s_ind_bt=minDist<gs%const%epsMinDistBound
endif
end

function pg_get_bound_s_bt(xx,yy,ibnd,smin,ttmin,only_bt1,bn)
!dec$ attributes dllexport:: pg_get_bound_s_bt
!вычисление индекса участка границы и дуговой абсциссы по координатам (x,y)
!с условием на тип границы
use pgmod 
real(8) smin  !найденная дуговая абсцисса ближайшей точки границы 
real(8) ttmin !угол наклона касательной ближайшей точки на гарнице
integer(4) i,ibnd,k,j
logical only_bt1 !искать только для boundLineType==1
type(TBound_near) bn,bnt !информация о ближайшей точке на границе
integer(4) pg_get_bound_s_bt !0 - область, 1 - граница, 2 - окрестность границы (расстояние меньше длины панели)
real(8) gam0,gam,tt,s,xx,yy,minDist,dist,DistanceToBound,vect_p !,px,py,DistanceToSegment
type(TBound), pointer :: b
type(TBoundline2), pointer :: bl2
complex(8) z,ztt
integer(4) in_bound
minDist=1.0d10
j=0
do k=1,gsarea%nb
  b=>gsarea%bnd(k)
  if (b%boundLineType==2) then
    if (only_bt1) cycle
    do i=1,b%nline
      bl2=>b%line2(i)  
	    z=dcmplx(xx,yy)-bl2%zc
      dist=dabs(cdabs(z)-bl2%r)
	    if (dist<minDist) then
        minDist=dist
        j=i
        gam0=zarg(z)
	      gam=(gam0-bl2%g0)*bl2%dir
		    if (gam<d0) gam=gam+pi2
        smin=gam*bl2%r
        ttmin=gam0+pi5*bl2%dir
	      ibnd=k
      endif
	  enddo
  else
    dist=DistanceToBound(xx,yy,k,s,tt,bnt)
    if (dist<minDist) then
      minDist=dist
	    smin=s
      ttmin=tt
	    ibnd=k
      bn=bnt
    endif
  endif
  if (minDist<gs%const%epsMinDistBound2) exit
enddo
b=>gsarea%bnd(ibnd)
in_bound=0
if (b%boundLineType==2) then
  bl2=>b%line2(j)  
  if (minDist/bl2%r<gs%const%epsMinDistBound3) in_bound=1
else
  if (gs_test_point_near_bound) then
    if (minDist/bn%panelL<gs%const%epsMinDistBound3) then
      in_bound=1
    elseif (minDist<bn%panelL) then
      in_bound=2
    endif
  else
    if (minDist<gs%const%epsMinDistBound) in_bound=1
  endif
  if (in_bound/=1) then !определяем в области или снаружи точка
    ztt=cdexp(ii*ttmin)
    z=dcmplx(xx,yy)-dcmplx(bn%x,bn%y)
    bn%in_domain=vect_p(ztt,z)>d0
  endif
endif
pg_get_bound_s_bt=in_bound
end

function DistanceToBound(xx,yy,ib,smin,ttmin,bn)
!расстояние от точки (xx,yy) до границы (ib) текущей области
use pgmod
real(8) xx,yy !координат точки
integer(4) ib !номер границы
real(8) smin  !дуговая абсцисса ближайшей точки на гарнице
real(8) ttmin !угол наклона касательной ближайшей точки на гарнице
type(TBound_near) bn !информация о ближайшей точке на границе
real(8) DistanceToBound
real(8) DistanceToSegment,DistanceToCircleArc,sred,DistanceToLine,get_panelL
type(TBound), pointer :: b
type(TBoundline_geomdetale), pointer :: gd
integer(4) i,j,imin,jmin
real(8) px,py,dist,minDist,s,gam
b=>gsarea%bnd(ib)
minDist=1.0d10
imin=0
do i=1,b%ngeom_detale
  gd=>b%geom_detale(i)
  select case (gd%mode)
  case (1)
    dist=DistanceToLine(xx,yy,b%x(gd%i_begin:gd%i_end),b%y(gd%i_begin:gd%i_end),gd%i_end-gd%i_begin,j,s)
  case (2)
    dist=DistanceToSegment(gd%x,gd%y,gd%x2,gd%y2,xx,yy,s,px,py)
  case (3)
    dist=DistanceToCircleArc(xx,yy,gd%x,gd%y,gd%r,gd%gam1,gd%gam2,s)
  endselect
  if (dist<minDist) then
    minDist=dist
	  smin=s
    if (gd%mode==1) jmin=j
    imin=i
  endif
  if (minDist<gs%const%epsMinDistBound2) exit
enddo
gd=>b%geom_detale(imin)
select case (gd%mode)
case (1)
  j=jmin+gd%i_begin-1
  bn%x=sred(d0,d1,b%x(j),b%x(j+1),smin)
  bn%y=sred(d0,d1,b%y(j),b%y(j+1),smin)
  smin=sred(d0,d1,b%s(j),b%s(j+1),smin)
  ttmin=zarg(b%ett(j))
  if (gs_test_point_near_bound) bn%panelL=b%l(j)
case (2)
  bn%x=sred(d0,d1,gd%x,gd%x2,smin)
  bn%y=sred(d0,d1,gd%y,gd%y2,smin)
  smin=sred(d0,d1,gd%s1,gd%s2,smin)
  ttmin=zarg(b%ett(gd%i_begin))
  if (gs_test_point_near_bound) bn%panelL=get_panelL(b,gd,smin,minDist)
case (3)
  gam=smin
  smin=sred(gd%gam1,gd%gam2,gd%s1,gd%s2,gam)
  bn%x=gd%x+gd%r*dcos(gam)
  bn%y=gd%y+gd%r*dsin(gam)
  if (gd%gam2>gd%gam1) then
    ttmin=gam+pi5
  else
    ttmin=gam-pi5
  endif
  if (gs_test_point_near_bound) bn%panelL=get_panelL(b,gd,smin,minDist)
endselect
bn%dist=minDist
DistanceToBound=minDist
end

function get_panelL(b,gd,s,Dist)
!найти длину ближайшей панели
use pgmod
real(8) get_panelL
type(TBound) b
type(TBoundline_geomdetale) gd
real(8) s !дуговая абсцисса
real(8) Dist !расстояние до границы (для проверки)
integer(4) i,j
if (Dist>gd%panelLmax) then
  get_panelL=gd%panelLmax
else
  j=-1
  do i=gd%i_begin+1,gd%i_end
    if (s<b%s(i)) then
      j=i
      exit
    endif
  enddo
  if (j<0) then
    j=gd%i_end-1
  else
    j=j-1
  endif
  get_panelL=b%l(j)
endif
end

function pg_get_bound_s(xx,yy,ibnd,s)
!dec$ attributes dllexport:: pg_get_bound_s
!вычисление дуговой абсциссы по координатам (x,y)
use pgmod 
real(8) s  !дуговая абсцисса
integer(4) ibnd
integer(4) pg_get_bound_s,pg_get_bound_s_bt
real(8) xx,yy,tt
type(TBound_near) bn
pg_get_bound_s=pg_get_bound_s_bt(xx,yy,ibnd,s,tt,.false.,bn)
end

function pg_get_bound_s_tt(xx,yy,ibnd,s,tt)
!dec$ attributes dllexport:: pg_get_bound_s_tt
!вычисление дуговой абсциссы и угла наклона касательной к панели по координатам (x,y)
use pgmod 
real(8) s  !дуговая абсцисса
real(8) tt  !угол наклона касательной
integer(4) ibnd
integer(4) pg_get_bound_s_tt,pg_get_bound_s_bt
real(8) xx,yy
type(TBound_near) bn
pg_get_bound_s_tt=pg_get_bound_s_bt(xx,yy,ibnd,s,tt,.false.,bn)
end

function pg_get_fun_xy2(xx,yy,nf,xder,yder,mode,in_bound)
!dec$ attributes dllexport:: pg_get_fun_xy2
!вычисление функции и производной в области или на границе по координатам (x,y)
!nf - 1 - глаавная функция (psi)
     !2 - производная главной функции по направлению (xder,yder)
	   !3 - вторая функция (om)
	   !4 - производная второй функции по направлению (xder,yder)
!mode - 0 - определить, область или граница
       !1 - область
	     !2 - граница
       !3 - определить, область или граница (только для boundLineType=1)
       !4 - граница (только для boundLineType=1)
use pgmod 
integer(4) nf,ibnd,mode
real(8) s,pg_get_fun_xy2,pg_get_fun_xy_solve,xx,yy,tt,xder,yder
type(TBound_near) bn
integer(4) in_bound
call pg_get_fun_xy_prepare(xx,yy,mode,in_bound,ibnd,s,tt,bn)
pg_get_fun_xy2=pg_get_fun_xy_solve(xx,yy,nf,xder,yder,in_bound,ibnd,s,tt,bn)
end

subroutine pg_get_fun_xy_prepare(xx,yy,mode,in_bound,ibnd,s,tt,bn)
!dec$ attributes dllexport:: pg_get_fun_xy_prepare
!подготовка к вычислению функции и производной в области или на границе по координатам (x,y)
!mode - 0 - определить, область или граница
       !1 - область
	     !2 - граница
       !3 - определить, область или граница (только для boundLineType=1)
       !4 - граница (только для boundLineType=1)
use pgmod 
integer(4) ibnd,mode
real(8) s,xx,yy,tt
type(TBound_near) bn
integer(4) in_bound,pg_get_bound_s_bt
in_bound=0
bn%in_domain=.true.
selectcase (mode)
case (0,2)
  in_bound=pg_get_bound_s_bt(xx,yy,ibnd,s,tt,.false.,bn)
case (3,4)
  in_bound=pg_get_bound_s_bt(xx,yy,ibnd,s,tt,.true.,bn)
endselect
end

subroutine get_fun_xy_prepare_bi(xx,yy,mode,bi)
!подготовка к вычислению функции и производной в области или на границе по координатам (x,y)
!mode - 0 - определить, область или граница
       !1 - область
	     !2 - граница
       !3 - определить, область или граница (только для boundLineType=1)
       !4 - граница (только для boundLineType=1)
use pgmod 
integer(4) mode
real(8) xx,yy
type(TBound_info) bi
call pg_get_fun_xy_prepare(xx,yy,mode,bi%in_bound,bi%ibnd,bi%s,bi%tt,bi%bn)
end

function pg_get_fun_xy_solve(xx,yy,nf,xder,yder,in_bound,ibnd,s,tt,bn)
!dec$ attributes dllexport:: pg_get_fun_xy_solve
!вычисление функции и производной в области или на границе по координатам (x,y)
!nf - 1 - глаавная функция (psi)
     !2 - производная главной функции по направлению (xder,yder)
	   !3 - вторая функция (om)
	   !4 - производная второй функции по направлению (xder,yder)
use pgmod 
integer(4) nf,ibnd
real(8) s,pg_get_fun_xy_solve,xx,yy,tt,xder,yder,ff(2)
type(TBound_near) bn
integer(4) in_bound
call pg_get_fun_xy_solve_gradient(xx,yy,nf,xder,yder,in_bound,ibnd,s,tt,ff,bn)
pg_get_fun_xy_solve=ff(1)
end

function get_fun_xy_solve_bi(xx,yy,nf,xder,yder,bi)
!вычисление функции и производной в области или на границе по координатам (x,y)
!nf - 1 - глаавная функция (psi)
     !2 - производная главной функции по направлению (xder,yder)
	   !3 - вторая функция (om)
	   !4 - производная второй функции по направлению (xder,yder)
use pgmod 
integer(4) nf
real(8) get_fun_xy_solve_bi,pg_get_fun_xy_solve,xx,yy,xder,yder
type(TBound_info) bi
get_fun_xy_solve_bi=pg_get_fun_xy_solve(xx,yy,nf,xder,yder,bi%in_bound,bi%ibnd,bi%s,bi%tt,bi%bn)
end

recursive subroutine pg_get_fun_xy_solve_gradient(xx,yy,nf,xder,yder,in_bound,ibnd,s,tt,ff,bn)
!dec$ attributes dllexport:: pg_get_fun_xy_solve_gradient
!вычисление функции и производной в области или на границе по координатам (x,y)
!nf - 1 - глаавная функция (psi)
     !2 - производная главной функции по направлению (xder,yder)
	   !3 - вторая функция (om)
	   !4 - производная второй функции по направлению (xder,yder)
     !5 - градиент главной функции
     !6 - градиент второй функции
use pgmod 
integer(4) nf,ibnd,k,i
real(8) pg_f_psioml,s,xx,yy,tt,xder,yder,eps,der1,der2,ff(2),val1,val2,ff1(2),xx1,yy1,tt1,sred,tt2,ll
type(TBound_near) bn
complex(8) lder,tau,nn,gg,zz1,zz2,csred
real(8) pg_get_fun,pg_get_fun_1,pg_get_fun2,pg_get_fun2_1,scal_p
integer(4) in_bound
logical q1,q2
type(TBound), pointer :: b
eps=1.0d-8
ff=d0
if (nf==2.or.nf==4) then
  lder=dcmplx(xder,yder)
  lder=lder/cdabs(lder)
endif
select case (in_bound)
case (1)
  !на границе
  b=>gsarea%bnd(ibnd)
  if (nf==1.or.nf==3) then
    ff(1)=pg_f_psioml(ibnd,s,nf,0)
  else
    tau=cdexp(ii*tt)
    nn=-ii*tau    !внешняя нормаль
    if(nf==2.or.nf==4) then  
      der1=scal_p(tau,lder) !скалярное произведение tau*lder
      der2=scal_p(nn,lder) !скалярное произведение nn*lder
      if (dabs(der1)>eps) ff(1)=ff(1)+pg_f_psioml(ibnd,s,nf-1,1)*der1
	    if (dabs(der2)>eps) ff(1)=ff(1)+pg_f_psioml(ibnd,s,nf,0)*der2
    else
      k=3
      if (nf==6) k=2
      gg=dcmplx(pg_f_psioml(ibnd,s,nf-k,0),pg_f_psioml(ibnd,s,nf-k-1,1))
      gg=gg*nn
      ff(1)=dreal(gg)
      ff(2)=dimag(gg)
    endif
  endif
case (0)
  !в области
  if (bn%in_domain) then
    if (nf==2.or.nf==4) then
      der1=dreal(lder)
      der2=dimag(lder)
    endif
    select case (nf)
    case (1)
      ff(1)=pg_get_fun(xx,yy)
    case (2)
      q1=dabs(der1)>eps
      q2=dabs(der2)>eps
      if (gs_use_dual_integral_solve.and.q1.and.q2) then
        call pg_get_fun_1_dual(xx,yy,val1,val2)
      else
        if (q1) val1=pg_get_fun_1(xx,yy,1)
	      if (q2) val2=pg_get_fun_1(xx,yy,2)
      endif
      if (q1) ff(1)=ff(1)+val1*der1
      if (q2) ff(1)=ff(1)+val2*der2
    case (3)
      ff(1)=pg_get_fun2(xx,yy)
    case (4)
      if (dabs(der1)>eps) ff(1)=ff(1)+pg_get_fun2_1(xx,yy,1)*der1
	    if (dabs(der2)>eps) ff(1)=ff(1)+pg_get_fun2_1(xx,yy,2)*der2
    case (5)
      if (gs_use_dual_integral_solve) then
        call pg_get_fun_1_dual(xx,yy,ff(1),ff(2))
      else
        ff(1)=pg_get_fun_1(xx,yy,1)
        ff(2)=pg_get_fun_1(xx,yy,2)
      endif
    case (6)
      if (gs_use_dual_integral_solve) then
        call pg_get_fun2_1_dual(xx,yy,ff(1),ff(2))
      else
        ff(1)=pg_get_fun2_1(xx,yy,1)
        ff(2)=pg_get_fun2_1(xx,yy,2)
      endif
    end select
  else
    call pg_get_fun_xy_solve_gradient(bn%x,bn%y,nf,xder,yder,1,ibnd,s,tt,ff,bn) !на границе
  endif
case (2)
  !в окрестности границы
  call pg_get_fun_xy_solve_gradient(bn%x,bn%y,nf,xder,yder,1,ibnd,s,tt,ff,bn) !на границе
  if (bn%in_domain) then
    tt1=tt+pi5 !направление нормали внутрь области
    ll=bn%panelL*0.2d0
    zz1=cdexp(ii*tt1)
    if (bn%dist>ll) then !подправляем угол, если направления нормалей tt1 и tt2 не совпадают
      tt2=datan2(yy-bn%y,xx-bn%x)
      zz2=cdexp(ii*tt2)
      if (scal_p(zz1,zz2)<d0) zz2=-zz2
      zz1=csred(ll,bn%panelL,zz1,zz2,bn%dist)
    endif
    xx1=bn%x+bn%panelL*dreal(zz1)
    yy1=bn%y+bn%panelL*dimag(zz1)
    call pg_get_fun_xy_solve_gradient(xx1,yy1,nf,xder,yder,0,ibnd,s,tt,ff1,bn) !в области
    do i=1,2
      ff(i)=sred(d0,bn%panelL,ff(i),ff1(i),bn%dist)
      if (nf<5) exit
    enddo
  endif
end select
end

subroutine get_fun_xy_solve_gradient_bi(xx,yy,nf,xder,yder,bi,ff)
!вычисление функции и производной в области или на границе по координатам (x,y)
!nf - 1 - глаавная функция (psi)
     !2 - производная главной функции по направлению (xder,yder)
	   !3 - вторая функция (om)
	   !4 - производная второй функции по направлению (xder,yder)
     !5 - градиент главной функции
     !6 - градиент второй функции
use pgmod 
integer(4) nf
real(8) xx,yy,xder,yder,ff(2)
type(TBound_info) bi
call pg_get_fun_xy_solve_gradient(xx,yy,nf,xder,yder,bi%in_bound,bi%ibnd,bi%s,bi%tt,ff,bi%bn)
end

function pg_get_funf_xy(xx,yy,nf,mode)
!dec$ attributes dllexport:: pg_get_funf_xy
use pgmod
integer(4) nf,mode
real(8) pg_get_funf_xy,pg_get_fun_xy,xx,yy
pg_get_funf_xy=pg_get_fun_xy(xx,yy,nf,d0,d0,mode)
end

function pg_get_fun_xy(xx,yy,nf,xder,yder,mode)
!dec$ attributes dllexport:: pg_get_fun_xy
integer(4) nf,mode
real(8) pg_get_fun_xy2,pg_get_fun_xy,xx,yy,xder,yder
integer(4) in_bound
pg_get_fun_xy=pg_get_fun_xy2(xx,yy,nf,xder,yder,mode,in_bound)
end

subroutine pg_get_fun_xy_gradient(xx,yy,nf,mode,ff)
!dec$ attributes dllexport:: pg_get_fun_xy_gradient
use pgmod
integer(4) nf,mode,ibnd
real(8) xx,yy,xder,yder,s,tt
real(8) ff(2) !градиент функции (nf=5 - главной функции, nf=6 - второй функции)
integer(4) in_bound
type(TBound_near) bn
call pg_get_fun_xy_prepare(xx,yy,mode,in_bound,ibnd,s,tt,bn)
call pg_get_fun_xy_solve_gradient(xx,yy,nf,xder,yder,in_bound,ibnd,s,tt,ff,bn)
end
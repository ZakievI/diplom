function intd_z(z0,ia,j,nf)
!интеграл, когда контрольная точка является серединой панели
use pgmod
!j - номер тек треугольника, по которому идет интегрирование
!nf - номер интеграла
!ia - номера области
integer(4) j,nf,ia
real(8) intd_z,intd_an
complex(8) zmz,z0
if (gs%a(ia)%a%npe_ar(j)==3) then
  intd_z=intd_an(z0,zmz(1,j,ia),zmz(2,j,ia),zmz(3,j,ia),nf)
else
  intd_z=intd_an(z0,zmz(1,j,ia),zmz(2,j,ia),zmz(3,j,ia),nf)+intd_an(z0,zmz(1,j,ia),zmz(3,j,ia),zmz(4,j,ia),nf)
endif
end

subroutine intd_z_2(z0,ia,j,nf,val,val2)
!интеграл, когда контрольная точка является серединой панели
use pgmod
!j - номер тек треугольника, по которому идет интегрирование
!nf - номер интеграла
!ia - номера области
integer(4) j,nf,ia
real(8) val,val2,val3,val4
complex(8) zmz,z0
if (gs%a(ia)%a%npe_ar(j)==3) then
  call intd_an_2(z0,zmz(1,j,ia),zmz(2,j,ia),zmz(3,j,ia),nf,val,val2,.true.)
else
  call intd_an_2(z0,zmz(1,j,ia),zmz(2,j,ia),zmz(3,j,ia),nf,val,val2,.true.)
  call intd_an_2(z0,zmz(1,j,ia),zmz(3,j,ia),zmz(4,j,ia),nf,val3,val4,.true.)
  val=val+val3
  val2=val2+val4
endif
end
  
function intd_m(ia,j,i,knd,nf)
!интеграл, когда контрольная точка является серединой панели
use pgmod
!j - номер тек треугольника, по которому идет интегрирование
!i - номер контрольной точки 
!nf - номер интеграла
!knd - номер границы с контрольной точкой
!ia - номера области
integer(4) j,i,nf,knd,ia
real(8) intd_m,intd_z,get_ab_cash_integral,val
type(TBound), pointer :: b
val=real8_inf
if (gs_use_cash) then
  if (gs_cash_presolve) then
    intd_m=get_ab_cash_integral_2(gs%a(ia),j,knd,i,nf)
    !DEC$ IF DEFINED (DEBUG)
    if (intd_m==real8_inf) call gs_print_stop('Error intd_m gs_cash_presolve!!!')
    !DEC$ ENDIF
    return
  endif
  val=get_ab_cash_integral(ia,j,knd,i,nf)
endif
if (val==real8_inf) then
  b=>gs%a(ia)%bnd(knd)
  val=intd_z(b%zc(i),ia,j,nf)
  if (gs_use_cash) call set_ab_cash_integral(ia,j,knd,i,nf,val)
endif
intd_m=val
end

function intd_ar(ia,j,i,nf)
!интеграл, когда контрольная точка является центром треугольника
use pgmod
!j - номер тек треугольника, по которому идет интегрирование
!i - номер треугольника с контрольной точки 
!nf - номер интеграла
!knd - номер границы с контрольной точкой
!ia - номера области
integer(4) j,i,nf,ia
real(8) intd_ar,intd_z,get_aa_cash_integral,val
type(areatype), pointer :: aa
val=real8_inf
if (gs_use_cash) then
  if (gs_cash_presolve) then
    intd_ar=get_aa_cash_integral_2(gs%a(ia),j,i,nf)
    !DEC$ IF DEFINED (DEBUG)
    if (intd_ar==real8_inf) call gs_print_stop('Error intd_ar gs_cash_presolve!!!')
    !DEC$ ENDIF
    return
  endif
  val=get_aa_cash_integral(ia,j,i,nf)
endif
if (val==real8_inf) then
  aa=>gs%a(ia)%a
  val=intd_z(aa%zmc(i),ia,j,nf)
  if (gs_use_cash) call set_aa_cash_integral(ia,j,i,nf,val)
endif
intd_ar=val
end

function intd_num(z0,z1,z2,z3,nf,isOss)
!численный интеграл от G_1 по треугольнику (z1,z2,z3)
!nf=0 - G_1
   !2 - G_2
   !4 - dG_1/dx
   !5 - dG_1/dy
   !8 - dG_2/dx
   !9 - dG_2/dy
!isOss=0 - G_1
      !1 - G_1/y
use gen_mod
real(8) intd_num,xmin,xmax,ymin,ymax,ds,dsx,dsy,xt,yt,s,r,ff,ff2
real(8) area1,area2
complex(8) z,z1,z2,z3,z0,dz
integer(4) nx,ny,i,j,n,nf,isOss
logical test_point_in_tr
ds=1.0d-3
xmin=min(dreal(z1),dreal(z2),dreal(z3))
xmax=max(dreal(z1),dreal(z2),dreal(z3))
ymin=min(dimag(z1),dimag(z2),dimag(z3))
ymax=max(dimag(z1),dimag(z2),dimag(z3))
nx=(xmax-xmin)/ds
ny=(ymax-ymin)/ds
dsx=(xmax-xmin)/nx
dsy=(ymax-ymin)/ny
s=d0
n=0
do i=1,nx
  xt=xmin+(i-d5)*dsx
  do j=1,ny
    yt=ymin+(j-d5)*dsy
    z=dcmplx(xt,yt)
	  if (test_point_in_tr(z,z1,z2,z3)) then
	    n=n+1
	    dz=z-z0
	    r=cdabs(dz)
	    if (r>1.0d-8) then
	      ff=d1
	      if (isOss==1) ff=d1/yt
        select case (nf)
		    case (0)
		      ff2=dlog(r)
        case (2)
          ff2=r*r*(dlog(r)-d1)*0.25d0
		    case (4)
		      ff2=-dreal(dz)/(r**2)
		    case (5)
		      ff2=-dimag(dz)/(r**2)
        case (8)
		      ff2=-dreal(dz)*(d2*dlog(r)-d1)*0.25d0
		    case (9)
		      ff2=-dimag(dz)*(d2*dlog(r)-d1)*0.25d0
        case default
          call gs_print_stop("Error intd_num!")
		    end select
	      s=s+ff2*ff
	    endif
	  endif
  enddo
enddo
area1=n*dsx*dsy
area2=d5*dabs(dimag(dconjg(z2-z1)*(z3-z1)))
intd_num=s*dsx*dsy
end

function intd_oss_m(ia,j,i,knd,nf)
!интеграл, когда контрольная точка является серединой панели
use pgmod
!j - номер тек треугольника, по которому идет интегрирование
!i - номер контрольной точки 
!nf - номер интеграла
!knd - номер границы с контрольной точкой
!ia - номера области
integer(4) j,i,nf,knd,ia
real(8) intd_oss_z,intd_oss_m,val,get_ab_cash_integral
type(TBound), pointer :: b
val=real8_inf
if (gs_use_cash) then
  if (gs_cash_presolve) then
    intd_oss_m=get_ab_cash_integral_2(gs%a(ia),j,knd,i,nf)
    !DEC$ IF DEFINED (DEBUG)
    if (intd_oss_m==real8_inf) call gs_print_stop('Error intd_oss_m gs_cash_presolve!!!')
    !DEC$ ENDIF
    return
  endif
  val=get_ab_cash_integral(ia,j,knd,i,nf)
endif
if (val==real8_inf) then
  b=>gs%a(ia)%bnd(knd)
  val=intd_oss_z(b%zc(i),ia,j,nf)
  if (gs_use_cash) call set_ab_cash_integral(ia,j,knd,i,nf,val)
endif
intd_oss_m=val
end

function intd_oss_ar(ia,j,i,nf)
!интеграл, когда контрольная точка является центром треугольника
use pgmod
!j - номер тек треугольника, по которому идет интегрирование
!i - номер треугольника с контрольной точкой
!nf - номер интеграла
!ia - номера области
integer(4) j,i,nf,ia
real(8) intd_oss_z,intd_oss_ar,val,get_aa_cash_integral
type(areatype), pointer :: aa
val=real8_inf
if (gs_use_cash) then
  if (gs_cash_presolve) then
    intd_oss_ar=get_aa_cash_integral_2(gs%a(ia),j,i,nf)
    !DEC$ IF DEFINED (DEBUG)
    if (intd_oss_ar==real8_inf) call gs_print_stop('Error intd_oss_ar gs_cash_presolve!!!')
    !DEC$ ENDIF
    return
  endif
  val=get_aa_cash_integral(ia,j,i,nf)
endif
if (val==real8_inf) then
  aa=>gs%a(ia)%a
  val=intd_oss_z(aa%zmc(i),ia,j,nf)
  if (gs_use_cash) call set_aa_cash_integral(ia,j,i,nf,val)
endif
intd_oss_ar=val
end

function intd_oss_z(z0,ia,j,nf)
!интеграл, когда контрольная точка является центром треугольника
use pgmod
!j - номер тек треугольника, по которому идет интегрирование
!nf - номер интеграла
!knd - номер границы с контрольной точкой
!ia - номера области
integer(4) j,nf,ia
real(8) intd_oss,intd_oss_z
complex(8) zmz,z0
if (gs%a(ia)%a%npe_ar(j)==3) then
  intd_oss_z=intd_oss(z0,zmz(1,j,ia),zmz(2,j,ia),zmz(3,j,ia),nf)
else
  intd_oss_z=intd_oss(z0,zmz(1,j,ia),zmz(2,j,ia),zmz(3,j,ia),nf)+intd_oss(z0,zmz(1,j,ia),zmz(3,j,ia),zmz(4,j,ia),nf)
endif
end

function intd_oss(z0,z1,z2,z3,nf)
!вычисление интеграла G_1/y в осесимметричном случае
complex(8) z1,z2,z3,z0
real(8) intd_num,intd_oss_,intd_oss
integer(4) use_num,nf
use_num=0
if (use_num==1) then
  intd_oss=intd_num(z0,z1,z2,z3,nf,1)
else
  intd_oss=intd_oss_(z0,z1,z2,z3,nf)
endif
end

RECURSIVE function intd_oss_(z0,z1,z2,z3,nf)
!вычисление интеграла G_1/y в осесимметричном случае
use gen_mod
complex(8) z,z1,z2,z3,z0
complex(8) z12,z23,z31
real(8) intd_oss_,intd_an
!real(8) ymin,ymax
integer(4) nf
!ymin=min(dimag(z1),dimag(z2),dimag(z3))
!ymax=max(dimag(z1),dimag(z2),dimag(z3))
!if (ymin<1.0d-6) ymin=1.0d-6
!if (ymax/ymin<2.0d0.or.ymax-ymin<1.0d-2) then
if (.true.) then
  z=(z1+z2+z3)/3.0d0
  intd_oss_=intd_an(z0,z1,z2,z3,nf)/dimag(z)
else
  z12=(z1+z2)*d5
  z23=(z2+z3)*d5
  z31=(z3+z1)*d5
  intd_oss_=intd_oss_(z0,z1,z12,z31,nf)+intd_oss_(z0,z12,z2,z23,nf)+intd_oss_(z0,z12,z23,z31,nf)+intd_oss_(z0,z23,z3,z31,nf)
endif
end

function test_point_in_tr(z,z1,z2,z3)
!обход точек треугольника такой, что его внутренность слева
use gen_mod
logical test_point_in_tr
complex(8) z,z1,z2,z3,a,b
real(8) v
test_point_in_tr=.false.
a=z2-z1
b=z-z1
v=dimag(dconjg(a)*b)
if (v<d0) return
a=z3-z2
b=z-z2
v=dimag(dconjg(a)*b)
if (v<d0) return
a=z1-z3
b=z-z3
v=dimag(dconjg(a)*b)
if (v<d0) return
test_point_in_tr=.true.
end

function intd_an(z0,z1,z2,z3,nf)
use gen_mod
integer(4) nf
real(8) intd_an,val,val2
complex(8) z0,z1,z2,z3
call intd_an_2(z0,z1,z2,z3,nf,val,val2,.false.)
intd_an=val
end

subroutine intd_an_2(z0,z1,z2,z3,nf,val,val2,use_dual)
!аналитический интеграл по треугольнику (z1,z2,z3) от функции G_1
!порядок вершин треугольника такой, чтобы при обходе внутренность оставалась слева
!nf= 0 - G_1
    !2 - G_2
    !4 - dG_1/dx0 
    !5 - dG_1/dy0
   !8 - dG_2/dx
   !9 - dG_2/dy
use gen_mod
integer(4) nf
real(8) intd_an_out,vp12,vp23,vp31,intd_an_corn,sp12,sp23,sp31,val,val2,val3,val4,val5,val6,vect_p
complex(8) z0,z1,z2,z3,a,b,c
logical use_dual
val=d0
val2=d0
!проверяеем, что треугольник не вырожденный
a=z1-z2
b=z1-z3
vp12=vect_p(a,b) 
if (dabs(vp12)<1.0d-8) return
!вычисляем для невырожденного треугольника
if (cdabs(z0-z1)<1.0d-8) then
  !z0 совпадает с z1
  if (use_dual) then
    call intd_an_corn_2(z1,z2,z3,nf,val,val2)
  else
    val=intd_an_corn(z1,z2,z3,nf)
  endif
elseif (cdabs(z0-z2)<1.0d-8) then
  !z0 совпадает с z2
  if (use_dual) then
    call intd_an_corn_2(z2,z3,z1,nf,val,val2)
  else
    val=intd_an_corn(z2,z3,z1,nf)
  endif
elseif (cdabs(z0-z3)<1.0d-8) then
  !z0 совпадает с z3
  if (use_dual) then
    call intd_an_corn_2(z3,z1,z2,nf,val,val2)
  else
    val=intd_an_corn(z3,z1,z2,nf)
  endif
else
  a=z0-z1
  b=z2-z1
  c=z0-z2
  vp12=dimag(dconjg(a)*b)
  sp12=dreal(dconjg(a)*c)
  if (cdabs(vp12/a)<1.0d-8.and.sp12<d0) then
    !z0 на отрезке (z1,z2)
    if (use_dual) then
      call intd_an_corn_2(z0,z3,z1,nf,val,val2)
      call intd_an_corn_2(z0,z2,z3,nf,val3,val4)
      val=val+val3
      val2=val2+val4
    else
	    val=intd_an_corn(z0,z3,z1,nf)+intd_an_corn(z0,z2,z3,nf)
    endif
  else
    a=z0-z2
    b=z3-z2
    c=z0-z3
    vp23=dimag(dconjg(a)*b)
    sp23=dreal(dconjg(a)*c)
	  if (cdabs(vp23/a)<1.0d-8.and.sp23<d0) then
	    !z0 на отрезке (z2,z3)
      if (use_dual) then
        call intd_an_corn_2(z0,z1,z2,nf,val,val2)
        call intd_an_corn_2(z0,z3,z1,nf,val3,val4)
        val=val+val3
        val2=val2+val4
      else
	      val=intd_an_corn(z0,z1,z2,nf)+intd_an_corn(z0,z3,z1,nf)
      endif
    else
      a=z0-z3
      b=z1-z3
      c=z0-z1
      vp31=dimag(dconjg(a)*b)
      sp31=dreal(dconjg(a)*c)
	    if (cdabs(vp31/a)<1.0d-8.and.sp31<d0) then
	      !z0 на отрезке (z3,z1)
        if (use_dual) then
          call intd_an_corn_2(z0,z2,z3,nf,val,val2)
          call intd_an_corn_2(z0,z1,z2,nf,val3,val4)
          val=val+val3
          val2=val2+val4
        else
		      val=intd_an_corn(z0,z2,z3,nf)+intd_an_corn(z0,z1,z2,nf)
        endif
      else
	      if (vp12<0.and.vp23<0.and.vp31<0) then
		      !z0 внутри треугольника
          if (use_dual) then
            call intd_an_corn_2(z0,z1,z2,nf,val,val2)
            call intd_an_corn_2(z0,z2,z3,nf,val3,val4)
            call intd_an_corn_2(z0,z3,z1,nf,val5,val6)
            val=val+val3+val5
            val2=val2+val4+val6
          else  
		        val=intd_an_corn(z0,z1,z2,nf)+intd_an_corn(z0,z2,z3,nf)+intd_an_corn(z0,z3,z1,nf)
          endif
		    else
		      !z0 снаружи треугольника
          if (use_dual) then
            call intd_an_out_2(z0,z1,z2,z3,nf,val,val2)
          else
            val=intd_an_out(z0,z1,z2,z3,nf)
          endif
		    endif
      endif
	  endif
  endif
endif
end

function intd_an_out(z0,z1,z2,z3,nf)
!аналитический интеграл по треугольнику (z1,z2,z3)
!z0 снаружи треугольника
!nf= 0 - G_1
    !2 - G_2
    !4 - dG_1/dx0 
    !5 - dG_1/dy0
   !8 - dG_2/dx
   !9 - dG_2/dy
use gen_mod
integer(4) nf
real(8) intd_an_out,intd_an0_out,intd_an45_out,intd_an2_out,intd_an89_out
complex(8) z0,z1,z2,z3
select case (nf)
case (0)
  intd_an_out=intd_an0_out(z0,z1,z2,z3)
case (2)
  intd_an_out=intd_an2_out(z0,z1,z2,z3)
case (4,5)
  intd_an_out=intd_an45_out(z0,z1,z2,z3,nf)
case (8,9)
  intd_an_out=intd_an89_out(z0,z1,z2,z3,nf)
case default
  call gs_print_stop("Error intd_an_out!")
end select
end 

subroutine intd_an_out_2(z0,z1,z2,z3,nf,val,val2)
!аналитический интеграл по треугольнику (z1,z2,z3)
!z0 снаружи треугольника
!nf= 0 - G_1
    !2 - G_2
    !4 - dG_1/dx0 
    !5 - dG_1/dy0
   !8 - dG_2/dx
   !9 - dG_2/dy
use gen_mod
integer(4) nf
real(8) val,val2
complex(8) z0,z1,z2,z3
select case (nf)
case (4)
  call intd_an45_out_2(z0,z1,z2,z3,val,val2)
case (8)
  call intd_an89_out_2(z0,z1,z2,z3,val,val2)
case default
  call gs_print_stop("Error intd_an_out_2!")
end select
end 

function intd_an_corn(z1,z2,z3,nf)
!аналитический интеграл по треугольнику (z1,z2,z3)
!z0 совпадает с z1
!nf= 0 - G_1
    !2 - G_2
    !4 - dG_1/dx0 
    !5 - dG_1/dy0
   !8 - dG_2/dx
   !9 - dG_2/dy
use gen_mod
integer(4) nf
real(8) intd_an_corn,intd_an0_corn,intd_an45_corn,intd_an2_corn,intd_an89_corn
complex(8) z1,z2,z3
select case (nf)
case (0)
  intd_an_corn=intd_an0_corn(z1,z2,z3)
case (2)
  intd_an_corn=intd_an2_corn(z1,z2,z3)
case (4,5)
  intd_an_corn=intd_an45_corn(z1,z2,z3,nf)
case (8,9)
  intd_an_corn=intd_an89_corn(z1,z2,z3,nf)
case default
  call gs_print_stop("Error intd_an_corn!")
end select
end

subroutine intd_an_corn_2(z1,z2,z3,nf,val,val2)
!аналитический интеграл по треугольнику (z1,z2,z3)
!z0 совпадает с z1
!nf= 0 - G_1
    !2 - G_2
    !4 - dG_1/dx0 
    !5 - dG_1/dy0
   !8 - dG_2/dx
   !9 - dG_2/dy
use gen_mod
integer(4) nf
real(8) val,val2
complex(8) z1,z2,z3
select case (nf)
case (4)
  call intd_an45_corn_2(z1,z2,z3,val,val2)
case (8)
  call intd_an89_corn_2(z1,z2,z3,val,val2)
case default
  call gs_print_stop("Error intd_an_corn_2!")
end select
end

function intd_an0_out(z0,z1,z2,z3)
!аналитический интеграл по треугольнику (z1,z2,z3) от функции G_1
!z0 снаружи треугольника
use gen_mod
real(8) intd_an0_out,tet,gam,delta_n12,delta_k12,delta_n13,delta_k13
complex(8) z0,z1,z2,z3,z21,z13,z32,lnn12,lnk12,lnn13,lnk13
complex(8) fzk12,fzn12,fzk13,fzn13,fz,z30,z20,z10
z21=z2-z1
z13=z1-z3
z32=z3-z2
z30=z3-z0
z10=z1-z0
z20=z2-z0
tet=zarg(z21)
gam=zarg(z32)
call get_delta_l(z0,z1,z2,delta_n12,delta_k12)
call get_delta_l(z0,z1,z3,delta_n13,delta_k13)
lnn12=dlog(cdabs(z10))+ii*delta_n12
lnn13=dlog(cdabs(z10))+ii*delta_n13
lnk12=dlog(cdabs(z20))+ii*delta_k12
lnk13=dlog(cdabs(z30))+ii*delta_k13
fzn12=z10**2*(d2*lnn12-3.0d0)
fzk12=z20**2*(d2*lnk12-3.0d0)
fzn13=z10**2*(d2*lnn13-3.0d0)
fzk13=z30**2*(d2*lnk13-3.0d0)
fz=z21/z13*(fzk13-fzn13)+(fzk12-fzn12)
intd_an0_out=-dreal(dabs(dsin(tet-gam))*cdexp(-ii*(tet+gam))*fz*0.25d0)
end

function intd_an0_corn(z1,z2,z3)
!аналитический интеграл по треугольнику (z1,z2,z3) от функции G_1
!z0 совпадает с z1
use gen_mod
real(8) intd_an0_corn,tet,gam,delta_n,delta_k
complex(8) z1,z2,z3,z21,z31,z32,lnn,lnk
complex(8) fzk,fzn,fz
z21=z2-z1
z31=z3-z1
z32=z3-z2
tet=zarg(z21)
gam=zarg(z32)
call get_delta_l(z1,z2,z3,delta_n,delta_k)
lnn=dlog(cdabs(z21))+ii*delta_n
lnk=dlog(cdabs(z31))+ii*delta_k
fzn=z21*(d2*lnn-3.0d0)
fzk=z31*(d2*lnk-3.0d0)
fz=z21*(fzk-fzn)
intd_an0_corn=dreal(dabs(dsin(tet-gam))*cdexp(-ii*(tet+gam))*fz*0.25d0)
end

function intd_an2_out(z0,z1,z2,z3)
!аналитический интеграл по треугольнику (z1,z2,z3) от функции G_2
!z0 снаружи треугольника
use gen_mod
real(8) intd_an2_out,gam1,gam2,gam3,delta_n12,delta_k12,delta_n13,delta_k13
complex(8) z0,z1,z2,z3,z21,z13,z32,lnn12,lnk12,lnn13,lnk13
complex(8) fzk,fzn,ii1,ii2,ii3,ii4,ii5,z30,z20,z10
z21=z2-z1
z13=z1-z3
z32=z3-z2
z30=z3-z0
z10=z1-z0
z20=z2-z0
gam1=zarg(z21)
gam2=zarg(z32)
gam3=zarg(z13)
call get_delta_l(z0,z1,z2,delta_n12,delta_k12)
call get_delta_l(z0,z1,z3,delta_n13,delta_k13)
lnn12=dlog(cdabs(z10))+ii*delta_n12
lnn13=dlog(cdabs(z10))+ii*delta_n13
lnk12=dlog(cdabs(z20))+ii*delta_k12
lnk13=dlog(cdabs(z30))+ii*delta_k13
fzn=dconjg(z10)*z10**3*(6.0d0*lnn13-11.0d0)
fzk=dconjg(z30)*z30**3*(6.0d0*lnk13-11.0d0)
ii1=fzk-fzn
fzn=z10**4*(12.0d0*lnn13-25.0d0)
fzk=z30**4*(12.0d0*lnk13-25.0d0)
!ii2=(cdexp(-d2*ii*gam2)+cdexp(-d2*ii*gam3))*0.125d0*(fzk-fzn)
ii2=cdexp(-ii*(gam2+gam3))*d2*dcos(gam2-gam3)*0.125d0*(fzk-fzn)
fzn=dconjg(z10)*z10**3*(6.0d0*lnn12-11.0d0)
fzk=dconjg(z20)*z20**3*(6.0d0*lnk12-11.0d0)
ii3=fzk-fzn
fzn=z10**4*(12.0d0*lnn12-25.0d0)
fzk=z20**4*(12.0d0*lnk12-25.0d0)
!ii4=(cdexp(-d2*ii*gam1)+cdexp(-d2*ii*gam2))*0.125d0*(fzk-fzn)
ii4=cdexp(-ii*(gam1+gam2))*d2*dcos(gam1-gam2)*0.125d0*(fzk-fzn)
ii5=-cdexp(-ii*(gam1+gam2))*(z21/z13*(ii1-ii2)+ii3-ii4)
intd_an2_out=dabs(dsin(gam2-gam1))*dreal(ii5)/144.0d0
end

function intd_an2_corn(z1,z2,z3)
!аналитический интеграл по треугольнику (z1,z2,z3) от функции G_2
!z0 совпадает с z1
use gen_mod
real(8) intd_an2_corn,gam1,gam2,gam3,delta_n,delta_k
complex(8) z1,z2,z3,z21,z13,z32,lnn,lnk,z31
complex(8) fzk,fzn,ii1,ii2,ii3
z21=z2-z1
z13=z1-z3
z31=-z13
z32=z3-z2
gam1=zarg(z21)
gam2=zarg(z32)
gam3=zarg(z13)
call get_delta_l(z1,z2,z3,delta_n,delta_k)
lnn=dlog(cdabs(z21))+ii*delta_n
lnk=dlog(cdabs(z31))+ii*delta_k
fzn=dconjg(z21)*z21**2*(6.0d0*lnn-11.0d0)
fzk=dconjg(z31)*z31**2*(6.0d0*lnk-11.0d0)
ii1=fzk-fzn
fzn=z21**3*(12.0d0*lnn-25.0d0)*cdexp(-ii*(gam1+gam2))*d2*dcos(gam1-gam2)
fzk=z31**3*(12.0d0*lnk-25.0d0)*cdexp(-ii*(gam2+gam3))*d2*dcos(gam2-gam3)
ii2=0.125d0*(fzk-fzn)
ii3=cdexp(-ii*(gam1+gam2))*z21*(ii1-ii2)
intd_an2_corn=dabs(dsin(gam2-gam1))*dreal(ii3)/144.0d0
end

function intd_an45_outz(z0,z1,z2,z3)
!аналитический интеграл по треугольнику (z1,z2,z3) от функции dG_1/dx0 или dG_1/dy0
!z0 снаружи треугольника
!nf= 4 - dG_1/dx0, 5 - dG_1/dy0
use gen_mod
real(8) tet,gam,delta_n12,delta_k12,delta_n13,delta_k13
complex(8) z0,z1,z2,z3,z21,z13,z32,lnn12,lnk12,lnn13,lnk13
complex(8) fzk12,fzn12,fzk13,fzn13,fz,z30,z20,z10,intd_an45_outz
z21=z2-z1
z13=z1-z3
z32=z3-z2
z30=z3-z0
z10=z1-z0
z20=z2-z0
tet=zarg(z21)
gam=zarg(z32)
call get_delta_l(z0,z1,z2,delta_n12,delta_k12)
call get_delta_l(z0,z1,z3,delta_n13,delta_k13)
lnn12=dlog(cdabs(z10))+ii*delta_n12
lnn13=dlog(cdabs(z10))+ii*delta_n13
lnk12=dlog(cdabs(z20))+ii*delta_k12
lnk13=dlog(cdabs(z30))+ii*delta_k13
fzn12=z10*(lnn12-d1)
fzk12=z20*(lnk12-d1)
fzn13=z10*(lnn13-d1)
fzk13=z30*(lnk13-d1)
fz=z21/z13*(fzk13-fzn13)+(fzk12-fzn12)
intd_an45_outz=dabs(dsin(tet-gam))*cdexp(-ii*(tet+gam))*fz  !1/4 сократилась с 4 из производной
end

function intd_an45_out(z0,z1,z2,z3,nf)
!аналитический интеграл по треугольнику (z1,z2,z3) от функции dG_1/dx0 или dG_1/dy0
!z0 снаружи треугольника
!nf= 4 - dG_1/dx0, 5 - dG_1/dy0
use gen_mod
integer(4) nf
real(8) intd_an45_out
complex(8) z0,z1,z2,z3
complex(8) fz,intd_an45_outz
fz=intd_an45_outz(z0,z1,z2,z3)
if (nf==4) then
  intd_an45_out=dreal(fz)  
else
  intd_an45_out=-dimag(fz)
endif
end

subroutine intd_an45_out_2(z0,z1,z2,z3,val,val2)
!аналитический интеграл по треугольнику (z1,z2,z3) от функции dG_1/dx0 или dG_1/dy0
!z0 снаружи треугольника
!nf= 4 - dG_1/dx0, 5 - dG_1/dy0
use gen_mod
real(8) val,val2
complex(8) z0,z1,z2,z3
complex(8) fz,intd_an45_outz
fz=intd_an45_outz(z0,z1,z2,z3)
val=dreal(fz)  
val2=-dimag(fz)
end

function intd_an45_cornz(z1,z2,z3)
!аналитический интеграл по треугольнику (z1,z2,z3) от функции от функции dG_1/dx0 или dG_1/dy0
!z0 совпадает с z1
!nf= 4 - dG_1/dx0, 5 - dG_1/dy0
!по упрощенным формулам, полученным при написании статьи
use gen_mod
real(8) tet,gam,delta_n,delta_k
complex(8) z1,z2,z3,z21,z13,z32,lnn,lnk
complex(8) fz,intd_an45_cornz
z21=z2-z1
z32=z3-z2
z13=z1-z3
tet=zarg(z21)
gam=zarg(z32)
call get_delta_l(z1,z2,z3,delta_n,delta_k)
lnn=dlog(cdabs(z21))+ii*delta_n
lnk=dlog(cdabs(z13))+ii*delta_k
fz=z21*(lnk-lnn)
intd_an45_cornz=-dabs(dsin(tet-gam))*cdexp(-ii*(tet+gam))*fz
end

function intd_an45_corn(z1,z2,z3,nf)
!аналитический интеграл по треугольнику (z1,z2,z3) от функции от функции dG_1/dx0 или dG_1/dy0
!z0 совпадает с z1
!nf= 4 - dG_1/dx0, 5 - dG_1/dy0
!по упрощенным формулам, полученным при написании статьи
use gen_mod
integer(4) nf
real(8) intd_an45_corn
complex(8) z1,z2,z3
complex(8) fz,intd_an45_cornz
fz=intd_an45_cornz(z1,z2,z3)
if (nf==4) then
  intd_an45_corn=dreal(fz)  
else
  intd_an45_corn=-dimag(fz)
endif
end

subroutine intd_an45_corn_2(z1,z2,z3,val,val2)
!аналитический интеграл по треугольнику (z1,z2,z3) от функции от функции dG_1/dx0 или dG_1/dy0
!z0 совпадает с z1
!nf= 4 - dG_1/dx0, 5 - dG_1/dy0
!по упрощенным формулам, полученным при написании статьи
use gen_mod
real(8) val,val2
complex(8) z1,z2,z3
complex(8) fz,intd_an45_cornz
fz=intd_an45_cornz(z1,z2,z3)
val=dreal(fz)  
val2=-dimag(fz)
end

subroutine intd_an89_outz(z0,z1,z2,z3,ii1,ii2,ii3,ii4,ii5,ii6,fz1,fz2,fr)
!аналитический интеграл по треугольнику (z1,z2,z3) от функции dG_2/dx0 или dG_2/dy0
!z0 снаружи треугольника
!nf= 8 - dG_2/dx0, 9 - dG_2/dy0
use gen_mod
real(8) gam1,gam2,gam3,delta_n12,delta_k12,delta_n13,delta_k13,fr
complex(8) z0,z1,z2,z3,z21,z13,z32,lnn12,lnk12,lnn13,lnk13
complex(8) fzk,fzn,ii1,ii2,ii3,ii4,ii5,ii6,z30,z20,z10,fz1,fz2
z21=z2-z1
z13=z1-z3
z32=z3-z2
z30=z3-z0
z10=z1-z0
z20=z2-z0
gam1=zarg(z21)
gam2=zarg(z32)
gam3=zarg(z13)
call get_delta_l(z0,z1,z2,delta_n12,delta_k12)
call get_delta_l(z0,z1,z3,delta_n13,delta_k13)
lnn12=dlog(cdabs(z10))+ii*delta_n12
lnn13=dlog(cdabs(z10))+ii*delta_n13
lnk12=dlog(cdabs(z20))+ii*delta_k12
lnk13=dlog(cdabs(z30))+ii*delta_k13
fzn=z10**3*(6.0d0*lnn13-11.0d0)
fzk=z30**3*(6.0d0*lnk13-11.0d0)
ii1=fzk-fzn
ii3=ii1*cdexp(-ii*(gam2+gam3))*d2*dcos(gam2-gam3)
fzn=dconjg(z10)*z10**2*(d2*lnn13-3.0d0)
fzk=dconjg(z30)*z30**2*(d2*lnk13-3.0d0)
ii2=(fzk-fzn)*9.0d0
fzn=z10**3*(6.0d0*lnn12-11.0d0)
fzk=z20**3*(6.0d0*lnk12-11.0d0)
ii4=fzk-fzn
ii6=ii4*cdexp(-ii*(gam1+gam2))*d2*dcos(gam1-gam2)
fzn=dconjg(z10)*z10**2*(d2*lnn12-3.0d0)
fzk=dconjg(z20)*z20**2*(d2*lnk12-3.0d0)
ii5=(fzk-fzn)*9.0d0
fz1=cdexp(-ii*(gam1+gam2))
fz2=z21/z13
fr=dabs(dsin(gam2-gam1))/144.0d0
!if (nf==8) then
!  ii7=cdexp(-ii*(gam1+gam2))*(z21/z13*(ii1+ii2-ii3)+ii4+ii5-ii6)
!  intd_an89_out=dabs(dsin(gam2-gam1))*dreal(ii7)/144.0d0  
!else
!  ii7=cdexp(-ii*(gam1+gam2))*(z21/z13*(ii1-ii2+ii3)+ii4-ii5+ii6)
!  intd_an89_out=dabs(dsin(gam2-gam1))*dimag(ii7)/144.0d0  
!endif
end

function intd_an89_out(z0,z1,z2,z3,nf)
!аналитический интеграл по треугольнику (z1,z2,z3) от функции dG_2/dx0 или dG_2/dy0
!z0 снаружи треугольника
!nf= 8 - dG_2/dx0, 9 - dG_2/dy0
use gen_mod
integer(4) nf
real(8) intd_an89_out,fr
complex(8) z0,z1,z2,z3
complex(8) ii1,ii2,ii3,ii4,ii5,ii6,ii7,fz1,fz2
call intd_an89_outz(z0,z1,z2,z3,ii1,ii2,ii3,ii4,ii5,ii6,fz1,fz2,fr)
if (nf==8) then
  ii7=fz1*(fz2*(ii1+ii2-ii3)+ii4+ii5-ii6)
  intd_an89_out=fr*dreal(ii7)
else
  ii7=fz1*(fz2*(ii1-ii2+ii3)+ii4-ii5+ii6)
  intd_an89_out=fr*dimag(ii7)
endif
end

subroutine intd_an89_out_2(z0,z1,z2,z3,val,val2)
!аналитический интеграл по треугольнику (z1,z2,z3) от функции dG_2/dx0 или dG_2/dy0
!z0 снаружи треугольника
!nf= 8 - dG_2/dx0, 9 - dG_2/dy0
use gen_mod
real(8) val,val2,fr
complex(8) z0,z1,z2,z3
complex(8) ii1,ii2,ii3,ii4,ii5,ii6,ii7,fz1,fz2
call intd_an89_outz(z0,z1,z2,z3,ii1,ii2,ii3,ii4,ii5,ii6,fz1,fz2,fr)
ii7=fz1*(fz2*(ii1+ii2-ii3)+ii4+ii5-ii6)
val=fr*dreal(ii7)
ii7=fz1*(fz2*(ii1-ii2+ii3)+ii4-ii5+ii6)
val2=fr*dimag(ii7)
end

subroutine intd_an89_cornz(z1,z2,z3,ii1,ii2,ii3,fz,z21,fr)
!аналитический интеграл по треугольнику (z1,z2,z3) от функции dG_2/dx0 или dG_2/dy0
!z0 совпадает с z1
!nf= 8 - dG_2/dx0, 9 - dG_2/dy0
use gen_mod
real(8) gam1,gam2,gam3,delta_n,delta_k,fr
complex(8) z1,z2,z3,z21,z13,z32,z31,lnn,lnk
complex(8) fzk,fzn,ii1,ii2,ii3,fz
z21=z2-z1
z31=z3-z1
z13=-z31
z32=z3-z2
gam1=zarg(z21)
gam2=zarg(z32)
gam3=zarg(z13)
call get_delta_l(z1,z2,z3,delta_n,delta_k)
lnn=dlog(cdabs(z21))+ii*delta_n
lnk=dlog(cdabs(z31))+ii*delta_k
fzn=z21**2*(6.0d0*lnn-11.0d0)
fzk=z31**2*(6.0d0*lnk-11.0d0)
ii1=fzk-fzn
fzn=dconjg(z21)*z21*(d2*lnn-3.0d0)
fzk=dconjg(z31)*z31*(d2*lnk-3.0d0)
ii2=(fzk-fzn)*9.0d0
fzn=z21**2*(6.0d0*lnn-11.0d0)*cdexp(-ii*(gam1+gam2))*d2*dcos(gam1-gam2)
fzk=z31**2*(6.0d0*lnk-11.0d0)*cdexp(-ii*(gam2+gam3))*d2*dcos(gam2-gam3)
ii3=(fzk-fzn)
fz=cdexp(-ii*(gam1+gam2))
fr=dabs(dsin(gam2-gam1))/144.0d0
!if (nf==8) then
!  ii4=cdexp(-ii*(gam1+gam2))*z21*(ii1+ii2-ii3)
!  intd_an89_corn=-dabs(dsin(gam2-gam1))*dreal(ii4)/144.0d0  
!else
!  ii4=cdexp(-ii*(gam1+gam2))*z21*(-ii1+ii2-ii3)
!  intd_an89_corn=dabs(dsin(gam2-gam1))*dimag(ii4)/144.0d0  
!endif
end

function intd_an89_corn(z1,z2,z3,nf)
!аналитический интеграл по треугольнику (z1,z2,z3) от функции dG_2/dx0 или dG_2/dy0
!z0 совпадает с z1
!nf= 8 - dG_2/dx0, 9 - dG_2/dy0
use gen_mod
integer(4) nf
real(8) intd_an89_corn,fr
complex(8) z1,z2,z3,z21
complex(8) ii1,ii2,ii3,ii4,fz
call intd_an89_cornz(z1,z2,z3,ii1,ii2,ii3,fz,z21,fr)
if (nf==8) then
  ii4=fz*z21*(ii1+ii2-ii3)
  intd_an89_corn=-fr*dreal(ii4)
else
  ii4=fz*z21*(-ii1+ii2-ii3)
  intd_an89_corn=fr*dimag(ii4)
endif
end

subroutine intd_an89_corn_2(z1,z2,z3,val,val2)
!аналитический интеграл по треугольнику (z1,z2,z3) от функции dG_2/dx0 или dG_2/dy0
!z0 совпадает с z1
!nf= 8 - dG_2/dx0, 9 - dG_2/dy0
use gen_mod
real(8) val,val2,fr
complex(8) z1,z2,z3,z21
complex(8) ii1,ii2,ii3,ii4,fz
call intd_an89_cornz(z1,z2,z3,ii1,ii2,ii3,fz,z21,fr)
ii4=fz*z21*(ii1+ii2-ii3)
val=-fr*dreal(ii4)
ii4=fz*z21*(-ii1+ii2-ii3)
val2=fr*dimag(ii4)
end

function calcint__(z0,z1,z2,nf)
!nf 0 - G_1
!   1 - G'_1=dG_1/dn
!   2 - G_2
!   3 - G'_2=dG_2/dn
!   4 - dG_1/dx_0
!   5 - dG_1/dy_0
!   6 - dG'_1/dx_0
!   7 - dG'_1/dy_0
!   8 - dG_2/dx_0
!   9 - dG_2/dy_0
!   10- dG'_2/dx_0
!   11- dG'_2/dy_0
use pgmod
integer(4) i,nf,nint
real(8) xx,yy,rr,xt,yt,s,ff,costt,sintt,calcint__
real(8) ds0,it
complex(8) fz1,zt,z0,z1,z2,dz,z21
z21=z2-z1
costt=dreal(z21)/cdabs(z21)
sintt=dimag(z21)/cdabs(z21)
ds0=1.0d-6
nint=cdabs(z21)/ds0
dz=z21/nint
s=d0
do i=1,nint
  it=i-d5
  zt=z1+dz*it
  xt=dreal(zt)
  yt=dimag(zt)
  xx=xt-dreal(z0)
  yy=yt-dimag(z0)
  rr=xx*xx+yy*yy
  if (rr<1.0d-16) cycle
  if (nf.eq.0) then
    ff=d5*dlog(rr)
  elseif (nf.eq.1) then
	ff=(xx*sintt-yy*costt)/rr
  elseif (nf.eq.2) then
    ff=rr*(d5*dlog(rr)-d1)*0.25d0
  elseif (nf.eq.3) then
    ff=(xx*sintt-yy*costt)*(dlog(rr)-d1)*0.25d0
  elseif (nf.eq.4) then
    ff=-xx/rr
  elseif (nf.eq.5) then
    ff=-yy/rr
  elseif (nf.eq.6) then
    zt=dcmplx(xx,yy)
    fz1=(costt+ii*sintt)/(zt**2)
    ff=dimag(fz1)
  elseif (nf.eq.7) then
    zt=dcmplx(xx,yy)
    fz1=(costt+ii*sintt)/(zt**2)
    ff=dreal(fz1)
  elseif (nf.eq.8) then
    ff=-xx*(dlog(rr)-d1)*0.25d0
  elseif (nf.eq.9) then
    ff=-yy*(dlog(rr)-d1)*0.25d0
  elseif (nf.eq.10) then
    ff=(-sintt*(dlog(rr)-d1+d2*xx*xx/rr)+costt*d2*xx*yy/rr)*0.25d0
  elseif (nf.eq.11) then
    ff=(costt*(dlog(rr)-d1+d2*yy*yy/rr)-sintt*d2*xx*yy/rr)*0.25d0
  else
    call gs_print_stop("Error calcint__!")
  endif
  s=s+ff
enddo
calcint__=s*cdabs(z21)/nint
end

subroutine test_intd
use gen_mod
complex(8) z,z1,z2,z3
real(8) f1,f2,intd_an,f3,ds,f4,f5,f6
real(8) intd_num
integer(4) nf,nf2
nf=2 !0,2
nf2=4
if (nf==2) nf2=8
ds=1.0d-6
z1=dcmplx(d0,d0)
z2=dcmplx(d1,d0)
z3=dcmplx(d0,d1)
z=dcmplx(0.2d0,0.5d0) 
!f1=intd_num(z,z1,z2,z3,nf,0)
!f2=intd_an(z,z1,z2,z3,nf)
!f3=f1-f2

f1=(intd_an(z+c1*ds,z1,z2,z3,nf)-intd_an(z,z1,z2,z3,nf))/ds
f2=intd_an(z,z1,z2,z3,nf2)
f3=f2-f1
f4=(intd_an(z+ii*ds,z1,z2,z3,nf)-intd_an(z,z1,z2,z3,nf))/ds
f5=intd_an(z,z1,z2,z3,nf2+1)
f6=f5-f4

f1=(intd_num(z+c1*ds,z1,z2,z3,nf,0)-intd_num(z,z1,z2,z3,nf,0))/ds
f2=intd_num(z,z1,z2,z3,nf2,0)
f3=f2-f1
f4=(intd_num(z+ii*ds,z1,z2,z3,nf,0)-intd_num(z,z1,z2,z3,nf,0))/ds
f5=intd_num(z,z1,z2,z3,nf2+1,0)
f6=f5-f4
end

subroutine test_intl
use gen_mod
complex(8) z,z1,z2,zt
real(8) f1,calcint__,f2,calcint_nf1011,f3
integer(4) i
z1=dcmplx(d0,d0)
z2=dcmplx(d0,d1)
z=dcmplx(0.0d0,0.4d0)
do i=1,3
  zt=z+c1*0.00001d0*(i-d2)
  f1=calcint__(zt,z1,z2,11)
  f2=calcint_nf1011(zt,z1,z2,11) 
  f3=f1-f2
  write(*,"(F9.5, ' ', F9.5, ' ', F9.5, ' ', E15.8)") dreal(zt),f1,f2,f3
enddo
end
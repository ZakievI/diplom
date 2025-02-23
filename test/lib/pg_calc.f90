function calcintn(k,x0,y0,nf,ia,knd)
!численное вычисление интегралов
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
!   12- G_3(Y0)  функция Бесселя
!   13- G_3'(Y0) производная по нормали от функции Бесселя
!   14- G_4(K0)  функция Бесселя
!   15- G_4'(K0) производная по нормали от функции Бесселя
!   16- dG_4/dx_0   
!   17- dG_4/dy_0        
!   18- dG_4'/dx_0
!   19- dG_4'/dy_0
use pgmod
integer(4) k,i,nf,knd,nint,ia
real(8) x0,y0,dx,dy,xx,yy,rr,rr2,xt,yt,s,ff,costt,sintt,calcintn
real(8) ds0,it,ff2,eta0 !,eta
complex(8) fz1,zt,z0
logical collenear3z,is_analit
type(TBound), pointer :: b
b=>gs%a(ia)%bnd(knd)
costt=dreal(b%ett(k))
sintt=dimag(b%ett(k))
!ds0=1.0d-3
ds0=gs%const%ds0numint
nint=b%l(k)/ds0
if (mod(nint,2).eq.1) nint=nint+1
dx=(b%x(k+1)-b%x(k))/nint
dy=(b%y(k+1)-b%y(k))/nint
s=d0
z0=dcmplx(x0,y0)
is_analit=.false.
eta0=dimag((z0-b%z(k))*dconjg(b%ett(k)))
select case (nf)
case (0)
  do i=1,nint
    it=i-d5
    xt=b%x(k)+dx*it
    yt=b%y(k)+dy*it
    xx=xt-x0
    yy=yt-y0
    rr=xx*xx+yy*yy
    if (dabs(rr)<1.0d-16) cycle
    ff=d5*dlog(rr)
    s=s+ff
  enddo
case (1)
  do i=1,nint
    it=i-d5
    xt=b%x(k)+dx*it
    yt=b%y(k)+dy*it
    xx=xt-x0
    yy=yt-y0
    rr=xx*xx+yy*yy
    if (dabs(rr)<1.0d-16) cycle
    ff=(xx*sintt-yy*costt)/rr
    s=s+ff
  enddo
case (2)
  do i=1,nint
    it=i-d5
    xt=b%x(k)+dx*it
    yt=b%y(k)+dy*it
    xx=xt-x0
    yy=yt-y0
    rr=xx*xx+yy*yy
    if (dabs(rr)<1.0d-16) cycle
    ff=rr*(d5*dlog(rr)-d1)*0.25d0
    s=s+ff
  enddo
case (3)
  do i=1,nint
    it=i-d5
    xt=b%x(k)+dx*it
    yt=b%y(k)+dy*it
    xx=xt-x0
    yy=yt-y0
    rr=xx*xx+yy*yy
    if (dabs(rr)<1.0d-16) cycle
    ff=(xx*sintt-yy*costt)*(dlog(rr)-d1)*0.25d0
    s=s+ff
  enddo
case (4)
  do i=1,nint
    it=i-d5
    xt=b%x(k)+dx*it
    yt=b%y(k)+dy*it
    xx=xt-x0
    yy=yt-y0
    rr=xx*xx+yy*yy
    if (dabs(rr)<1.0d-16) cycle
    ff=-xx/rr
    s=s+ff
  enddo
case (5)
  do i=1,nint
    it=i-d5
    xt=b%x(k)+dx*it
    yt=b%y(k)+dy*it
    xx=xt-x0
    yy=yt-y0
    rr=xx*xx+yy*yy
    if (dabs(rr)<1.0d-16) cycle
    ff=-yy/rr
    s=s+ff
  enddo
case (6)
  do i=1,nint
    it=i-d5
    xt=b%x(k)+dx*it
    yt=b%y(k)+dy*it
    xx=xt-x0
    yy=yt-y0
    rr=xx*xx+yy*yy
    if (dabs(rr)<1.0d-16) cycle
    zt=dcmplx(xx,yy)
    fz1=b%ett(k)/(zt**2)
    ff=dimag(fz1)
    s=s+ff
  enddo
case (7)
  do i=1,nint
    it=i-d5
    xt=b%x(k)+dx*it
    yt=b%y(k)+dy*it
    xx=xt-x0
    yy=yt-y0
    rr=xx*xx+yy*yy
    if (dabs(rr)<1.0d-16) cycle
    zt=dcmplx(xx,yy)
    fz1=b%ett(k)/(zt**2)
    ff=dreal(fz1)
    s=s+ff
  enddo
case (8)
  do i=1,nint
    it=i-d5
    xt=b%x(k)+dx*it
    yt=b%y(k)+dy*it
    xx=xt-x0
    yy=yt-y0
    rr=xx*xx+yy*yy
    if (dabs(rr)<1.0d-16) cycle
    ff=-xx*(dlog(rr)-d1)*0.25d0
    s=s+ff
  enddo
case (9)
  do i=1,nint
    it=i-d5
    xt=b%x(k)+dx*it
    yt=b%y(k)+dy*it
    xx=xt-x0
    yy=yt-y0
    rr=xx*xx+yy*yy
    if (dabs(rr)<1.0d-16) cycle
    ff=-yy*(dlog(rr)-d1)*0.25d0
    s=s+ff
  enddo
case (10)
  do i=1,nint
    it=i-d5
    xt=b%x(k)+dx*it
    yt=b%y(k)+dy*it
    xx=xt-x0
    yy=yt-y0
    rr=xx*xx+yy*yy
    if (dabs(rr)<1.0d-16) cycle
    ff=(-sintt*(dlog(rr)-d1+d2*xx*xx/rr)+costt*d2*xx*yy/rr)*0.25d0
    s=s+ff
  enddo
case (11)
  do i=1,nint
    it=i-d5
    xt=b%x(k)+dx*it
    yt=b%y(k)+dy*it
    xx=xt-x0
    yy=yt-y0
    rr=xx*xx+yy*yy
    if (dabs(rr)<1.0d-16) cycle
    ff=(costt*(dlog(rr)-d1+d2*yy*yy/rr)-sintt*d2*xx*yy/rr)*0.25d0
    s=s+ff
  enddo
case (12)
  do i=1,nint
    it=i-d5
    xt=b%x(k)+dx*it
    yt=b%y(k)+dy*it
    xx=xt-x0
    yy=yt-y0
    rr=xx*xx+yy*yy
    if (dabs(rr)<1.0d-16) cycle
    ff=dbsy0(gsareaconst%k_helm*dsqrt(rr))
    s=s+ff
  enddo
  s=s*pi5
case (13)
  if (gs_test_point_collenear_for_BesselG.and.collenear3z(b%z(k),b%z(k+1),z0)) then
    is_analit=.true.
  else
    do i=1,nint
      it=i-d5
      xt=b%x(k)+dx*it
      yt=b%y(k)+dy*it
      xx=xt-x0
      yy=yt-y0
      rr=xx*xx+yy*yy
      if (dabs(rr)<1.0d-16) cycle
      rr2=dsqrt(rr)
      !eta=xx*sintt-yy*costt
      !ff=dbsy1(gsareaconst%k_helm*rr2)*(xx*sintt-yy*costt)/rr2
      ff=dbsy1(gsareaconst%k_helm*rr2)*eta0/rr2
      s=s+ff
    enddo
    s=-s*pi5*gsareaconst%k_helm
  endif
case (14)
  do i=1,nint
    it=i-d5
    xt=b%x(k)+dx*it
    yt=b%y(k)+dy*it
    xx=xt-x0
    yy=yt-y0
    rr=xx*xx+yy*yy
    if (dabs(rr)<1.0d-16) cycle
    ff=dbsk0(gsareaconst%k_helm*dsqrt(rr))
    s=s+ff
  enddo
  s=-s
case (15)
  if (gs_test_point_collenear_for_BesselG.and.collenear3z(b%z(k),b%z(k+1),z0)) then
    is_analit=.true.
  else
    do i=1,nint
      it=i-d5
      xt=b%x(k)+dx*it
      yt=b%y(k)+dy*it
      xx=xt-x0
      yy=yt-y0
      rr=xx*xx+yy*yy
      if (dabs(rr)<1.0d-16) cycle
      rr2=dsqrt(rr)
      !ff=dbsk1(gsareaconst%k_helm*rr2)*(xx*sintt-yy*costt)/rr2
      ff=dbsk1(gsareaconst%k_helm*rr2)*eta0/rr2
      s=s+ff
    enddo
    s=s*gsareaconst%k_helm
  endif
case (16)
  if (gs_test_point_collenear_for_BesselG.and.collenear3z(b%z(k),b%z(k+1),z0).and.dabs(costt)<1.0d-8) then
    is_analit=.true.
  else
    if (gs_test_point_collenear_for_BesselG.and.dabs(sintt)<1.0d-8) then
      s=(dbsk0(gsareaconst%k_helm*cdabs(b%z(k+1)-z0))-dbsk0(gsareaconst%k_helm*cdabs(b%z(k)-z0)))*costt
      is_analit=.true.
    else
      do i=1,nint
        it=i-d5
        xt=b%x(k)+dx*it
        yt=b%y(k)+dy*it
        xx=xt-x0
        yy=yt-y0
        rr=xx*xx+yy*yy
        if (dabs(rr)<1.0d-16) cycle
        rr2=dsqrt(rr)
        ff=dbsk1(gsareaconst%k_helm*rr2)*xx/rr2
        s=s+ff
      enddo
      s=-s*gsareaconst%k_helm
    endif
  endif
case (17)
  if (gs_test_point_collenear_for_BesselG.and.collenear3z(b%z(k),b%z(k+1),z0).and.dabs(sintt)<1.0d-8) then
    is_analit=.true.
  else
    if (gs_test_point_collenear_for_BesselG.and.dabs(costt)<1.0d-8) then
      s=(dbsk0(gsareaconst%k_helm*cdabs(b%z(k+1)-z0))-dbsk0(gsareaconst%k_helm*cdabs(b%z(k)-z0)))*sintt
      is_analit=.true.
    else
      do i=1,nint
        it=i-d5
        xt=b%x(k)+dx*it
        yt=b%y(k)+dy*it
        xx=xt-x0
        yy=yt-y0
        rr=xx*xx+yy*yy
        if (dabs(rr)<1.0d-16) cycle
        rr2=dsqrt(rr)
        ff=dbsk1(gsareaconst%k_helm*rr2)*yy/rr2
        s=s+ff
      enddo
      s=-s*gsareaconst%k_helm
    endif
  endif
case (18)
  if (gs_test_point_collenear_for_BesselG.and.collenear3z(b%z(k),b%z(k+1),z0)) then
    if (dabs(sintt)<1.0d-8) then
      is_analit=.true.
    else
      do i=1,nint
        it=i-d5
        xt=b%x(k)+dx*it
        yt=b%y(k)+dy*it
        xx=xt-x0
        yy=yt-y0
        rr=xx*xx+yy*yy
        if (dabs(rr)<1.0d-16) cycle
        rr2=dsqrt(rr)
        ff=dbsk1(gsareaconst%k_helm*rr2)/rr2
        s=s+ff
      enddo
      s=-s*gsareaconst%k_helm*sintt
    endif
  else
    do i=1,nint
      it=i-d5
      xt=b%x(k)+dx*it
      yt=b%y(k)+dy*it
      xx=xt-x0
      yy=yt-y0
      rr=xx*xx+yy*yy
      if (dabs(rr)<1.0d-16) cycle
      rr2=dsqrt(rr)
	    ff2=dbsk1(gsareaconst%k_helm*rr2)/rr2
      !ff=(gsareaconst%k_helm*dbsk0(gsareaconst%k_helm*rr2)+d2*ff2)*xx/rr*(xx*sintt-yy*costt)-ff2*sintt
      ff=(gsareaconst%k_helm*dbsk0(gsareaconst%k_helm*rr2)+d2*ff2)*xx/rr*eta0-ff2*sintt
      s=s+ff
    enddo
    s=s*gsareaconst%k_helm
  endif
case (19)
  if (gs_test_point_collenear_for_BesselG.and.collenear3z(b%z(k),b%z(k+1),z0)) then
    if (dabs(costt)<1.0d-8) then
      is_analit=.true.
    else
      do i=1,nint
        it=i-d5
        xt=b%x(k)+dx*it
        yt=b%y(k)+dy*it
        xx=xt-x0
        yy=yt-y0
        rr=xx*xx+yy*yy
        if (dabs(rr)<1.0d-16) cycle
        rr2=dsqrt(rr)
        ff=dbsk1(gsareaconst%k_helm*rr2)/rr2
        s=s+ff
      enddo
      s=s*gsareaconst%k_helm*costt
    endif
  else
    do i=1,nint
      it=i-d5
      xt=b%x(k)+dx*it
      yt=b%y(k)+dy*it
      xx=xt-x0
      yy=yt-y0
      rr=xx*xx+yy*yy
      if (dabs(rr)<1.0d-16) cycle
      rr2=dsqrt(rr)
	    ff2=dbsk1(gsareaconst%k_helm*rr2)/rr2
      !ff=(gsareaconst%k_helm*dbsk0(gsareaconst%k_helm*rr2)+d2*ff2)*yy/rr*(xx*sintt-yy*costt)+ff2*costt
      ff=(gsareaconst%k_helm*dbsk0(gsareaconst%k_helm*rr2)+d2*ff2)*yy/rr*eta0+ff2*costt
      s=s+ff
    enddo
    s=s*gsareaconst%k_helm
  endif
case default
  call gs_print_stop("Error calcintn!")
end select
if (is_analit) then
  calcintn=s
else
  calcintn=s*b%l(k)/nint
endif
end

subroutine calcintn_dual(k,x0,y0,nf,ia,knd,s,s2)
!численное вычисление интегралов
!nf 
!   16- dG_4/dx_0   
!   17- dG_4/dy_0        
!   18- dG_4'/dx_0
!   19- dG_4'/dy_0
use pgmod
integer(4) k,i,nf,knd,nint,ia
real(8) x0,y0,dx,dy,xx,yy,rr,rr2,xt,yt,s,s2,ff,costt,sintt,ss,ff2,ff3
real(8) ds0,it,eta0 !,eta
complex(8) z0
logical collenear3z,is_analit,is_analit2
type(TBound), pointer :: b
b=>gs%a(ia)%bnd(knd)
costt=dreal(b%ett(k))
sintt=dimag(b%ett(k))
!ds0=1.0d-3
ds0=gs%const%ds0numint
nint=b%l(k)/ds0
if (mod(nint,2).eq.1) nint=nint+1
dx=(b%x(k+1)-b%x(k))/nint
dy=(b%y(k+1)-b%y(k))/nint
s=d0
s2=d0
z0=dcmplx(x0,y0)
is_analit=.false.
is_analit2=.false.
select case (nf)
case (16)
  if (gs_test_point_collenear_for_BesselG) then
    if (collenear3z(b%z(k),b%z(k+1),z0)) then
      is_analit=dabs(costt)<1.0d-8
      is_analit2=dabs(sintt)<1.0d-8
    else
      if (dabs(sintt)<1.0d-8) then
        s=(dbsk0(gsareaconst%k_helm*cdabs(b%z(k+1)-z0))-dbsk0(gsareaconst%k_helm*cdabs(b%z(k)-z0)))*costt
        is_analit=.true.
      endif
      if (dabs(costt)<1.0d-8) then
        s2=(dbsk0(gsareaconst%k_helm*cdabs(b%z(k+1)-z0))-dbsk0(gsareaconst%k_helm*cdabs(b%z(k)-z0)))*sintt
        is_analit2=.true.
      endif
    endif
  endif
  if ((.not.is_analit).or.(.not.is_analit2)) then
    do i=1,nint
      it=i-d5
      xt=b%x(k)+dx*it
      yt=b%y(k)+dy*it
      xx=xt-x0
      yy=yt-y0
      rr=xx*xx+yy*yy
      if (dabs(rr)<1.0d-16) cycle
      rr2=dsqrt(rr)
      ff=dbsk1(gsareaconst%k_helm*rr2)/rr2
      if (.not.is_analit) s=s+ff*xx
      if (.not.is_analit2) s2=s2+ff*yy
    enddo
    if (.not.is_analit) s=-s*gsareaconst%k_helm
    if (.not.is_analit2) s2=-s2*gsareaconst%k_helm
  endif
case (18)
  if (gs_test_point_collenear_for_BesselG.and.collenear3z(b%z(k),b%z(k+1),z0)) then
    ss=d0
    do i=1,nint
      it=i-d5
      xt=b%x(k)+dx*it
      yt=b%y(k)+dy*it
      xx=xt-x0
      yy=yt-y0
      rr=xx*xx+yy*yy
      if (dabs(rr)<1.0d-16) cycle
      rr2=dsqrt(rr)
      ff=dbsk1(gsareaconst%k_helm*rr2)/rr2
      ss=ss+ff
    enddo
    if (dabs(sintt)<1.0d-8) then
      is_analit=.true.
    else
      s=-ss*gsareaconst%k_helm*sintt
    endif
    if (dabs(costt)<1.0d-8) then
      is_analit2=.true.
    else
      s2=ss*gsareaconst%k_helm*costt
    endif
  else
    eta0=dimag((z0-b%z(k))*dconjg(b%ett(k)))
    do i=1,nint
      it=i-d5
      xt=b%x(k)+dx*it
      yt=b%y(k)+dy*it
      xx=xt-x0
      yy=yt-y0
      rr=xx*xx+yy*yy
      if (dabs(rr)<1.0d-16) cycle
      rr2=dsqrt(rr)
	    ff2=dbsk1(gsareaconst%k_helm*rr2)/rr2
      ff3=(gsareaconst%k_helm*dbsk0(gsareaconst%k_helm*rr2)+d2*ff2)/rr*eta0
      s=s+ff3*xx-ff2*sintt
      s2=s2+ff3*yy+ff2*costt
    enddo
    s=s*gsareaconst%k_helm
    s2=s2*gsareaconst%k_helm
  endif
case default
  call gs_print_stop("Error calcintn_dual!")
end select
if (.not.is_analit) s=s*b%l(k)/nint
if (.not.is_analit2) s2=s2*b%l(k)/nint
end

function calcintn2(k,i,nf,ia,knd,knd0)
!интегралы, допускающие только численное вычисление
use pgmod
real(8) calcintn2
real(8) calcintn,val,get_bb_cash_integral
integer(4) k,nf,knd,i,knd0,ia
type(TBound), pointer :: b0
val=real8_inf
if (gs_use_cash) then
  if (gs_cash_presolve) then
    calcintn2=get_bb_cash_integral_2(gs%a(ia),knd,k,knd0,i,nf)
    !DEC$ IF DEFINED (DEBUG)
    if (calcintn2==real8_inf) call gs_print_stop('Error calcintn2 gs_cash_presolve!!!')
    !DEC$ ENDIF
    return
  endif
  val=get_bb_cash_integral(ia,knd,k,knd0,i,nf)
endif
if (val==real8_inf) then
  b0=>gs%a(ia)%bnd(knd0)
  val=calcintn(k,b0%xc(i),b0%yc(i),nf,ia,knd)
  if (gs_use_cash) call set_bb_cash_integral(ia,knd,k,knd0,i,nf,val)
endif
calcintn2=val
end

function calcintn_ar(k,i,nf,ia,knd)
use pgmod
real(8) calcintn_ar
real(8) calcintn,xmc,ymc,get_ba_cash_integral,val
integer(4) k,nf,knd,i,ia
val=real8_inf
if (gs_use_cash) then
  if (gs_cash_presolve) then
    calcintn_ar=get_ba_cash_integral_2(gs%a(ia),knd,k,i,nf)
    !DEC$ IF DEFINED (DEBUG)
    if (calcintn_ar==real8_inf) call gs_print_stop('Error calcintn_ar gs_cash_presolve!!!')
    !DEC$ ENDIF
    return
  endif
  val=get_ba_cash_integral(ia,knd,k,i,nf)
endif
if (val==real8_inf) then
  val=calcintn(k,xmc(i,ia),ymc(i,ia),nf,ia,knd)
  if (gs_use_cash) call set_ba_cash_integral(ia,knd,k,i,nf,val)
endif
calcintn_ar=val
end

function calcint(k,x0,y0,nf,ia,knd)
!вычисление интеграла в области
use pgmod
real(8) calcint,x0,y0,calcinta
real(8) calcintn
integer(4) k,nf,knd,ia
if (gs%const%use_numerical_int) then
  calcint=calcintn(k,x0,y0,nf,ia,knd)
else
  calcint=calcinta(k,x0,y0,nf,ia,knd)
endif
end

function calcinta(k,x0,y0,nf,ia,knd)
!аналитическое вычисление интегралов
use pgmod
real(8) calcinta,x0,y0,calcint_nf0,calcint_nf1,calcint_nf2,calcint_nf3,calcint_nf45,calcint_nf67,calcint_nf89,calcint_nf1011
integer(4) k,nf,knd,ia
complex(8) z0,z1,z2
real(8),PARAMETER :: eps=1.0d-4
type(TBound), pointer :: b
b=>gs%a(ia)%bnd(knd)
z0=dcmplx(x0,y0)
z1=b%z(k)
z2=b%z(k+1)
select case (nf)
case (0)
  calcinta=calcint_nf0(z0,z1,z2)
case (1)
  calcinta=calcint_nf1(z0,z1,z2)
case (2)
  calcinta=calcint_nf2(z0,z1,z2)
case (3)
  calcinta=calcint_nf3(z0,z1,z2)
case (4,5)
  calcinta=calcint_nf45(z0,z1,z2,nf)
case (6,7)
  calcinta=calcint_nf67(z0,z1,z2,nf)
case (8,9)
  calcinta=calcint_nf89(z0,z1,z2,nf)
case (10,11)
  calcinta=calcint_nf1011(z0,z1,z2,nf)
case default
  call gs_print_stop("Error calcinta!")
endselect
end

subroutine calcinta_2(k,x0,y0,nf,ia,knd,val,val2)
!аналитическое вычисление интегралов
use pgmod
real(8) val,val2,x0,y0
integer(4) k,nf,knd,ia
complex(8) z0,z1,z2
real(8),PARAMETER :: eps=1.0d-4
type(TBound), pointer :: b
b=>gs%a(ia)%bnd(knd)
z0=dcmplx(x0,y0)
z1=b%z(k)
z2=b%z(k+1)
select case (nf)
case (4)
  call calcint_nf45_2(z0,z1,z2,val,val2)
case (6)
  call calcint_nf67_2(z0,z1,z2,val,val2)
case (8)
  call calcint_nf89_2(z0,z1,z2,val,val2)
case (10)
  call calcint_nf1011_2(z0,z1,z2,val,val2)
case default
  call gs_print_stop("Error calcinta_2!")
endselect
end

function calcint_nf0(z0,z1,z2)
!аналитическое вычисление интеграла от G_1
use gen_mod
real(8) calcint_nf0,delta_n,delta_k,tt
complex(8) z1,z2,z0,zn,zk,fz,lnn,lnk
call get_delta_l(z0,z1,z2,delta_n,delta_k)
zn=z1-z0
zk=z2-z0
tt=zarg(z2-z1)
lnn=dlog(cdabs(zn))+ii*delta_n
lnk=dlog(cdabs(zk))+ii*delta_k
fz=zk*(lnk-c1)-zn*(lnn-c1)
calcint_nf0=dreal(cdexp(-ii*tt)*fz)
end

function calcint_nf1(z0,z1,z2)
!аналитическое вычисление интеграла от G'_1
use gen_mod
real(8) calcint_nf1,delta_n,delta_k
complex(8) z1,z2,z0,zn,zk,fz
call get_delta_l(z0,z1,z2,delta_n,delta_k)
zn=z1-z0
zk=z2-z0
fz=dlog(cdabs(zk/zn))+ii*(delta_k-delta_n)
calcint_nf1=dimag(fz)
end

function calcint_nf2(z0,z1,z2)
!аналитическое вычисление интеграла от G_2
use gen_mod
real(8) calcint_nf2,delta_n,delta_k,tt
complex(8) z1,z2,z0,zn,zk,fz1k,fz1n,lnk,lnn,fz2k,fz2n,fz
call get_delta_l(z0,z1,z2,delta_n,delta_k)
zn=z1-z0
zk=z2-z0
tt=zarg(z2-z1)
lnn=dlog(cdabs(zn))+ii*delta_n
lnk=dlog(cdabs(zk))+ii*delta_k
fz1n=zn**2*(d2*lnn-3.0d0)*dconjg(zn)
fz1k=zk**2*(d2*lnk-3.0d0)*dconjg(zk)
fz2n=zn**3*(6.0d0*lnn-11.0d0)
fz2k=zk**3*(6.0d0*lnk-11.0d0)
fz=fz1k-fz1n-(fz2k-fz2n)*cdexp(-d2*ii*tt)/9.0d0
calcint_nf2=dreal(cdexp(-ii*tt)*fz/16.0d0)
end

function calcint_nf3(z0,z1,z2)
!аналитическое вычисление интеграла от G'_2
use gen_mod
real(8) calcint_nf3,delta_n,delta_k,tt
complex(8) z1,z2,z0,zn,zk,fz,lnk,lnn,fz1k,fz1n,fz2k,fz2n
call get_delta_l(z0,z1,z2,delta_n,delta_k)
zn=z1-z0
zk=z2-z0
tt=zarg(z2-z1)
lnn=dlog(cdabs(zn))+ii*delta_n
lnk=dlog(cdabs(zk))+ii*delta_k
fz1n=zn*dconjg(zn)*(lnn-d1)
fz1k=zk*dconjg(zk)*(lnk-d1)
fz2n=zn**2*(d2*lnn-3.0d0)
fz2k=zk**2*(d2*lnk-3.0d0)
fz=(fz1k-fz1n)-d5*(fz2k-fz2n)*cdexp(-d2*ii*tt)
calcint_nf3=dimag(0.25d0*fz)
end

function calcint_nf45z(z0,z1,z2)
!аналитическое вычисление интеграла от dG_1/dx_0, dG_1/dy_0
use gen_mod
complex(8) calcint_nf45z
real(8) delta_n,delta_k,tt
complex(8) z1,z2,z0,zn,zk
call get_delta_l(z0,z1,z2,delta_n,delta_k)
zn=z1-z0
zk=z2-z0
tt=zarg(z2-z1)
calcint_nf45z=(dlog(cdabs(zk/zn))+ii*(delta_k-delta_n))*cdexp(-ii*tt)
end

function calcint_nf45(z0,z1,z2,nf)
!аналитическое вычисление интеграла от dG_1/dx_0, dG_1/dy_0
use gen_mod
integer(4) nf
real(8) calcint_nf45
complex(8) z1,z2,z0,fz,calcint_nf45z
fz=calcint_nf45z(z0,z1,z2)
if (nf==4) then
  calcint_nf45=-dreal(fz)
else
  calcint_nf45=dimag(fz)
endif
end

subroutine calcint_nf45_2(z0,z1,z2,val,val2)
!аналитическое вычисление интеграла от dG_1/dx_0, dG_1/dy_0
use gen_mod
real(8) val,val2
complex(8) z1,z2,z0,fz,calcint_nf45z
fz=calcint_nf45z(z0,z1,z2)
val=-dreal(fz)
val2=dimag(fz)
end

function calcint_nf67z(z0,z1,z2)
!аналитическое вычисление интеграла от dG'_1/dx_0, dG'_1/dy_0
use gen_mod
complex(8) calcint_nf67z,z1,z2,z0,zn,zk
zn=z1-z0
zk=z2-z0
calcint_nf67z=c1/zk-c1/zn
end

function calcint_nf67(z0,z1,z2,nf)
!аналитическое вычисление интеграла от dG'_1/dx_0, dG'_1/dy_0
use gen_mod
integer(4) nf
real(8) calcint_nf67
complex(8) z1,z2,z0,fz,calcint_nf67z
fz=calcint_nf67z(z0,z1,z2)
if (nf==6) then
  calcint_nf67=-dimag(fz)
else
  calcint_nf67=-dreal(fz)
endif
end

subroutine calcint_nf67_2(z0,z1,z2,val,val2)
!аналитическое вычисление интеграла от dG'_1/dx_0, dG'_1/dy_0
use gen_mod
real(8) val,val2
complex(8) z1,z2,z0,fz,calcint_nf67z
fz=calcint_nf67z(z0,z1,z2)
val=-dimag(fz)
val2=-dreal(fz)
end

subroutine calcint_nf89z(z0,z1,z2,fz1,fz2,fz3)
!аналитическое вычисление интеграла от dG_2/dx_0, dG_2/dy_0
use gen_mod
real(8) delta_n,delta_k,tt
complex(8) z1,z2,z0,zn,zk,fz1k,fz1n,lnk,lnn,fz2k,fz2n,fz1,fz2,fz3
call get_delta_l(z0,z1,z2,delta_n,delta_k)
zn=z1-z0
zk=z2-z0
tt=zarg(z2-z1)
lnn=dlog(cdabs(zn))+ii*delta_n
lnk=dlog(cdabs(zk))+ii*delta_k
fz1n=zn*(lnn-d1)*dconjg(zn)
fz1k=zk*(lnk-d1)*dconjg(zk)
fz2n=zn**2*(d2*lnn-3.0d0)
fz2k=zk**2*(d2*lnk-3.0d0)
fz1=-4.0d0*(fz1k-fz1n)
fz2=fz2k-fz2n
fz3=cdexp(-ii*tt)
!if (nf==8) then
!  fz=-4.0d0*(fz1k-fz1n)+(fz2k-fz2n)*(cdexp(-d2*ii*tt)-d1)
!  calcint_nf89=dreal(cdexp(-ii*tt)*fz/16.0d0)
!else
!  fz=-4.0d0*(fz1k-fz1n)+(fz2k-fz2n)*(cdexp(-d2*ii*tt)+d1)
!  calcint_nf89=-dimag(cdexp(-ii*tt)*fz/16.0d0)
!endif
end

function calcint_nf89(z0,z1,z2,nf)
!аналитическое вычисление интеграла от dG_2/dx_0, dG_2/dy_0
use gen_mod
integer(4) nf
real(8) calcint_nf89
complex(8) z1,z2,z0,fz,fz1,fz2,fz3
call calcint_nf89z(z0,z1,z2,fz1,fz2,fz3)
if (nf==8) then
  fz=fz1+fz2*(fz3**2-c1)
  calcint_nf89=dreal(fz3*fz/16.0d0)
else
  fz=fz1+fz2*(fz3**2+c1)
  calcint_nf89=-dimag(fz3*fz/16.0d0)
endif
end

subroutine calcint_nf89_2(z0,z1,z2,val,val2)
!аналитическое вычисление интеграла от dG_2/dx_0, dG_2/dy_0
use gen_mod
real(8) val,val2
complex(8) z1,z2,z0,fz,fz1,fz2,fz3
call calcint_nf89z(z0,z1,z2,fz1,fz2,fz3)
fz=fz1+fz2*(fz3**2-c1)
val=dreal(fz3*fz/16.0d0)
fz=fz1+fz2*(fz3**2+c1)
val2=-dimag(fz3*fz/16.0d0)
end

subroutine calcint_nf1011z(z0,z1,z2,fz1,fz2,fz3)
!аналитическое вычисление интеграла от dG'_2/dx_0, dG'_2/dy_0
use gen_mod
real(8) delta_n,delta_k,tt
complex(8) z1,z2,z0,zn,zk,fz1k,fz1n,lnk,lnn,fz2k,fz2n,fz1,fz2,fz3
call get_delta_l(z0,z1,z2,delta_n,delta_k)
zn=z1-z0
zk=z2-z0
tt=zarg(z2-z1)
lnn=dlog(cdabs(zn))+ii*delta_n
lnk=dlog(cdabs(zk))+ii*delta_k
fz1n=lnn*dconjg(zn)
fz1k=lnk*dconjg(zk)
fz2n=zn*(lnn-d1)
fz2k=zk*(lnk-d1)
fz1=-(fz1k-fz1n)
fz2=fz2k-fz2n
fz3=d2*cdexp(-d2*ii*tt)
!if (nf==10) then
!  fz=-(fz1k-fz1n)+(fz2k-fz2n)*(d2*cdexp(-d2*ii*tt)-d1)
!  calcint_nf1011=dimag(fz)*0.25d0
!else
!  fz=-(fz1k-fz1n)+(fz2k-fz2n)*(d2*cdexp(-d2*ii*tt)+d1)
!  calcint_nf1011=dreal(fz)*0.25d0
!endif
end

function calcint_nf1011(z0,z1,z2,nf)
!аналитическое вычисление интеграла от dG'_2/dx_0, dG'_2/dy_0
use gen_mod
integer(4) nf
real(8) calcint_nf1011
complex(8) z1,z2,z0,fz,fz1,fz2,fz3
call calcint_nf1011z(z0,z1,z2,fz1,fz2,fz3)
if (nf==10) then
  fz=fz1+fz2*(fz3-c1)
  calcint_nf1011=dimag(fz)*0.25d0
else
  fz=fz1+fz2*(fz3+c1)
  calcint_nf1011=dreal(fz)*0.25d0
endif
end

subroutine calcint_nf1011_2(z0,z1,z2,val,val2)
!аналитическое вычисление интеграла от dG'_2/dx_0, dG'_2/dy_0
use gen_mod
real(8) val,val2
complex(8) z1,z2,z0,fz,fz1,fz2,fz3
call calcint_nf1011z(z0,z1,z2,fz1,fz2,fz3)
fz=fz1+fz2*(fz3-c1)
val=dimag(fz)*0.25d0
fz=fz1+fz2*(fz3+c1)
val2=dreal(fz)*0.25d0
end

subroutine get_delta_l(z,tn,tk,delta_n,delta_k)
use gen_mod
complex(8) z,tn,tk,a,b,c
real(8) eps,vp,delta_n,delta_k,ddelta !,sp
eps=1.0d-8
a=tk-tn
b=tn-z
c=tk-z
delta_n=zarg(b)
delta_k=zarg(c)
vp=dimag(dconjg(a)*b)
if (cdabs(vp/a)<eps) then
  !на одной линии
  !sp=dreal(dconjg(b)*c)
  !if (sp<d0) then
  !  ddelta=-pi
  !else
  !  ddelta=d0
  !endif
  ddelta=d0
else
  ddelta=delta_k-delta_n
  if (vp<d0) then
    !z слева
	if (ddelta<d0) ddelta=ddelta+pi2
  else
    !z справа
	if (ddelta>d0) ddelta=ddelta-pi2
  endif
endif
delta_k=delta_n+ddelta
end

function calcintm(j,i,nf,ia,knd,knd0)
!интеграл, когда контрольная точка является серединой панели
!интегралы, допускающие аналитическое вычисление
use pgmod
!j - номер тек панели, по которой идет интегрирование
!i - номер контрольной точки (связан с номером уравнения)
!nf - номер интеграла
!ia0, knd0 - номера области и границы с контрольной точкой
!ia, knd - номера области и границы, где находится панель интегрирования
integer(4) j,i,nf,knd,knd0,ia
real(8) calcintm,calcintm_,val,get_bb_cash_integral
val=real8_inf
if (gs_use_cash) then
  if (gs_cash_presolve) then
    calcintm=get_bb_cash_integral_2(gs%a(ia),knd,j,knd0,i,nf)
    !DEC$ IF DEFINED (DEBUG)
    if (calcintm==real8_inf) call gs_print_stop('Error calcintm gs_cash_presolve!!!')
    !DEC$ ENDIF
    return
  endif
  val=get_bb_cash_integral(ia,knd,j,knd0,i,nf)
endif
if (val==real8_inf) then
  val=calcintm_(j,i,nf,ia,knd,knd0)
  if (gs_use_cash) call set_bb_cash_integral(ia,knd,j,knd0,i,nf,val)
endif
calcintm=val
end

function calcintm_(j,i,nf,ia,knd,knd0)
use pgmod
integer(4) j,i,nf,knd,knd0,ia
real(8) calcintn,calcinta,calcintm_
type(TBound), pointer :: b0
b0=>gs%a(ia)%bnd(knd0)
if (gs%const%use_numerical_int) then
  calcintm_=calcintn(j,b0%xc(i),b0%yc(i),nf,ia,knd)
else
  calcintm_=calcinta(j,b0%xc(i),b0%yc(i),nf,ia,knd)
endif
end

subroutine calcintm_2(j,i,nf,ia,knd,knd0,val,val2)
use pgmod
integer(4) j,i,nf,knd,knd0,ia
real(8) val,val2
type(TBound), pointer :: b0
b0=>gs%a(ia)%bnd(knd0)
call calcinta_2(j,b0%xc(i),b0%yc(i),nf,ia,knd,val,val2)
end

function calcint_ar(j,i,nf,ia,knd)
!интеграл по панели, когда контрольная точка является серединой треугольника
use pgmod
!j - номер тек панели, по которой идет интегрирование
!i - номер контрольной точки - треугольника
!nf - номер интеграла
!ia, knd - номера области и границы, где находится панель интегрирования
integer(4) j,i,nf,knd,ia
real(8) calcint_ar,get_ba_cash_integral,val,calcint_ar_
val=real8_inf
if (gs_use_cash) then
  if (gs_cash_presolve) then
    calcint_ar=get_ba_cash_integral_2(gs%a(ia),knd,j,i,nf)
    !DEC$ IF DEFINED (DEBUG)
    if (calcint_ar==real8_inf) call gs_print_stop('Error calcint_ar gs_cash_presolve!!!')
    !DEC$ ENDIF
    return
  endif
  val=get_ba_cash_integral(ia,knd,j,i,nf)
endif
if (val==real8_inf) then
  val=calcint_ar_(j,i,nf,ia,knd)
  if (gs_use_cash) call set_ba_cash_integral(ia,knd,j,i,nf,val)
endif
calcint_ar=val
end

function calcint_ar_(j,i,nf,ia,knd)
use pgmod
integer(4) j,i,nf,knd,ia
real(8) calcint_ar_,calcinta,calcintn,xmc,ymc
if (gs%const%use_numerical_int) then
  calcint_ar_=calcintn(j,xmc(i,ia),ymc(i,ia),nf,ia,knd)
else
  calcint_ar_=calcinta(j,xmc(i,ia),ymc(i,ia),nf,ia,knd)
endif
end

subroutine calcint_ar_2(j,i,nf,ia,knd,val,val2)
use pgmod
integer(4) j,i,nf,knd,ia
real(8) xmc,ymc,val,val2
call calcinta_2(j,xmc(i,ia),ymc(i,ia),nf,ia,knd,val,val2)
end
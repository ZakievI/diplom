!интегралы по круговым панелям

function cp_calcint__(zc,r,z0,nf,g0,dir)
!численное вычисление интегралов gj кругу
!zc - центр круга
!r - радиус
!   0 - G_1
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
!   20 - G_1*sin(g)
!   21 - G_1*cos(g)
!   22 - G'_1*sin(g)
!   23 - G'_1*cos(g)
!   24 - G_2*sin(g)
!   25 - G_2*cos(g)
!   26 - G'_2*sin(g)
!   27 - G'_2*cos(g)
!   28 - d(G_1*sin(g))/dx_0
!   29 - d(G_1*sin(g))/dy_0
!   30 - d(G_1*cos(g))/dx_0
!   31 - d(G_1*cos(g))/dy_0
!   32 - d(G'_1*sin(g))/dx_0
!   33 - d(G'_1*sin(g))/dy_0
!   34 - d(G'_1*cos(g))/dx_0
!   35 - d(G'_1*cos(g))/dy_0
!   36 - d(G_2*sin(g))/dx_0
!   37 - d(G_2*sin(g))/dy_0
!   38 - d(G_2*cos(g))/dx_0
!   39 - d(G_2*cos(g))/dy_0
!   40 - d(G'_2*sin(g))/dx_0
!   41 - d(G'_2*sin(g))/dy_0
!   42 - d(G'_2*cos(g))/dx_0
!   43 - d(G'_2*cos(g))/dy_0
use gen_mod
integer(4) i,nint,nf
real(8) dg,s,g,ff,pi32,tt,sintt,costt,cp_calcint__
real(8) xt,yt,xx,yy,rr,xc,yc,r,x0,y0
real(8) g0,g1
integer(4) dir
complex(8) zc,z0,zt,fz1
xc=dreal(zc)
yc=dimag(zc)
x0=dreal(z0)
y0=dimag(z0)
nint=100
pi32=pi*1.5d0
s=d0
dg=pi2/nint
do i=1,nint
  g=(i-d5)*dg
  g1=g0+g*dir
  xt=xc+r*dcos(g1)
  yt=yc+r*dsin(g1)
  xx=xt-x0
  yy=yt-y0
  rr=xx*xx+yy*yy
  if (dabs(rr)<1.0d-16) then
    nint=nint
    cycle
  endif
  tt=g1+pi5*dir  !fi
  sintt=dsin(tt)
  costt=dcos(tt)
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
    fz1=cdexp(ii*tt)/(zt**2)
    ff=dimag(fz1)
    !ff=((xx**2-yy**2)*sintt-d2*xx*yy*costt)/rr**2
  elseif (nf.eq.7) then
    zt=dcmplx(xx,yy)
    fz1=cdexp(ii*tt)/(zt**2)
    ff=dreal(fz1)
    !ff=((xx**2-yy**2)*costt+d2*xx*yy*sintt)/rr**2
  elseif (nf.eq.8) then
    ff=-xx*(dlog(rr)-d1)*0.25d0
  elseif (nf.eq.9) then
    ff=-yy*(dlog(rr)-d1)*0.25d0
  elseif (nf.eq.10) then
    ff=(-sintt*(dlog(rr)-d1+d2*xx*xx/rr)+costt*d2*xx*yy/rr)*0.25d0
  elseif (nf.eq.11) then
    ff=(costt*(dlog(rr)-d1+d2*yy*yy/rr)-sintt*d2*xx*yy/rr)*0.25d0
  elseif (nf.eq.20) then
    ff=d5*dlog(rr)*dsin(g)
  elseif (nf.eq.21) then
    ff=d5*dlog(rr)*dcos(g)
  elseif (nf.eq.22) then
	  ff=(xx*sintt-yy*costt)/rr*dsin(g)
  elseif (nf.eq.23) then
	  ff=(xx*sintt-yy*costt)/rr*dcos(g)
  elseif (nf.eq.24) then
    ff=rr*(d5*dlog(rr)-d1)*dsin(g)*0.25d0
  elseif (nf.eq.25) then
    ff=rr*(d5*dlog(rr)-d1)*dcos(g)*0.25d0
  elseif (nf.eq.26) then
    ff=(xx*sintt-yy*costt)*(dlog(rr)-d1)*dsin(g)*0.25d0
  elseif (nf.eq.27) then
    ff=(xx*sintt-yy*costt)*(dlog(rr)-d1)*dcos(g)*0.25d0
  elseif (nf.eq.28) then
    ff=-xx/rr*sin(g)
  elseif (nf.eq.29) then
    ff=-yy/rr*sin(g)
  elseif (nf.eq.30) then
    ff=-xx/rr*cos(g)
  elseif (nf.eq.31) then
    ff=-yy/rr*cos(g)  
  elseif (nf.eq.32) then
    zt=dcmplx(xx,yy)
    fz1=cdexp(ii*tt)/(zt**2)
    ff=dimag(fz1)*dsin(g)
    !ff=((xx**2-yy**2)*sintt-d2*xx*yy*costt)/rr**2*dsin(g)
  elseif (nf.eq.33) then
    zt=dcmplx(xx,yy)
    fz1=cdexp(ii*tt)/(zt**2)
    ff=dreal(fz1)*dsin(g)
    !ff=((xx**2-yy**2)*costt+d2*xx*yy*sintt)/rr**2*dsin(g)
  elseif (nf.eq.34) then
    zt=dcmplx(xx,yy)
    fz1=cdexp(ii*tt)/(zt**2)
    ff=dimag(fz1)*dcos(g)
    !ff=((xx**2-yy**2)*sintt-d2*xx*yy*costt)/rr**2*dcos(g)
  elseif (nf.eq.35) then
    zt=dcmplx(xx,yy)
    fz1=cdexp(ii*tt)/(zt**2)
    ff=dreal(fz1)*dcos(g)
    !ff=((xx**2-yy**2)*costt+d2*xx*yy*sintt)/rr**2*dcos(g)
  elseif (nf.eq.36) then
    ff=-xx*(dlog(rr)-d1)*0.25d0*dsin(g)
  elseif (nf.eq.37) then
    ff=-yy*(dlog(rr)-d1)*0.25d0*dsin(g)
  elseif (nf.eq.38) then
    ff=-xx*(dlog(rr)-d1)*0.25d0*dcos(g)
  elseif (nf.eq.39) then
    ff=-yy*(dlog(rr)-d1)*0.25d0*dcos(g)
  elseif (nf.eq.40) then
    ff=(-sintt*(dlog(rr)-d1+d2*xx*xx/rr)+costt*d2*xx*yy/rr)*0.25d0*dsin(g)
  elseif (nf.eq.41) then
    ff=(costt*(dlog(rr)-d1+d2*yy*yy/rr)-sintt*d2*xx*yy/rr)*0.25d0*dsin(g)
  elseif (nf.eq.42) then
    ff=(-sintt*(dlog(rr)-d1+d2*xx*xx/rr)+costt*d2*xx*yy/rr)*0.25d0*dcos(g)
  elseif (nf.eq.43) then
    ff=(costt*(dlog(rr)-d1+d2*yy*yy/rr)-sintt*d2*xx*yy/rr)*0.25d0*dcos(g)
  else
    call gs_print_stop("Error cp_calcint__!")
  endif
  s=s+ff
enddo
cp_calcint__=s*dg*r
end

function cp_calcint_(ztc,r,nf,g0,dir)
!численное вычисление интегралов gj кругу
!ztc=zc-z0
!zc - центр круга
!z0 - точка, в которой строится функция
!r - радиус
use pgmod
integer(4) i,nint,nf
real(8) dg,s,g,ff,cp_calcint_
real(8) rr,r,r_,xy
complex(8) ztc,ksi,zeta,ett
real(8) g0,g1
integer(4) dir
nint=1000
s=d0
dg=pi2/nint
gs_use_parallel_num_int=1
!$omp parallel do if (gs_use_parallel_num_int==1) reduction(+ : s) Schedule(dynamic,1000) private (i,g,ff,rr,r_,xy,zeta,ksi,g1,ett) shared (nint,nf,ztc,r,g0,dir,dg) 
do i=1,nint
  g=(i-d5)*dg
  g1=g0+g*dir
  ett=cdexp(ii*g1)
  ksi=r*ett !радиус вектор точки, бегущей по кругу
  ett=ett*ii*dir
  zeta=ztc+ksi !=z1-z0
  r_=cdabs(zeta)
  if (r_/r<1.0d-8) cycle
  rr=r_*r_
  xy=dir*dreal(dconjg(zeta)*ksi/r) !=xx*sintt-yy*costt
  if (nf.eq.0) then
	  ff=dlog(r_)
  elseif (nf.eq.1) then
	  ff=xy/rr
  elseif (nf.eq.2) then
    ff=rr*(dlog(r_)-d1)*0.25d0
  elseif (nf.eq.3) then
    ff=xy*(dlog(rr)-d1)*0.25d0
  elseif (nf.eq.4) then
    ff=-dreal(zeta)/rr
  elseif (nf.eq.5) then
    ff=-dimag(zeta)/rr
  elseif (nf.eq.6) then
    ff=dimag(ett/zeta**2)
  elseif (nf.eq.7) then
    ff=dreal(ett/zeta**2)
  elseif (nf.eq.8) then
    ff=-dreal(zeta)*(dlog(rr)-d1)*0.25d0
  elseif (nf.eq.9) then
    ff=-dimag(zeta)*(dlog(rr)-d1)*0.25d0
  elseif (nf.eq.10) then
    ff=(-dimag(ett)*(dlog(rr)-d1+d2*dreal(zeta)**2/rr)+dreal(ett)*dimag(zeta**2)/rr)*0.25d0
  elseif (nf.eq.11) then
    ff=(dreal(ett)*(dlog(rr)-d1+d2*dimag(zeta)**2/rr)-dimag(ett)*dimag(zeta**2)/rr)*0.25d0
  elseif (nf.eq.20) then
    ff=dlog(r_)*dsin(g)
  elseif (nf.eq.21) then
    ff=dlog(r_)*dcos(g)
  elseif (nf.eq.22) then
	  ff=xy/rr*dsin(g)
  elseif (nf.eq.23) then
	  ff=xy/rr*dcos(g)
  elseif (nf.eq.24) then
    ff=rr*(dlog(r_)-d1)*dsin(g)*0.25d0
  elseif (nf.eq.25) then
    ff=rr*(dlog(r_)-d1)*dcos(g)*0.25d0
  elseif (nf.eq.26) then
    ff=xy*(dlog(rr)-d1)*dsin(g)*0.25d0
  elseif (nf.eq.27) then
    ff=xy*(dlog(rr)-d1)*dcos(g)*0.25d0
  elseif (nf.eq.28) then
    ff=-dreal(zeta)/rr*sin(g)
  elseif (nf.eq.29) then
    ff=-dimag(zeta)/rr*sin(g)
  elseif (nf.eq.30) then
    ff=-dreal(zeta)/rr*cos(g)
  elseif (nf.eq.31) then
    ff=-dimag(zeta)/rr*cos(g)  
  elseif (nf.eq.32) then
    ff=dimag(ett/zeta**2)*sin(g)
  elseif (nf.eq.33) then
    ff=dreal(ett/zeta**2)*sin(g)
  elseif (nf.eq.34) then
    ff=dimag(ett/zeta**2)*cos(g)
  elseif (nf.eq.35) then
    ff=dreal(ett/zeta**2)*cos(g)
  elseif (nf.eq.36) then
    ff=-dreal(zeta)*(dlog(rr)-d1)*0.25d0*dsin(g)
  elseif (nf.eq.37) then
    ff=-dimag(zeta)*(dlog(rr)-d1)*0.25d0*dsin(g)
  elseif (nf.eq.38) then
    ff=-dreal(zeta)*(dlog(rr)-d1)*0.25d0*dcos(g)
  elseif (nf.eq.39) then
    ff=-dimag(zeta)*(dlog(rr)-d1)*0.25d0*dcos(g)
  elseif (nf.eq.40) then
    ff=(-dimag(ett)*(dlog(rr)-d1+d2*dreal(zeta)**2/rr)+dreal(ett)*dimag(zeta**2)/rr)*0.25d0*dsin(g)
  elseif (nf.eq.41) then
    ff=(dreal(ett)*(dlog(rr)-d1+d2*dimag(zeta)**2/rr)-dimag(ett)*dimag(zeta**2)/rr)*0.25d0*dsin(g)
  elseif (nf.eq.42) then
    ff=(-dimag(ett)*(dlog(rr)-d1+d2*dreal(zeta)**2/rr)+dreal(ett)*dimag(zeta**2)/rr)*0.25d0*dcos(g)
  elseif (nf.eq.43) then
    ff=(dreal(ett)*(dlog(rr)-d1+d2*dimag(zeta)**2/rr)-dimag(ett)*dimag(zeta**2)/rr)*0.25d0*dcos(g)
  else
    call gs_print_stop("Error cp_calcint_!")
  endif
  s=s+ff
enddo
!omp end parallel do nowait
cp_calcint_=s*dg*r
end

function cp_calcint(ztc,r,nf,g0,dir,mode)
!численное вычисление интегралов gj кругу
!ztc=zc-z0
!zc - центр круга
!z0 - точка, в которой строится функция
!r - радиус
!mode
!0 - z0 снаружи круга
!1 - z0 на круге
!2 - z0 внутри круга
use pgmod
integer(4) nf
real(8) cp_calcint_,cp_calcint,cp_calcinta
real(8) r
complex(8) ztc
real(8) g0
integer(4) dir,mode
!DEC$ IF DEFINED (DEBUG) !!!сюда не должны приходить
if (.not.((nf>=0.and.nf<=11).or.(nf>=20.and.nf<=43))) call gs_print_stop("Wrong parameters for boundLineType=2 in cp_calcint. nf="//itoa(nf))
!DEC$ ENDIF
!if (.true.) then
if (gs%const%use_numerical_int) then
  cp_calcint=cp_calcint_(ztc,r,nf,g0,dir)
else
  cp_calcint=cp_calcinta(ztc,r,nf,g0,dir,mode)
endif
end

function cp_calcinta(ztc,r,nf,g0,dir,mode)
!аналитическое вычисление интегралов gj кругу
!ztc=zc-z0
!zc - центр круга
!z0 - точка, в которой строится функция
!r - радиус
!mode
!0 - z0 снаружи круга
!1 - z0 на круге
!2 - z0 внутри круга
use pgmod
integer(4) nf
!real(8) cp_calcint_
real(8) cp_calcinta,cp_calcint_nf1,cp_calcint_nf3,cp_calcint_nf2021,cp_calcint_nf2223,cp_calcint_nf0
real(8) r,cp_calcint_nf2425,cp_calcint_nf2627,cp_calcint_nf1011,cp_calcint_nf2831,cp_calcint_nf3235
real(8) cp_calcint_nf3639,cp_calcint_nf4043,cp_calcint_nf45,cp_calcint_nf2,cp_calcint_nf89
complex(8) ztc
real(8) g0
integer(4) dir,mode
select case (nf)
case (0)
  cp_calcinta=cp_calcint_nf0(ztc,r,mode)
case (1)
  cp_calcinta=cp_calcint_nf1(mode)
case (2)
  cp_calcinta=cp_calcint_nf2(ztc,r,mode)
case (3)
  cp_calcinta=cp_calcint_nf3(ztc,r,mode)
case (4,5)
  cp_calcinta=cp_calcint_nf45(ztc,r,mode,nf)
case (6,7)
  cp_calcinta=d0
case (8,9)
  cp_calcinta=cp_calcint_nf89(ztc,r,mode,nf)
case (10,11)
  cp_calcinta=cp_calcint_nf1011(ztc,r,mode,nf)
case (20,21)
  cp_calcinta=cp_calcint_nf2021(ztc,r,mode,nf)
case (22,23)
  cp_calcinta=cp_calcint_nf2223(ztc,r,mode,nf)
case (24,25)
  cp_calcinta=cp_calcint_nf2425(ztc,r,mode,nf)
case (26,27)
  cp_calcinta=cp_calcint_nf2627(ztc,r,mode,nf)
case (28:31)
  cp_calcinta=cp_calcint_nf2831(ztc,r,mode,nf)
case (32:35)
  cp_calcinta=cp_calcint_nf3235(ztc,r,mode,nf)
case (36:39)
  cp_calcinta=cp_calcint_nf3639(ztc,r,mode,nf)
case (40:43)
  cp_calcinta=cp_calcint_nf4043(ztc,r,mode,nf)
case default
  call gs_print_stop("Error cp_calcinta!")
  !cp_calcinta=cp_calcint_(ztc,r,nf,g0,dir)
end select
return
dir=dir
g0=g0
end

subroutine cp_calcinta_2(ztc,r,nf,g0,dir,mode,val,val2)
!аналитическое вычисление интегралов gj кругу
!ztc=zc-z0
!zc - центр круга
!z0 - точка, в которой строится функция
!r - радиус
!mode
!0 - z0 снаружи круга
!1 - z0 на круге
!2 - z0 внутри круга
use pgmod
integer(4) nf
!real(8) cp_calcint_
real(8) val,val2,r
complex(8) ztc
real(8) g0
integer(4) dir,mode
select case (nf)
case (4,5)
  call cp_calcint_nf45_2(ztc,r,mode,val,val2)
case (6,7)
  val=d0
  val2=d0
case (8,9)
  call cp_calcint_nf89_2(ztc,r,mode,val,val2)
case (10,11)
  call cp_calcint_nf1011_2(ztc,r,mode,val,val2)
case (20,21)
  call cp_calcint_nf2021_2(ztc,r,mode,val,val2)
case (22,23)
  call cp_calcint_nf2223_2(ztc,r,mode,val,val2)
case (24,25)
  call cp_calcint_nf2425_2(ztc,r,mode,val,val2)
case (26,27)
  call cp_calcint_nf2627_2(ztc,r,mode,val,val2)
case default
  call gs_print_stop("Error cp_calcinta_2!")
end select
return
dir=dir
g0=g0
end

function cp_calcint_nf0(ztc,r,mode)
!аналитическое вычисление интеграла от G_1
use gen_mod
integer(4) mode
real(8) cp_calcint_nf0
real(8) r,rr
complex(8) ztc
select case (mode)
case (0)     
  !z0 снаружи круга
  rr=cdabs(ztc)
  cp_calcint_nf0=pi2*r*dlog(rr)
case default
  !z0 на круге и внутри круга
  cp_calcint_nf0=pi2*r*dlog(r)
endselect
end

function cp_calcint_nf1(mode)
!аналитическое вычисление интеграла от G'_1
use gen_mod
integer(4) mode
real(8) cp_calcint_nf1
select case (mode)
case (0)     
  !z0 снаружи круга
  cp_calcint_nf1=d0
case (1)
  !z0 на круге
  cp_calcint_nf1=-pi
case default
  !z0 внутри круга
  cp_calcint_nf1=-pi2
endselect
end

function cp_calcint_nf2(ztc,r,mode) result (kk)
!аналитическое вычисление интеграла от G_2
use gen_mod
real(8) r,rr,r2,rr2
complex(8) ztc
integer(4) mode
real(8) kk
select case (mode)
case (0)
  !z0 снаружи круга
  rr=cdabs(ztc)
  r2=r**2
  kk=pi5*r*((rr**2+r2)*(dlog(rr)-d1)+r2)
case (1)
  !z0 на круге
  kk=pi5*r**3*(d2*dlog(r)-d1)
case default
  !z0 внутри круга
  rr=cdabs(ztc)
  rr2=rr**2
  kk=pi5*r*((rr2+r**2)*(dlog(r)-d1)+rr2)
endselect
end

function cp_calcint_nf3(ztc,r,mode) result (kk)
!аналитическое вычисление интеграла от G'_2
use gen_mod
real(8) r,rr
complex(8) ztc
integer(4) mode
real(8) kk
select case (mode)
case (0)
  !z0 снаружи круга
  rr=cdabs(ztc)
  kk=-pi*r**2*dlog(rr)
case (1)
  !z0 на круге
  kk=-pi*r**2*dlog(r)
case default
  !z0 внутри круга
  rr=cdabs(ztc)
  kk=-pi5*(rr**2+r**2*(d2*dlog(r)-d1))
endselect
end

function cp_calcint_nf45z(ztc,r,mode)
!аналитическое вычисление интеграла от
!nf=
!   4- dG_1/dx_0
!   5- dG_1/dy_0
use pgmod
real(8) r
complex(8) ztc,kk
integer(4) mode
complex(8) cp_calcint_nf45z
select case (mode)
case (0)   
  !z0 снаружи круга
  kk=-pi2*r/dconjg(ztc)
case (1) 
  !z0 на круге
  kk=c0
  !DEC$ IF DEFINED (DEBUG)  !!!сюда не должны приходить
  if (gs_stop_on_int_dxdy_inbound) call gs_print_stop("Wrong parameters for boundLineType=2 in cp_calcint")
  !DEC$ ENDIF
case default                   
  !z0 внутри круга
  kk=c0
endselect
cp_calcint_nf45z=kk
end

function cp_calcint_nf45(ztc,r,mode,nf)
!аналитическое вычисление интеграла от
!nf=
!   4- dG_1/dx_0
!   5- dG_1/dy_0
use pgmod
real(8) r
complex(8) ztc,cp_calcint_nf45z,kk
integer(4) mode,nf
real(8) cp_calcint_nf45
kk=cp_calcint_nf45z(ztc,r,mode)
if (nf==4) then
  cp_calcint_nf45=dreal(kk)
else
  cp_calcint_nf45=dimag(kk)
endif
end

subroutine cp_calcint_nf45_2(ztc,r,mode,val,val2)
!аналитическое вычисление интеграла от
!nf=
!   4- dG_1/dx_0
!   5- dG_1/dy_0
use pgmod
real(8) r,val,val2
complex(8) ztc,cp_calcint_nf45z,kk
integer(4) mode
kk=cp_calcint_nf45z(ztc,r,mode)
val=dreal(kk)
val2=dimag(kk)
end

function cp_calcint_nf2021z(ztc,r,mode)
!аналитическое вычисление интеграла от
!nf=
!   20 - G_1*sin(g)
!   21 - G_1*cos(g)
use gen_mod
real(8) r
complex(8) ztc,kk
integer(4) mode
complex(8) cp_calcint_nf2021z
if (mode==0) then   
  !z0 снаружи круга
  kk=pi*r**2/dconjg(ztc)
else                
  !z0 на круге или z0 внутри круга
  kk=pi*ztc
endif
cp_calcint_nf2021z=dconjg(kk)
end

function cp_calcint_nf2021(ztc,r,mode,nf)
!аналитическое вычисление интеграла от
!nf=
!   20 - G_1*sin(g)
!   21 - G_1*cos(g)
use gen_mod
real(8) r
complex(8) ztc,kk,cp_calcint_nf2021z
integer(4) mode,nf
real(8) cp_calcint_nf2021
kk=cp_calcint_nf2021z(ztc,r,mode)
if (nf==20) then
  cp_calcint_nf2021=dimag(kk)
else
  cp_calcint_nf2021=dreal(kk)
endif 
end

subroutine cp_calcint_nf2021_2(ztc,r,mode,val,val2)
!аналитическое вычисление интеграла от
!nf=
!   20 - G_1*sin(g)
!   21 - G_1*cos(g)
use gen_mod
real(8) r
complex(8) ztc,kk,cp_calcint_nf2021z
integer(4) mode
real(8) val,val2
kk=cp_calcint_nf2021z(ztc,r,mode)
val=dimag(kk)
val2=dreal(kk)
end

function cp_calcint_nf2021_(ztc,r,mode,nf)
!аналитическое вычисление интеграла от
!nf=
!   20 - G_1*sin(g)
!   21 - G_1*cos(g)
use gen_mod
real(8) r,rr,kk,gg
complex(8) ztc
integer(4) mode,nf
real(8) cp_calcint_nf2021,cp_calcint_nf2021_
rr=cdabs(ztc)
gg=zarg(ztc)+pi
if (mode==0) then   
  !z0 снаружи круга
  kk=-pi*r**2/rr
else                
  !z0 на круге или z0 внутри круга
  kk=-pi*rr
endif
if (nf==20) then
  cp_calcint_nf2021=-kk*dsin(gg)
else
  cp_calcint_nf2021=kk*dcos(gg)
endif
cp_calcint_nf2021_=cp_calcint_nf2021
end

function cp_calcint_nf2223z(ztc,r,mode)
!аналитическое вычисление интеграла от
!nf=
!   22 - G'_1*sin(g)
!   23 - G'_1*cos(g)
use gen_mod
real(8) r
complex(8) ztc,kk
integer(4) mode
complex(8) cp_calcint_nf2223z
select case (mode)
case (0)     
  !z0 снаружи круга
  kk=-pi*r/dconjg(ztc)
case (1) 
  !z0 на круге
  kk=0
case default
  !z0 внутри круга
  kk=pi*ztc/r
endselect
cp_calcint_nf2223z=dconjg(kk)
end

function cp_calcint_nf2223(ztc,r,mode,nf)
!аналитическое вычисление интеграла от
!nf=
!   22 - G'_1*sin(g)
!   23 - G'_1*cos(g)
use gen_mod
real(8) r
complex(8) ztc,kk
integer(4) mode,nf
real(8) cp_calcint_nf2223
complex(8) cp_calcint_nf2223z
kk=cp_calcint_nf2223z(ztc,r,mode)
if (nf==22) then
  cp_calcint_nf2223=dimag(kk)
else
  cp_calcint_nf2223=dreal(kk)
endif
end

subroutine cp_calcint_nf2223_2(ztc,r,mode,val,val2)
!аналитическое вычисление интеграла от
!nf=
!   22 - G'_1*sin(g)
!   23 - G'_1*cos(g)
use gen_mod
real(8) r
complex(8) ztc,kk
integer(4) mode
real(8) val,val2
complex(8) cp_calcint_nf2223z
kk=cp_calcint_nf2223z(ztc,r,mode)
val=dimag(kk)
val2=dreal(kk)
end

function cp_calcint_nf2223_(ztc,r,mode,nf) result (cp_calcint_nf2223)
!аналитическое вычисление интеграла от
!nf=
!   22 - G'_1*sin(g)
!   23 - G'_1*cos(g)
use gen_mod
real(8) r,rr,kk,gg
complex(8) ztc
integer(4) mode,nf
real(8) cp_calcint_nf2223
rr=cdabs(ztc)
gg=zarg(ztc)+pi
select case (mode)
case (0)
  !z0 снаружи круга
  kk=pi*r/rr
case (1)
  !z0 на круге
  kk=0
case default
  !z0 внутри круга
  kk=-pi*rr/r
endselect
if (nf==22) then
  cp_calcint_nf2223=-kk*dsin(gg)
else
  cp_calcint_nf2223=kk*dcos(gg)
endif
end

function cp_calcint_nf2425z(ztc,r,mode)
!аналитическое вычисление интеграла от
!nf=
!   24 - G_2*sin(g)
!   25 - G_2*cos(g)
use gen_mod
real(8) r,ro_c2
complex(8) ztc,kk,ztc1
integer(4) mode
complex(8) cp_calcint_nf2425z
select case (mode)
case (0)            
  !z0 снаружи круга
  ztc1=dconjg(ztc)
  ro_c2=dreal(ztc*ztc1)
  kk=pi*r**2*(d2*ro_c2*(dlog(ro_c2)-d1)+r**2)/(8.0d0*ztc1)
case (1)         
  !z0 на круге
  kk=-pi*ztc*r**2*(d1-4.0d0*dlog(r))/8.0d0
case default                         
  !z0 внутри круга
  ztc1=dconjg(ztc)
  ro_c2=dreal(ztc*ztc1)
  kk=pi*ztc*(ro_C2+d2*r**2*(d2*dlog(r)-d1))/8.0d0
endselect
cp_calcint_nf2425z=dconjg(kk)
end

function cp_calcint_nf2425(ztc,r,mode,nf)
!аналитическое вычисление интеграла от
!nf=
!   24 - G_2*sin(g)
!   25 - G_2*cos(g)
use gen_mod
real(8) r
complex(8) ztc,kk
integer(4) mode,nf
real(8) cp_calcint_nf2425
complex(8) cp_calcint_nf2425z
kk=cp_calcint_nf2425z(ztc,r,mode)
if (nf==24) then
  cp_calcint_nf2425=dimag(kk)
else
  cp_calcint_nf2425=dreal(kk)
endif
end

subroutine cp_calcint_nf2425_2(ztc,r,mode,val,val2)
!аналитическое вычисление интеграла от
!nf=
!   24 - G_2*sin(g)
!   25 - G_2*cos(g)
use gen_mod
real(8) r
complex(8) ztc,kk
integer(4) mode
real(8) val,val2
complex(8) cp_calcint_nf2425z
kk=cp_calcint_nf2425z(ztc,r,mode)
val=dimag(kk)
val2=dreal(kk)
end

function cp_calcint_nf2425_(ztc,r,mode,nf) result (cp_calcint_nf2425)
!аналитическое вычисление интеграла от
!nf=
!   24 - G_2*sin(g)
!   25 - G_2*cos(g)
use gen_mod
real(8) r,rr,kk,gg
complex(8) ztc
integer(4) mode,nf
real(8) cp_calcint_nf2425
rr=cdabs(ztc)
gg=zarg(ztc)+pi
select case (mode)
case (0)            
  !z0 снаружи круга
  kk=-pi*r**2*(d2*rr**2*(d2*dlog(rr)-d1)+r**2)/(8.0d0*rr)
case (1)          
  !z0 на круге
  kk=pi*r**3*(d1-4.0d0*dlog(r))/8.0d0
case default                          
  !z0 внутри круга
  kk=-pi*rr*(rr**2+d2*r**2*(d2*dlog(r)-d1))/8.0d0
endselect
if (nf==24) then
  cp_calcint_nf2425=-kk*dsin(gg)
else
  cp_calcint_nf2425=kk*dcos(gg)
endif
end

function cp_calcint_nf2627z(ztc,r,mode)
!аналитическое вычисление интеграла от
!nf=
!   26 - G'_2*sin(g)
!   27 - G'_2*cos(g)
use gen_mod
real(8) r,ro_c2
complex(8) ztc,kk,ztc1
integer(4) mode
complex(8) cp_calcint_nf2627z
ztc1=dconjg(ztc)
ro_c2=dreal(ztc*ztc1)
select case (mode)
case (0)      
  !z0 снаружи круга
  kk=-pi*r*(d2*ztc*(dlog(ro_c2)-d1)+3.0d0*r**2/ztc1)/8.0d0
case (1)   
  !z0 на круге
  kk=-pi*r*(d1+4.0d0*dlog(r))*ztc/8.0d0
case default                    
  !z0 внутри круга
  kk=pi*(ro_c2-d2*r**2*(d2*dlog(r)+d1))*ztc/(8.0d0*r)
endselect
cp_calcint_nf2627z=dconjg(kk)
end

function cp_calcint_nf2627(ztc,r,mode,nf)
!аналитическое вычисление интеграла от
!nf=
!   26 - G'_2*sin(g)
!   27 - G'_2*cos(g)
use gen_mod
real(8) r
complex(8) ztc,kk
integer(4) mode,nf
real(8) cp_calcint_nf2627
complex(8) cp_calcint_nf2627z
kk=cp_calcint_nf2627z(ztc,r,mode)
if (nf==26) then
  cp_calcint_nf2627=dimag(kk)
else
  cp_calcint_nf2627=dreal(kk)
endif
end

subroutine cp_calcint_nf2627_2(ztc,r,mode,val,val2)
!аналитическое вычисление интеграла от
!nf=
!   26 - G'_2*sin(g)
!   27 - G'_2*cos(g)
use gen_mod
real(8) r
complex(8) ztc,kk
integer(4) mode
real(8) val,val2
complex(8) cp_calcint_nf2627z
kk=cp_calcint_nf2627z(ztc,r,mode)
val=dimag(kk)
val2=dreal(kk)
end

function cp_calcint_nf2627_(ztc,r,mode,nf) result(cp_calcint_nf2627)
!аналитическое вычисление интеграла от
!nf=
!   26 - G'_2*sin(g)
!   27 - G'_2*cos(g)
use gen_mod
real(8) r,rr,kk,gg
complex(8) ztc
integer(4) mode,nf
real(8) cp_calcint_nf2627
rr=cdabs(ztc)
gg=zarg(ztc)+pi
select case (mode)
case (0)        
  !z0 снаружи круга
  kk=pi*r*(d2*rr**2*(d2*dlog(rr)-d1)+3.0d0*r**2)/(8.0d0*rr)
case (1)   
  !z0 на круге
  kk=pi*r**2*(d1+4.0d0*dlog(r))/8.0d0
case default                    
  !z0 внутри круга
  kk=-pi*rr*(rr**2-d2*r**2*(d2*dlog(r)+d1))/(8.0d0*r)
endselect
if (nf==26) then
  cp_calcint_nf2627=-kk*dsin(gg)
else
  cp_calcint_nf2627=kk*dcos(gg)
endif
end

function cp_calcint_nf89z(ztc,r,mode)
!аналитическое вычисление интеграла от
!nf=
!   8 - dG_2/dx_0
!   9 - dG_2/dy_0
use pgmod
real(8) r,ro_c2
complex(8) ztc,kk,ztc1
integer(4) mode
complex(8) cp_calcint_nf89z
ztc1=dconjg(ztc)
ro_c2=dreal(ztc*ztc1)
select case (mode)
case (0)      
  !z0 снаружи круга
  kk=-pi5*r*(ztc*(dlog(ro_c2)-d1)+r**2/ztc1)
case (1)   
  !z0 на круге
  kk=c0
  !DEC$ IF DEFINED (DEBUG)  !!!сюда не должны приходить
  if (gs_stop_on_int_dxdy_inbound) call gs_print_stop("Wrong parameters for boundLineType=2 in cp_calcint")
  !DEC$ ENDIF
case default                    
  !z0 внутри круга
  kk=-pi*r*ztc*dlog(r)
endselect
cp_calcint_nf89z=kk
end

function cp_calcint_nf89(ztc,r,mode,nf)
!аналитическое вычисление интеграла от
!nf=
!   8 - dG_2/dx_0
!   9 - dG_2/dy_0
use pgmod
real(8) r
complex(8) ztc,kk,cp_calcint_nf89z
integer(4) mode,nf
real(8) cp_calcint_nf89
kk=cp_calcint_nf89z(ztc,r,mode)
if (nf==8) then
  cp_calcint_nf89=dreal(kk)
else
  cp_calcint_nf89=dimag(kk)
endif
end

subroutine cp_calcint_nf89_2(ztc,r,mode,val,val2)
!аналитическое вычисление интеграла от
!nf=
!   8 - dG_2/dx_0
!   9 - dG_2/dy_0
use pgmod
real(8) r
complex(8) ztc,kk,cp_calcint_nf89z
integer(4) mode
real(8) val,val2
kk=cp_calcint_nf89z(ztc,r,mode)
val=dreal(kk)
val2=dimag(kk)
end

function cp_calcint_nf1011z(ztc,r,mode)
!аналитическое вычисление интеграла от
!nf=
!   10- dG'_2/dx_0
!   11- dG'_2/dy_0
use pgmod
real(8) r
complex(8) ztc,kk
integer(4) mode
complex(8) cp_calcint_nf1011z
select case (mode)
case (0)   
  !z0 снаружи круга
  kk=pi*r**2/dconjg(ztc)
case (1) 
  !z0 на круге
  kk=c0
  !DEC$ IF DEFINED (DEBUG)  !!!сюда не должны приходить
  if (gs_stop_on_int_dxdy_inbound) call gs_print_stop("Wrong parameters for boundLineType=2 in cp_calcint")
  !DEC$ ENDIF
case default                   
  !z0 внутри круга
  kk=pi*ztc
endselect
cp_calcint_nf1011z=kk
end

function cp_calcint_nf1011(ztc,r,mode,nf)
!аналитическое вычисление интеграла от
!nf=
!   10- dG'_2/dx_0
!   11- dG'_2/dy_0
use pgmod
real(8) r
complex(8) ztc,kk,cp_calcint_nf1011z
integer(4) mode,nf
real(8) cp_calcint_nf1011
kk=cp_calcint_nf1011z(ztc,r,mode)
if (nf==10) then
  cp_calcint_nf1011=dreal(kk)
else
  cp_calcint_nf1011=dimag(kk)
endif
end

subroutine cp_calcint_nf1011_2(ztc,r,mode,val,val2)
!аналитическое вычисление интеграла от
!nf=
!   10- dG'_2/dx_0
!   11- dG'_2/dy_0
use pgmod
real(8) r,val,val2
complex(8) ztc,kk,cp_calcint_nf1011z
integer(4) mode
kk=cp_calcint_nf1011z(ztc,r,mode)
val=dreal(kk)
val2=dimag(kk)
end

function cp_calcint_nf2831(ztc,r,mode,nf)
!аналитическое вычисление интеграла от
!nf=
!   28 - d(G_1*sin(g))/dx_0
!   29 - d(G_1*sin(g))/dy_0
!   30 - d(G_1*cos(g))/dx_0
!   31 - d(G_1*cos(g))/dy_0
use pgmod
real(8) r
complex(8) ztc,kk,dzc
integer(4) mode,nf
real(8) cp_calcint_nf2831
select case (mode)
case (0)     
  !z0 снаружи круга
  if (nf==28.or.nf==30) then
    dzc=-c1
  else
    dzc=ii
  endif
  kk=-pi*r**2/dconjg(ztc)**2*dzc
case (1) 
  !z0 на круге
  !DEC$ IF DEFINED (DEBUG)  !!!сюда не должны приходить
  kk=c0
  if (gs_stop_on_int_dxdy_inbound) call gs_print_stop("Wrong parameters for boundLineType=2 in cp_calcint")
  !DEC$ ENDIF
case default                  
  !z0 внутри круга
  if (nf==28.or.nf==30) then
    dzc=-c1
  else
    dzc=-ii
  endif
  kk=pi*dzc
endselect
kk=dconjg(kk)
if (nf==30.or.nf==31) then
  cp_calcint_nf2831=dreal(kk)
else
  cp_calcint_nf2831=dimag(kk)
endif
end

subroutine cp_calcint_nf2831_2(ztc,r,mode,val,val2,val3,val4)
!аналитическое вычисление интеграла от
!nf=
!   28 - d(G_1*sin(g))/dx_0
!   29 - d(G_1*sin(g))/dy_0
!   30 - d(G_1*cos(g))/dx_0
!   31 - d(G_1*cos(g))/dy_0
use pgmod
real(8) r
complex(8) ztc,kk,kk1,kk2
integer(4) mode
real(8) val,val2,val3,val4
select case (mode)
case (0)     
  !z0 снаружи круга
  kk=-pi*r**2/dconjg(ztc)**2
  kk1=-kk
  kk2=kk*ii
case (1) 
  !z0 на круге
  !DEC$ IF DEFINED (DEBUG)  !!!сюда не должны приходить
  kk1=c0
  kk2=c0
  if (gs_stop_on_int_dxdy_inbound) call gs_print_stop("Wrong parameters for boundLineType=2 in cp_calcint")
  !DEC$ ENDIF
case default                  
  !z0 внутри круга
  kk=pi
  kk1=-kk
  kk2=-kk*ii
endselect
kk1=dconjg(kk1)
kk2=dconjg(kk2)
val=dimag(kk1)
val2=dimag(kk2)
val3=dreal(kk1)
val4=dreal(kk2)
end

function cp_calcint_nf3235(ztc,r,mode,nf)
!аналитическое вычисление интеграла от
!nf=
!   32 - d(G'_1*sin(g))/dx_0
!   33 - d(G'_1*sin(g))/dy_0
!   34 - d(G'_1*cos(g))/dx_0
!   35 - d(G'_1*cos(g))/dy_0
use pgmod
real(8) r
complex(8) ztc,kk,dzc
integer(4) mode,nf
real(8) cp_calcint_nf3235
select case (mode)
case (0)     
  !z0 снаружи круга
  if (nf==32.or.nf==34) then
    dzc=-c1
  else
    dzc=ii
  endif
  kk=pi*r/dconjg(ztc)**2*dzc
case (1) 
  !z0 на круге
  kk=c0
  !DEC$ IF DEFINED (DEBUG)  !!!сюда не должны приходить
  if (gs_stop_on_int_dxdy_inbound) call gs_print_stop("Wrong parameters for boundLineType=2 in cp_calcint")
  !DEC$ ENDIF
case default                  
  !z0 внутри круга
  if (nf==32.or.nf==34) then
    dzc=-c1
  else
    dzc=-ii
  endif
  kk=pi/r*dzc
endselect
kk=dconjg(kk)
if (nf==34.or.nf==35) then
  cp_calcint_nf3235=dreal(kk)
else
  cp_calcint_nf3235=dimag(kk)
endif
end

subroutine cp_calcint_nf3235_2(ztc,r,mode,val,val2,val3,val4)
!аналитическое вычисление интеграла от
!nf=
!   32 - d(G'_1*sin(g))/dx_0
!   33 - d(G'_1*sin(g))/dy_0
!   34 - d(G'_1*cos(g))/dx_0
!   35 - d(G'_1*cos(g))/dy_0
use pgmod
real(8) r
complex(8) ztc,kk,kk1,kk2
integer(4) mode
real(8) val,val2,val3,val4
select case (mode)
case (0)     
  !z0 снаружи круга
  kk=pi*r/dconjg(ztc)**2
  kk1=-kk
  kk2=kk*ii
case (1) 
  !z0 на круге
  kk1=c0
  kk2=c0
  !DEC$ IF DEFINED (DEBUG)  !!!сюда не должны приходить
  if (gs_stop_on_int_dxdy_inbound) call gs_print_stop("Wrong parameters for boundLineType=2 in cp_calcint")
  !DEC$ ENDIF
case default                  
  !z0 внутри круга
  kk=pi/r
  kk1=-kk
  kk2=-kk*ii
endselect
kk1=dconjg(kk1)
kk2=dconjg(kk2)
val=dimag(kk1)
val2=dimag(kk2)
val3=dreal(kk1)
val4=dreal(kk2)
end

function cp_calcint_nf3639(ztc,r,mode,nf)
!аналитическое вычисление интеграла от
!nf=
!   36 - d(G_2*sin(g))/dx_0
!   37 - d(G_2*sin(g))/dy_0
!   38 - d(G_2*cos(g))/dx_0
!   39 - d(G_2*cos(g))/dy_0
use pgmod
real(8) r,ro_c2,r2
complex(8) ztc,ztc1,kk,dzc,dzc1,kk1,kk2
integer(4) mode,nf
real(8) cp_calcint_nf3639
if (nf==36.or.nf==38) then
  dzc=-c1
  dzc1=-c1
else
  dzc=-ii
  dzc1=ii
endif
ztc1=dconjg(ztc)
ro_c2=dreal(ztc*ztc1)
r2=r**2
select case (mode)
case (0)     
  !z0 снаружи круга
  kk1=d2*dlog(ro_c2)*c1
  kk2=(d2*ro_c2-r2)/ztc1**2
  kk=pi*r2/8.0d0*(kk1*dzc+kk2*dzc1)
case (1) 
  !z0 на круге
  kk=c0
  !DEC$ IF DEFINED (DEBUG)  !!!сюда не должны приходить
  if (gs_stop_on_int_dxdy_inbound) call gs_print_stop("Wrong parameters for boundLineType=2 in cp_calcint")
  !DEC$ ENDIF
case default                  
  !z0 внутри круга
  kk1=d2*(ro_c2+r2*(dlog(r2)-d1))*c1
  kk2=ztc**2
  kk=pi/8.0d0*(kk1*dzc+kk2*dzc1)
endselect
kk=dconjg(kk)
if (nf==38.or.nf==39) then
  cp_calcint_nf3639=dreal(kk)
else
  cp_calcint_nf3639=dimag(kk)
endif
end

subroutine cp_calcint_nf3639_2(ztc,r,mode,val,val2,val3,val4)
!аналитическое вычисление интеграла от
!nf=
!   36 - d(G_2*sin(g))/dx_0
!   37 - d(G_2*sin(g))/dy_0
!   38 - d(G_2*cos(g))/dx_0
!   39 - d(G_2*cos(g))/dy_0
use pgmod
real(8) r,ro_c2,r2
complex(8) ztc,ztc1,kk1,kk2,kk_1,kk_2
integer(4) mode
real(8) val,val2,val3,val4,kk
ztc1=dconjg(ztc)
ro_c2=dreal(ztc*ztc1)
r2=r**2
select case (mode)
case (0)     
  !z0 снаружи круга
  kk1=d2*dlog(ro_c2)*c1
  kk2=(d2*ro_c2-r2)/ztc1**2
  kk=pi*r2/8.0d0
  kk_1=-kk*(kk1+kk2)
  kk_2=ii*kk*(-kk1+kk2)
case (1) 
  !z0 на круге
  kk_1=c0
  kk_2=c0
  !DEC$ IF DEFINED (DEBUG)  !!!сюда не должны приходить
  if (gs_stop_on_int_dxdy_inbound) call gs_print_stop("Wrong parameters for boundLineType=2 in cp_calcint")
  !DEC$ ENDIF
case default                  
  !z0 внутри круга
  kk1=d2*(ro_c2+r2*(dlog(r2)-d1))*c1
  kk2=ztc**2
  kk=pi/8.0d0
  kk_1=-kk*(kk1+kk2)
  kk_2=ii*kk*(-kk1+kk2)
endselect
kk_1=dconjg(kk_1)
kk_2=dconjg(kk_2)
val =dimag(kk_1)
val2=dimag(kk_2)
val3=dreal(kk_1)
val4=dreal(kk_2)
end

function cp_calcint_nf4043(ztc,r,mode,nf)
!аналитическое вычисление интеграла от
!nf=
!   40 - d(G'_2*sin(g))/dx_0
!   41 - d(G'_2*sin(g))/dy_0
!   42 - d(G'_2*cos(g))/dx_0
!   43 - d(G'_2*cos(g))/dy_0
use pgmod
real(8) r,ro_c2,r2
complex(8) ztc,ztc1,kk,dzc,dzc1,kk1,kk2
integer(4) mode,nf
real(8) cp_calcint_nf4043
if (nf==40.or.nf==42) then
  dzc=-c1
  dzc1=-c1
else
  dzc=-ii
  dzc1=ii
endif
ztc1=dconjg(ztc)
ro_c2=dreal(ztc*ztc1)
r2=r**2
select case (mode)
case (0)     
  !z0 снаружи круга
  kk1=d2*dlog(ro_c2)*c1
  kk2=(d2*ro_c2-3.0d0*r2)/ztc1**2
  kk=-pi*r/8.0d0*(kk1*dzc+kk2*dzc1)
case (1) 
  !z0 на круге
  kk=c0
  !DEC$ IF DEFINED (DEBUG)  !!!сюда не должны приходить
  if (gs_stop_on_int_dxdy_inbound) call gs_print_stop("Wrong parameters for boundLineType=2 in cp_calcint")
  !DEC$ ENDIF
case default                  
  !z0 внутри круга
  kk1=d2*(ro_c2-r2*(dlog(r2)+d1))*c1
  kk2=ztc**2
  kk=pi/(8.0d0*r)*(kk1*dzc+kk2*dzc1)
endselect
kk=dconjg(kk)
if (nf==42.or.nf==43) then
  cp_calcint_nf4043=dreal(kk)
else
  cp_calcint_nf4043=dimag(kk)
endif
end

subroutine cp_calcint_nf4043_2(ztc,r,mode,val,val2,val3,val4)
!аналитическое вычисление интеграла от
!nf=
!   40 - d(G'_2*sin(g))/dx_0
!   41 - d(G'_2*sin(g))/dy_0
!   42 - d(G'_2*cos(g))/dx_0
!   43 - d(G'_2*cos(g))/dy_0
use pgmod
real(8) r,ro_c2,r2
complex(8) ztc,ztc1,kk1,kk2,kk_1,kk_2
integer(4) mode
real(8) val,val2,val3,val4,kk
ztc1=dconjg(ztc)
ro_c2=dreal(ztc*ztc1)
r2=r**2
select case (mode)
case (0)     
  !z0 снаружи круга
  kk1=d2*dlog(ro_c2)*c1
  kk2=(d2*ro_c2-3.0d0*r2)/ztc1**2
  kk=-pi*r/8.0d0
  kk_1=-kk*(kk1+kk2)
  kk_2=ii*kk*(-kk1+kk2)
case (1) 
  !z0 на круге
  kk_1=c0
  kk_2=c0
  !DEC$ IF DEFINED (DEBUG)  !!!сюда не должны приходить
  if (gs_stop_on_int_dxdy_inbound) call gs_print_stop("Wrong parameters for boundLineType=2 in cp_calcint")
  !DEC$ ENDIF
case default                  
  !z0 внутри круга
  kk1=d2*(ro_c2-r2*(dlog(r2)+d1))*c1
  kk2=ztc**2
  kk=pi/(8.0d0*r)
  kk_1=-kk*(kk1+kk2)
  kk_2=ii*kk*(-kk1+kk2)
endselect
kk_1=dconjg(kk_1)
kk_2=dconjg(kk_2)
val =dimag(kk_1)
val2=dimag(kk_2)
val3=dreal(kk_1)
val4=dreal(kk_2)
end

function cp_calcintm(j,i,nf,ia,knd,knd0)
!интеграл по кругу, когда контрольная точка на обычной панели
use pgmod
!j - номер участка границы (круга), по которой идет интегрирование
!i - номер контрольной точки (середины панели) (связан с номером уравнения)
!nf - номер интеграла
!ia0, knd0 - номера области и границы с контрольной точкой
!ia, knd - номера области и границы, где находится линия интегрирования
integer(4) j,i,nf,knd,knd0,ia
real(8) cp_calcintm,cp_calcint
complex(8) ztc
type(TBound), pointer :: b0
type(TBoundline2), pointer :: bl
b0=>gs%a(ia)%bnd(knd0)
bl=>gs%a(ia)%bnd(knd)%line2(j)
ztc=bl%zc-b0%zc(i)
cp_calcintm=cp_calcint(ztc,bl%r,nf,bl%g0,bl%dir,0)
end

subroutine cp_calcintm_2(j,i,nf,ia,knd,knd0,val,val2)
!интеграл по кругу, когда контрольная точка на обычной панели
use pgmod
!j - номер участка границы (круга), по которой идет интегрирование
!i - номер контрольной точки (середины панели) (связан с номером уравнения)
!nf - номер интеграла
!ia0, knd0 - номера области и границы с контрольной точкой
!ia, knd - номера области и границы, где находится линия интегрирования
integer(4) j,i,nf,knd,knd0,ia
real(8) val,val2
complex(8) ztc
type(TBound), pointer :: b0
type(TBoundline2), pointer :: bl
b0=>gs%a(ia)%bnd(knd0)
bl=>gs%a(ia)%bnd(knd)%line2(j)
ztc=bl%zc-b0%zc(i)
call cp_calcinta_2(ztc,bl%r,nf,bl%g0,bl%dir,0,val,val2)
end

function cp_calcintm1(j,i,i1,k,nf,ia,knd,knd0)
!интеграл по панели, когда контрольная точка на круге
use pgmod
!j - номер панели, по которой идет интегрирование
!i - номер участка на круге с контрольной точкой
!i1 - номер массива контрольных точек на круге 
!k - номер контрольной точки (связан с номером уравнения)
!nf - номер интеграла
!ia0, knd0 - номера области и границы с контрольной точкой
!ia, knd - номера области и границы, где находится панель интегрирования
integer(4) j,i,nf,knd,knd0,ia,i1,k
real(8) cp_calcintm1,calcint
type(TBoundline_Collocate), pointer :: gc
gc=>gs%a(ia)%bnd(knd0)%line2(i)%gc(i1)
cp_calcintm1=calcint(j,dreal(gc%z(k)),dimag(gc%z(k)),nf,ia,knd)
end

subroutine cp_calcintm1_2(j,i,i1,k,nf,ia,knd,knd0,val,val2)
!интеграл по панели, когда контрольная точка на круге
use pgmod
!j - номер панели, по которой идет интегрирование
!i - номер участка на круге с контрольной точкой
!i1 - номер массива контрольных точек на круге 
!k - номер контрольной точки (связан с номером уравнения)
!nf - номер интеграла
!ia0, knd0 - номера области и границы с контрольной точкой
!ia, knd - номера области и границы, где находится панель интегрирования
integer(4) j,i,nf,knd,knd0,ia,i1,k
real(8) val,val2
type(TBoundline_Collocate), pointer :: gc
gc=>gs%a(ia)%bnd(knd0)%line2(i)%gc(i1)
call calcinta_2(j,dreal(gc%z(k)),dimag(gc%z(k)),nf,ia,knd,val,val2)
end

function cp_calcintm2(j,i,i1,k,nf,ia,knd,knd0)
!интеграл по кругу, когда контрольная точка на круге
use pgmod
!j - номер участка границы (круга), по которой идет интегрирование
!i - номер участка на круге с контрольной точкой
!i1 - номер массива контрольных точек на круге 
!k - номер контрольной точки (связан с номером уравнения)
!nf - номер интеграла
!ia0, knd0 - номера области и границы с контрольной точкой
!ia, knd - номера области и границы, где находится линия интегрирования
integer(4) j,i,nf,knd,knd0,ia,i1,k,mode
real(8) cp_calcintm2,cp_calcint
complex(8) ztc
type(TBoundline_Collocate), pointer :: gc
type(TBoundline2), pointer :: bl
gc=>gs%a(ia)%bnd(knd0)%line2(i)%gc(i1)
bl=>gs%a(ia)%bnd(knd)%line2(j)
mode=0
if (knd0==knd.and.i==j) then
  ztc=gc%ztc(k)
  mode=1
else
  ztc=bl%zc-gc%z(k)
endif
cp_calcintm2=cp_calcint(ztc,bl%r,nf,bl%g0,bl%dir,mode)
end

subroutine cp_calcintm2_2(j,i,i1,k,nf,ia,knd,knd0,val,val2)
!интеграл по кругу, когда контрольная точка на круге
use pgmod
!j - номер участка границы (круга), по которой идет интегрирование
!i - номер участка на круге с контрольной точкой
!i1 - номер массива контрольных точек на круге 
!k - номер контрольной точки (связан с номером уравнения)
!nf - номер интеграла
!ia0, knd0 - номера области и границы с контрольной точкой
!ia, knd - номера области и границы, где находится линия интегрирования
integer(4) j,i,nf,knd,knd0,ia,i1,k,mode
real(8) val,val2
complex(8) ztc
type(TBoundline_Collocate), pointer :: gc
type(TBoundline2), pointer :: bl
gc=>gs%a(ia)%bnd(knd0)%line2(i)%gc(i1)
bl=>gs%a(ia)%bnd(knd)%line2(j)
mode=0
if (knd0==knd.and.i==j) then
  ztc=gc%ztc(k)
  mode=1
else
  ztc=bl%zc-gc%z(k)
endif
call cp_calcinta_2(ztc,bl%r,nf,bl%g0,bl%dir,mode,val,val2)
end


function cp_mode(zc,r,z0)
!0 - z0 снаружи круга
!1 - z0 на круге
!2 - z0 внутри круга
use gen_mod
complex(8) zc,z0
real(8) r,rr
integer(4) cp_mode
rr=cdabs(zc-z0)
if (dabs(r-rr)<1.0d-8) then
  cp_mode=1
elseif (rr<r) then 
  cp_mode=2
else
  cp_mode=0
endif
end

subroutine test_cp_cint
use pgmod
integer(4),parameter :: n=1000
real(8) rz(3),cp_calcinta
real(8) ff(n*3,3),cp_calcint_,ff1(n*3,3),errff(n*3,3),rrr(3),maxerr(0:2),err,maxerr0
!real(8) cp_calcint__
integer(4) mode,cp_mode,nf,i,j,k,mode_(n*3,3),j1
complex(8) z,z_(3,n),ztc
data rz /d5,1.0d0,2.0d0/
data rrr /d5,1.0d0,1.5d0/
complex(8) zc
real(8) g0
integer(4) dir
real(8) rtc,time,time1
gs_stop_on_int_dxdy_inbound=.false.
g0=d0
dir=-1
zc=dcmplx(0.0d0,0.0d0)
nf=9
mode_=-1
maxerr=d0
maxerr0=d0
time=rtc()
do k=1,3
  do i=1,3
    do j=1,n
      j1=j+(k-1)*n
      z_(i,j)=rz(i)*cdexp(ii*(j-d1)/n*pi2)
      z=z_(i,j)
	    mode=cp_mode(zc,rrr(k),z)
      if (int_type(nf)>0.and.mode==1) cycle
      mode_(j1,i)=mode
      ztc=zc-z
      !ff(j1,i)=cp_calcint__(zc,rrr(k),z,nf,g0,dir)
      ff(j1,i)=cp_calcinta(ztc,rrr(k),nf,g0,dir,mode)
	    ff1(j1,i)=cp_calcint_(ztc,rrr(k),nf,g0,dir)
      err=dabs(ff(j1,i)-ff1(j1,i))
      if (err>maxerr(mode)) maxerr(mode)=err
      if (err>maxerr0) maxerr0=err
	    errff(j1,i)=err
    enddo
  enddo
enddo
time1=rtc()
time=time1-time
print*,'time=  ',time 
end

subroutine test_cp_cint2
use pgmod
integer(4) nf,mode,i,j,k,n
complex(8) zc,z0,ztc
real(8) f(0:3),f1(0:3),err(0:3),r,cp_calcint_nf2831
n=1000000
nf=28
zc=c1+ii
r=1
do mode=0,2,2
  if (mode==0) then
    z0=zc*d2
  else
    z0=zc*1.5d0
  endif
  ztc=zc-z0
  call init_gs_time
  do j=1,n
    do i=0,3
      f(i)=cp_calcint_nf2831(ztc,r,mode,nf+i)
    enddo
  enddo
  call print_total_time
  call init_gs_time
  do j=1,n
    call cp_calcint_nf2831_2(ztc,r,mode,f1(0),f1(1),f1(2),f1(3))
  enddo
  call print_total_time
  err=f1-f
  write(*,"(4(ES13.3))") (err(k),k=0,3)
enddo
end
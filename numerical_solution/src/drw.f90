subroutine drw_solv3
!вывод в файл решени€ в квадрате
use pgmod2
integer(4) i,j,nnr,nng,ibnd,ix,iy,ia
real(8) xx,yy,f,f_om,fa,foma,pg_f_psioml,s
logical havea,in_bound
integer(4) pg_get_bound_s,q
real(8) pg_get_fun,get_funb,get_fund,get_fund2,pg_get_fun2,get_fun_xy_solve_bi,pg_get_fun_xy,get_funbxy !,pg_get_fun_1
type(TBound_info) bi
real(8) dfxy(2),dfxya(2)
call gs_print('write 2.dat')
call init_gs_time
nng=10
nnr=10
OPEN (1,FILE='2.dat')
write(1,*) 'TITLE = "velosity"'
write(1,*) 'VARIABLES = "X", "Y", "psi", "om", "errpsi", "errom", "dfdx", "dfdy"'
write(1,"('ZONE T=""area"", I=', i0, ', J=', i0, ', F=POINT')") nnr+1,nng+1
do i=1,nng+1
  xx=(i-d1)/nng
  do j=1,nnr+1
    yy=(j-d1)/nnr+square_shift_dy
	  havea=.false.
    dfxy=d0
	  !аналитическое решение
    selectcase (nmain)
    case (1:3,8,15,20,25)
	    fa=get_funb(xx,yy)
	    foma=d0
      dfxya(1)=get_funbxy(xx,yy,1)
      dfxya(2)=get_funbxy(xx,yy,2)
	    havea=.true.
    case (4:7,23,24)
	    fa=get_fund(xx,yy)
	    foma=get_fund2(xx,yy)
	    havea=.true.
	  endselect
	  !численное решение
    in_bound=(i==1).or.(i==nng+1).or.(j==1).or.(j==nnr+1)
    if (in_bound.and.havea) then
      f=fa
      f_om=foma
      dfxy=dfxya
    else
      if (use_subdomain) then
        ix=xx*nsx+1
        if (ix<1) ix=1
        if (ix>nsx) ix=nsx
        iy=yy*nsy+1
        if (iy<1) iy=1
        if (iy>nsy) iy=nsy
        ia=(ix-1)*nsy+iy
        call pg_bind_domain(ia)
        call get_fun_xy_prepare_bi(xx,yy,0,bi)
        f=get_fun_xy_solve_bi(xx,yy,1,d0,d0,bi)
        f_om=get_fun_xy_solve_bi(xx,yy,3,d0,d0,bi)
      else
	      if (in_bound) then
	        q=pg_get_bound_s(xx,yy,ibnd,s)
	        f=pg_f_psioml(ibnd,s,1,0)
	        f_om=pg_f_psioml(ibnd,s,3,0)
          dfxy(1)=pg_get_fun_xy(xx,yy,2,d1,d0,2)
          dfxy(2)=pg_get_fun_xy(xx,yy,2,d0,d1,2)
        else
	        f=pg_get_fun(xx,yy)
	        f_om=pg_get_fun2(xx,yy)
          !dfxy(1)=pg_get_fun_1(xx,yy,1)
          !dfxy(2)=pg_get_fun_1(xx,yy,2)
          call pg_get_fun_xy_gradient(xx,yy,5,1,dfxy)
        endif
      endif
    endif
    WRITE(1,"(F9.5, ' ', F9.5, ' ', F9.5, ' ', F9.5, ' ', E13.5, ' ', E13.5, ' ', F9.5, ' ', F9.5)") xx,yy,f,f_om,dabs(f-fa),dabs(f_om-foma),dfxy(1),dfxy(2)
  enddo
enddo
close(1)
call print_total_time
!2quad.lay
end

subroutine drw_solv1
use pgmod2
real(8) r_circle,par(5)
integer(4) nng
external r_circle
real(8) bndg(nmax)
par(1)=d1
OPEN (1,FILE='2.dat')
write(1,*) 'TITLE = "velosity"'
write(1,*) 'VARIABLES = "X", "Y", "psi", "om", "u", "v", "errpsi", "errom", "erru", "errv"'
nng=nj_   
call drw_solv_cell_e(r_circle,par,nng,.false.,bndg)
close(1)
end

subroutine drw_solv2circ(nng)
use pgmod2
real(8) r_circle,par(5)
external r_circle
integer(4) nng
par(1)=d1
call drw_solv2(r_circle,par,nng,.false.)
end

subroutine drw_solv2mode(mode,par,nng,square_cell)
use pgmod2
integer(4) nng,mode
logical square_cell
real(8) par(5),r_ellipse,r_circle,r_virus,r_poly
external r_ellipse,r_circle,r_virus,r_poly
if (mode==0) then
  call drw_solv2(r_circle,par,nng,square_cell)
elseif (mode==1) then
  call drw_solv2(r_ellipse,par,nng,square_cell)
elseif (mode==2) then
  call drw_solv2(r_virus,par,nng,square_cell)
else
  call drw_solv2(r_poly,par,nng,square_cell)
endif
end

subroutine drw_solv2(fr,par,nng,square_cell)
use pgmod2
real(8) fr,par(5)
external fr
integer(4) nng
real(8) bndg(nmax)
logical square_cell
OPEN (1,FILE='2.dat')
write(1,*) 'TITLE = "velosity"'
write(1,*) 'VARIABLES = "X", "Y", "psi", "om", "u", "v", "errpsi", "errom", "erru", "errv"'
call drw_solv_cell_e(fr,par,nng,square_cell,bndg)
close(1)
end

subroutine drw_solv4circ
use pgmod2
real(8) r_circle,par(5)
external r_circle
par(1)=d1
call drw_solv4(r_circle,par,nj_,.false.)
end

subroutine drw_solv4mode(mode,par,nng,square_cell)
use pgmod2
integer(4) nng,mode
logical square_cell
real(8) par(5),r_ellipse,r_circle,r_virus,r_poly
external r_ellipse,r_circle,r_virus,r_poly
if (mode==0) then
  call drw_solv4(r_circle,par,nng,square_cell)
elseif (mode==1) then
  call drw_solv4(r_ellipse,par,nng,square_cell)
elseif (mode==2) then
  call drw_solv4(r_virus,par,nng,square_cell)
else
  call drw_solv4(r_poly,par,nng,square_cell)
endif
end

subroutine drw_solv4(fr,par,nng,square_cell)
use pgmod2
real(8) fr,par(5)
external fr
logical square_cell
integer(4) nng
real(8) bndg(nmax)
OPEN (1,FILE='4.dat')
write(1,*) 'TITLE = "velosity"'
write(1,*) 'VARIABLES = "X", "Y", "psi", "om", "u", "v", "errpsi", "errom", "erru", "errv"'
call drw_solv_cell_e(fr,par,nng,square_cell,bndg)
call drw_solv_cell_i(fr,par,nng,bndg,2)
close(1)
end

subroutine drw_solv_per_darci2
!вывод в файл решени€ уравнени€ ƒарси с переменной проницаемостью в круговой €чейке
use pgmod2
real(8) r_circle,par(5)
external r_circle
integer(4) i,nng
real(8) bndg(nmax)
OPEN (1,FILE='4.dat')
write(1,*) 'TITLE = "velosity"'
write(1,*) 'VARIABLES = "X", "Y", "psi", "om", "u", "v", "errpsi", "errom", "erru", "errv"'
par(1)=d1
nng=nj_
do i=1,nng+1
  bndg(i)=(i-d1)/nng*pi
enddo
call drw_solv_cell_i(r_circle,par,nng,bndg,1)
close(1)
end

subroutine drw_solv_cell_e(fr,par,nng,square_cell,bndg)
!вывод в файл решени€ во внешней €чейке (кругла€ и квадратна€ €чейка)
use pgmod2
real(8) fr,par(5)
external fr
integer(4) i,j,nng,nnr
logical square_cell,in_bound
real(8) g,xx,yy,rr(nmax)
real(8) bndg(nmax),bndrv(nmax),rrr,rmin
real(8) sred
character(75) str

if (.not.only_analit) call pg_bind_domain(1)
if (drw_e_tr) then
  call drw_solv_cell_e_tr
  return
endif
call get_rr(rr,pi/nng,nnr)
nnr=nnr-1 !число отрезков
if (square_cell) then
  call init_vnesh(nng,bndg,bndrv)
else
  do i=1,nng+1
    bndg(i)=(i-d1)/nng*pi
	  bndrv(i)=h
  enddo
endif
write(1,"('ZONE T=""area_e"", I=', i0, ', J=', i0, ', F=POINT')") nnr+1,nng+1
call init_gs_time
do i=1,nng+1
  write(str,"('Write file area_e. ', I0, ' (from ', I0, ')')") i, nng+1
  call write_by_time(str)
  g=bndg(i)
  rmin=fr(g,par)
  do j=1,nnr+1
    rrr=sred(h2,h,rmin,bndrv(i),rr(j))
    xx=rrr*dcos(g)
 	  yy=rrr*dsin(g)
	  in_bound=(i==1).or.(i==nng+1).or.(j==1).or.(j==nnr+1)
    call drw_solv_cell_e_point(xx,yy,rrr,g,in_bound,gs_test_point_near_bound)
  enddo
enddo
call print_total_time
end

subroutine drw_solv_cell_e_tr
!вывод в файл решени€ во внешней €чейке на треугольной сетке
use pgmod2
integer(4) i,j,jprev
real(8) xx,yy,rrr,g
logical in_bound
character(75) str
if (gsarea%a%npe==3) then
  write(1,"('ZONE T=""area_e"" N=', i0, ', E=', i0, ', F=FEPOINT, ET=TRIANGLE')") gsarea%a%n, gsarea%a%ntr
else
  write(1,"('ZONE T=""area_e"" N=', i0, ', E=', i0, ', F=FEPOINT, ET=QUADRILATERAL')") gsarea%a%n, gsarea%a%ntr
endif
jprev=-1
call init_gs_time
do i=1,gsarea%a%n
  j=i*100/gsarea%a%n
  if (j>jprev) then
    write(str,"('Write file area_e. ', I0, '%')") j
	call write_by_time(str)
	jprev=j
  endif
  xx=dreal(gsarea%a%zm(i))
  yy=dimag(gsarea%a%zm(i))
  rrr=dsqrt(xx**2+yy**2)
  g=datan2(yy,xx)
  if (g<-1.0d-3) g=g+pi2
  call drw_solv_cell_e_point(xx,yy,rrr,g,in_bound,.true.)
enddo
do i=1,gsarea%a%ntr
  WRITE(1,"(4(' ',i0))") (gsarea%a%trm(j,i),j=1,gsarea%a%npe)
enddo
call print_total_time
end

subroutine drw_solv_cell_e_point(xx,yy,rrr,g,in_bound,find_inbound)
!вывод в файл решени€ во внешней €чейке (кругла€ и квадратна€ €чейка)
!текуща€ точка
use pgmod2
logical in_bound,find_inbound
real(8) g,xx,yy,xx1,yy1,f,fom
real(8) get_funa,get_func,get_funa_oss,get_fun_xy_solve_bi
real(8) psi1,om,u,v !,fcell
real(8) rrr
real(8) fu,fv,psi_inf,ff(2)
type(TBound_info) bi
integer(4) mode
	if (nmain==16) then
	  xx1=dlog(rrr)
	  yy1=g
	else
	  xx1=xx
	  yy1=yy
	endif
	!численное решение
  
  if (.not.only_analit) then
    mode=0
    if (.not.find_inbound) then
      mode=1
      if (in_bound) mode=2
    endif
    call get_fun_xy_prepare_bi(xx,yy,mode,bi)
    f=get_fun_xy_solve_bi(xx,yy,1,d0,d0,bi)
    if (.not.solve_only_psi) then
      fom=-get_fun_xy_solve_bi(xx,yy,3,d0,d0,bi)
      call get_fun_xy_solve_gradient_bi(xx,yy,5,d0,d0,bi,ff)
      fu=ff(2)
      fv=-ff(1)
      !call get_fun_xy_solve_gradient_bi(xx,yy,6,d0,d0,bi,ff)
      !fcell=-ff(1)*fu-ff(2)*fv
    endif
  	if (nmain==21) then
  	  if (yy<1.0d-8) then
  	    fu=d1
  		  fv=d0
  	  else
        fu=fu/yy1
        fv=fv/yy1
  	  endif
  	elseif (nmain==22) then !внешн€€ задача
  	  f=f+psi_inf(xx1,yy1,1,0,1)
  	  fu=fu+psi_inf(xx1,yy1,1,1,1)
  	  fv=fv+psi_inf(xx1,yy1,1,2,1)
    endif
  endif
  !аналитическое решение
	psi1=f
  om=fom
  u=fu
	v=fv
	!if (.not.in_bound) then
 	  if (nmain==11) then
	    psi1=get_funa(rrr,g)
	  elseif (nmain==12.and.nfun==33) then
	    call Get_kuv_uvom_psi(xx,yy,u,v,om,psi1)
	  elseif (nmain==13) then
	    psi1=get_func(dcmplx(xx,yy))
	  elseif (nmain==14.and.nfun==34) then
	    call Get_kuv_uvom_psi(xx,yy,u,v,om,psi1)
    elseif (nmain==16.or.nmain==17) then
	    call Get_kuv_uvom_psi(xx,yy,u,v,om,psi1)
	  elseif (nmain==18.and.nfun==35) then
	    call Get_kuv_uvom_psi(xx,yy,u,v,om,psi1)
	  elseif (nmain==19.and.nfun==36) then
  		call Get_kuv_uvom_psi_bri(xx,yy,u,v,om,psi1)
	  elseif (nmain==21) then
	    psi1=get_funa_oss(rrr,g,u,v)
    elseif (nmain==22.and.nfun==40) then
	    call get_fune(xx1,yy1,psi1,u,v)
	  endif
  !endif
  if (only_analit) f=psi1
	WRITE(1,"(F13.9, ' ', F13.9, ' ', E15.8, ' ', E15.8, ' ', E15.8, ' ', E15.8, ' ', F13.9, ' ', F13.9, ' ', F13.9, ' ', F13.9)") xx1,yy1,f,fom,fu,fv,dabs(f-psi1),dabs(fom-om),dabs(fu-u),dabs(fv-v)
  !WRITE(1,"(F13.9, ' ', F13.9, ' ', E15.8, ' ', E15.8, ' ', E15.8, ' ', E15.8, ' ', F13.9, ' ', F13.9, ' ', F13.9, ' ', F13.9)") xx1,yy1,f,fom,fu,fv,fcell,dabs(fom-om),dabs(fu-u),dabs(fv-v)
end

subroutine drw_solv_cell_i(fr,par,nng,bndg,ia)
!вывод в файл решени€ во внутренней €чейке
use pgmod2
real(8) fr,par(5)
external fr
integer(4) i,j,nng,nnr,ibnd,ia
logical in_bound
integer(4) in_bound1,pg_get_bound_s_tt
real(8) g,xx,yy,f,pg_get_fun,fom,pg_get_fun2,rr(nmax),pg_f_psioml,s !,pg_get_fun_1
real(8) get_func,sred,tt
real(8) psi1,om,u,v
real(8) bndg(nmax),rrr,rmin
real(8) fu,fv,ff(2)
complex(8) uvz,ttl
character(75) str
if (.not.only_analit) then
  call pg_bind_domain(ia)
  call pg_bind_bound(1)
endif
nnr=nng/pi
do i=1,nnr+1
  rr(i)=(i-d1)/nnr
enddo
write(1,"('ZONE T=""area_i"", I=', i0, ', J=', i0, ', F=POINT')") nnr+1,nng+1
call init_gs_time
do i=1,nng+1
  write(str,"('Write file area_i. ', I0, ' (from ', I0, ')')") i, nng+1
  call write_by_time(str)  
  g=bndg(i)
  rmin=fr(g,par)
  do j=1,nnr+1
    rrr=sred(d0,d1,d0,rmin,rr(j))
    xx=rrr*dcos(g)
 	  yy=rrr*dsin(g)
    if (.not.only_analit) then
	    !численное решение
	    in_bound=(i==1).or.(i==nng+1).or.(j==1).or.(j==nnr+1)
	    if (in_bound) then
        in_bound1=pg_get_bound_s_tt(xx,yy,ibnd,s,tt)
	      ttl=cdexp(ii*tt)
	      f=pg_f_psioml(ibnd,s,1,0)
        if (.not.solve_only_psi) then
	        fom=-pg_f_psioml(ibnd,s,3,0)
	        uvz=-dcmplx(pg_f_psioml(ibnd,s,2,0),pg_f_psioml(ibnd,s,1,1))*ttl
	        fu=dreal(uvz)
          fv=dimag(uvz)
        endif
	    else
	      f=pg_get_fun(xx,yy)
        if (.not.solve_only_psi) then
	        fom=-pg_get_fun2(xx,yy)
	        !fu=pg_get_fun_1(xx,yy,2)
	        !fv=-pg_get_fun_1(xx,yy,1)
          call pg_get_fun_xy_gradient(xx,yy,5,1,ff)
          fu=ff(2)
          fv=-ff(1)
        endif
      endif
    endif
    !аналитическое решение
	  psi1=f
    om=fom
    u=fu
	  v=fv
	  !if (.not.in_bound) then
	    if (nmain==13) then
	      psi1=get_func(dcmplx(xx,yy))
	    elseif (nmain==14.and.nfun==34) then
	      call Get_kuv_uvpsi_c(xx,yy,u,v,psi1)
	      om=d0
	    elseif (nmain==18.and.nfun==35) then
	      call Get_kuv_uvom_psi_c_bri(xx,yy,u,v,om,psi1)
      elseif (nmain==19.and.nfun==36) then
	      call Get_kuv_uvom_psi_c_bri(xx,yy,u,v,om,psi1)
      elseif (nmain==10.and.nfun==32) then
	  	  call Get_kuv_uvom_psi_c_per_darci(xx,yy,u,v,om,psi1)
      endif
    !endif
    if (only_analit) f=psi1
	  WRITE(1,"(F13.9, ' ', F13.9, ' ', E15.8, ' ', E15.8, ' ', E15.8, ' ', E15.8, ' ', F13.9, ' ', F13.9, ' ', F13.9, ' ', F13.9)") xx,yy,f,fom,fu,fv,dabs(f-psi1),dabs(fom-om),dabs(fu-u),dabs(fv-v)
  enddo
enddo
call print_total_time
end

subroutine save_vx
use pgmod2
integer(4) i,n,mode
real(8) y,pg_get_fun_xy
OPEN (1,FILE='vx.dat')
n=100
do i=0,n
  mode=1
  if (i==0.or.i==n) mode=2
  y=(i+d0)/n*(h-d1)+d1
  WRITE(1,"(F9.5, ' ', F9.5)") y,pg_get_fun_xy(d0,y,2,d0,d1,mode)
enddo
close(1)
end

subroutine drw_p
use pgmod2
real(8), allocatable :: p(:), s(:)
real(8) ds,eps,k,pg_int_psioml,yy,dp,ga_gradp_oneline_ds,ds2
integer(4) n,n1,i,j
k=s_darci**-2
eps=1.0d-6
ds=0.01d0
n=(h+eps)/ds
ds=h/n
allocate(p(n+1),s(n+1))
p(1)=d0
s(1)=d0
n1=-1
do i=1,n+1
  s(i)=(i-1)*ds
  if (s(i)>d1-eps.and.n1<0) then
    s(i)=d1
    n1=i
  endif
enddo
call pg_bind_domain(2)
call pg_bind_bound(1)
ds2=gsbnd%s(gsbnd%npanel+1)-d1
do i=1,n1-1
  p(i+1)=p(i)+pg_int_psioml(1,s(i)+ds2,s(i+1)+ds2,2)/k
enddo
call pg_bind_domain(1)
do i=n1,n
  p(i+1)=p(i)-pg_int_psioml(1,s(i)-d1,s(i+1)-d1,4)
enddo
j=(n1+n+1)/2
yy=1.5d0
ds=1.0d-3
dp=ga_gradp_oneline_ds(d0,yy,s(j),yy,ds)
dp=dp+ga_gradp_oneline_ds(s(j),yy,s(j),d0,ds)
dp=-dp
dp=dp-p(j)
do i=n1+1,n+1
  p(i)=p(i)+dp
enddo
OPEN (1,FILE='px.dat')
write(1,"('ZONE T=""p=p"", I=', i4, ', F=POINT')") 2*n+1
do i=n+1,2,-1
  write(1,"(2(E15.6))") -s(i),-p(i)
enddo
do i=1,n+1
  write(1,"(2(E15.6))") s(i),p(i)
enddo
deallocate(p,s)
close(1)
end
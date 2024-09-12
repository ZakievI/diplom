program main
use mod
nf=2
emptyBoundEq=.true. !уравнения для граничных ребер не составляются
use_subdomain=.false. !использовать подобласти
call main1
call deallocate_mem
end

subroutine deallocate_mem
use mod
if (allocated(ff)) deallocate(ff)
if (allocated(ffknow)) deallocate(ffknow)
if (allocated(err)) deallocate(err)
if (allocated(cell_domain)) deallocate(cell_domain)
end

subroutine main1
use mod
real(8) dfdl,func,ds2,dfdl2,bndg(200),bndrv(200),h
integer(4) i,n,ng
external dfdl,dfdl2
gsareapart=>ap

ng=60
h=5*d1

ds=0.1d0
ds2=0.01d0
call ga_init_vneshg(ng,bndg,bndrv,h,h,8)
call ga_init_mesh_rcell_quads(bndrv,bndg,ng+1,ng+1)
call ga_mesh_square(d1,h,ng,bndrv,bndg)
call ga_drw_trmesh(1)

!call ga_init_mesh_quadcell_ds(d0,d1,d0,d1,ds,4)

n=nint(d1/ds)
allocate(ff(ap%n),ffknow(ap%n),err(ap%n))
ff=d0
ffknow=.false.
do i=1,1+n
  ffknow(i)=.true.
  ff(i)=func(dreal(ap%zm(i)),dimag(ap%zm(i)))
enddo
call intf_dxdy(ap%zm,ap%trm,ap%n,ap%ntr,ap%npe,ff,ffknow,dfdl,ds2)
call getErr
end

subroutine getErr
use mod
integer(4) i
real(8) func
maxErr=d0
do i=1,ap%n
  err(i)=dabs(ff(i)-func(dreal(ap%zm(i)),dimag(ap%zm(i))))
  if (err(i)>maxErr) maxErr=err(i)
enddo
print *,maxErr
end

function dfdl(x,y,xder,yder)
use mod
real(8) dfdl,x,y,xder,yder,der1,der2,eps,f,dfdxy
!dfdl=dfdl2(x,y,xder,yder,0,0,0,0,.false.)
complex(8) lder
eps=1.0d-8
lder=dcmplx(xder,yder)
lder=lder/cdabs(lder)
der1=dreal(lder)
der2=dimag(lder)
f=d0
if (dabs(der1)>eps) f=dfdxy(x,y,1)*der1
if (dabs(der2)>eps) f=f+dfdxy(x,y,2)*der2
dfdl=f
end

function dfdl2(x,y,xder,yder,itr,ib,k1,k2,is_bound)
use mod
real(8) dfdl2,x,y,xder,yder,der1,der2,eps,f,dfdxy
complex(8) lder
integer(4) itr,ib,k1,k2
logical is_bound
if (emptyBoundEq.and.is_bound.and.y<1.0d-5) then
  dfdl2=ifdxdy_emptyval
  return
endif
eps=1.0d-8
lder=dcmplx(xder,yder)
lder=lder/cdabs(lder)
der1=dreal(lder)
der2=dimag(lder)
f=d0
if (dabs(der1)>eps) f=dfdxy(x,y,1)*der1
if (dabs(der2)>eps) f=f+dfdxy(x,y,2)*der2
dfdl2=f
end

function func(x,y)
use mod
real(8) func,x,y
if (nf==1) then
  func=x*y
elseif (nf==2) then
  func=dsin(x)*y
endif
end

function dfdxy(x,y,k)
use mod
real(8) dfdxy,x,y
integer(4) k
if (nf==1) then
  if (k==1) then
    dfdxy=y
  else
    dfdxy=x
  endif
elseif (nf==2) then
  if (k==1) then
    dfdxy=y*dcos(x)
  else
    dfdxy=dsin(x)
  endif
endif
end
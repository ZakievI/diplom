subroutine ani_getnd(am,nd,labelD,nmax)
!получить количество областей из треугольной сетки
use ani_mod
type(Tani_mesh) am
integer(4) nd  !количество областей
integer(4) nmax
integer(4) labelD(nmax)  !индексы областей
integer(4) imax,i
logical, allocatable :: id(:)
imax=maxval(am%labelT(1:am%nt))
allocate(id(imax))
id=.false.
do i=1,am%nt
  id(am%labelT(i))=.true.
enddo
nd=0
do i=1,imax
  if (id(i)) then
    nd=nd+1
    labelD(nd)=i
  endif
enddo
deallocate(id)
end

subroutine pgani_addareapart(am, dind)
!создать текущий участок границы из треугольной сетки am с индексом подобласти dind
use pgmod
use ani_mod
type(Tani_mesh) am    !треугольная сетка
integer(4) dind       !текущий индекс подобласти, которую надо добавить в виде участка области
integer(4) ntr,i
integer(4), allocatable :: tr(:,:)
ntr=0
do i=1,am%nt
  if (am%labelT(i)==dind) ntr=ntr+1
enddo
allocate(tr(3,ntr))
do i=1,am%nt
  if (am%labelT(i)==dind) then
    tr(1,i)=am%tri(3,i)
    tr(2,i)=am%tri(2,i)
    tr(3,i)=am%tri(1,i)
  endif
enddo
call pg_allocate_areapart_geom(am%nv,ntr,3)
call pg_init_areapart_geom_xy(am%vrt(1,:),am%vrt(2,:))
call pg_init_areapart_geom_tr(tr)
deallocate(tr)
end

subroutine pgani_createarea(ntmax)
!построить сетку для текущей gsarea
!по границе
use pgmod
use ani_mod
integer(4) ntmax,nvmax,nbmax
integer(4) i,j,j1,nvr
type(Tbound), pointer :: b
type(Tani_bound) ab
type(Tani_mesh) am
nvr=0
do i=1,gsarea%nb
  b=>gsarea%bnd(i)
  nvr=nvr+b%npanel
enddo
nvr=nvr+gsarea%nb
call allocate_ani_bound(ab,0,nvr)
j=1
do i=1,gsarea%nb
  b=>gsarea%bnd(i)
  !j - индекс первой точки границы
  j1=j+b%npanel !индекс последней точки границы
  ab%vbr(1,j1:j:-1)=b%x(:)
  ab%vbr(2,j1:j:-1)=b%y(:)
  ab%vbr(1,j)=b%x(1)
  ab%vbr(2,j)=b%y(1)
  j=j1+1
enddo
call nbnv_fromntmax(nvmax,ntmax,nbmax)
call allocate_ani_mesh(am,nvmax,ntmax,nbmax,.false.)
call ani_aft2dfront(ab, am)
call pg_allocate_area(1)
call pg_bind_areapart(1)
call pgani_addareapart(am, 1)
call deallocate_ani_bound(ab)
call deallocate_ani_mesh(am)
end

subroutine pgani_createarea_a(ntmax,h)
!построить сетку для текущей gsarea
!по границе через аналитическое задание
use pgmod
use ani_mod
real(8) h !размер ячейки сетки
integer(4) ntmax,nvmax,nbmax
integer(4) i,j,j1,nbv,nbl,k
type(Tbound), pointer :: b
type(Tani_bound_a) aba
type(Tani_mesh) am
nbv=0
do i=1,gsarea%nb
  b=>gsarea%bnd(i)
  nbv=nbv+b%npanel
enddo
nbl=nbv
call allocate_ani_bound_a(aba,Nbv,Nbl)
j=1
do i=1,gsarea%nb
  b=>gsarea%bnd(i)
  !j - индекс первой точки границы
  j1=j+b%npanel-1 !индекс предпоследней точки границы
  aba%bv(1,j1:j:-1)=b%x(1:b%npanel)
  aba%bv(2,j1:j:-1)=b%y(1:b%npanel)
  do k=j,j1
    aba%bl(1,k)=k
    aba%bl(2,k)=k+1
    aba%bl(3,k)=0
    aba%bl(4,k)=-1
    aba%bl(5,k)=i
    aba%bl(6,k)=1
    aba%bl(7,k)=0
  enddo
  aba%bl(2,j1)=j !последняя вершина совпадает с первой
  j=j1+1
enddo
call nbnv_fromntmax(nvmax,ntmax,nbmax)
call allocate_ani_mesh(am,nvmax,ntmax,nbmax,.true.)
call ani_aft2dboundary(aba, h, am)
call pg_allocate_area(1)
call pg_bind_areapart(1)
call pgani_addareapart(am, 1)
call deallocate_ani_bound_a(aba)
call deallocate_ani_mesh(am)
end

subroutine nbnv_fromntmax(nvmax,ntmax,nbmax)
use gen_mod
integer(4) nvmax,ntmax,nbmax
nvmax=ntmax/2
nbmax=int(4*dsqrt(ntmax*d1))
end

subroutine pgani_meshfromarea(am,ntmax)
!скопировать текущую сетку gsarea в am
use pgmod
use ani_mod
type(Tani_mesh) am
type(areatype), pointer :: a
integer(4) nvmax,ntmax,nbmax
a=>gsarea%a
if (a%npe>3) call gs_print_stop("Error pgani_meshfromarea. Only triangle mesh")
call nbnv_fromntmax(nvmax,ntmax,nbmax)
call allocate_ani_mesh(am,nvmax,ntmax,nbmax,.true.)
am%nv=a%n
am%vrt(1,1:am%nv)=dreal(a%zm(1:am%nv))
am%vrt(2,1:am%nv)=dimag(a%zm(1:am%nv))
am%nt=a%ntr
am%tri(:,1:a%ntr)=a%trm(:,1:a%ntr)
am%labelT(1:a%ntr)=1
end
!сохранение расчета

subroutine pg_save_solve(fname,save_comment)
!dec$ attributes dllexport:: pg_save_solve
!сохранение расчета
use pgmod
integer(4), parameter :: vers=3  !версия выгрузки
                                 !1 - первая версия
								 !2 - добавлены boundLineType=2 
                 !3 - constvals перенесены в TGreenSolve
character(*) fname, save_comment
integer(4) i0,i,j,j1,k,k1
type(TGreenSolve), pointer :: gs1  !текущая задача
type(TArea), pointer :: gsarea1    !текущая область текущей задачи
type(TBound), pointer :: gsbnd1    !текущая граница    
type(TBoundline), pointer :: gsbndl1    !текущий участок границы
type(TBoundline2), pointer :: gsbndl_2    !текущий участок границы2
type(TBounLineFuncApprox), pointer :: ga
integer(4) val(8)
call DATE_AND_TIME(VALUES=val)
OPEN (1,FILE=fname)
write(1,"('-- pg save file ', I2.2, '.', I2.2, '.', I4.4, ' ',I2.2, ':', I2.2, ':', I2.2)") val(3),val(2),val(1),val(5),val(6),val(7)
write(1,"(i0,' -- version')") vers
write(1,*)"--: ",save_comment
write(1,"('----------------------------')")
write(1,"(i0,' -- problems count')") gsMain%ngs
do i0=1,gsMain%ngs
  write(1,"('-- PROBLEM ', i0)") i0
  gs1=>gsMain%ggs(i0)
  write(1,"(i0,' -- domains count')") gs1%na
  do i=1,gs1%na
    write(1,"('-- domain ', i0)") i
    gsarea1=>gs1%a(i)
	  write(1,"(i0,' -- domain equation')") gsarea1%type_eq(1)
	  if (gsarea1%type_eq(1)==8.or.gsarea1%type_eq(1)==11) write(1,"(i0,' -- domain equation2')") gsarea1%type_eq(2)
	  write(1,"('-- GEOMETRY')")
    write(1,"(i0,' -- bounds count')") gsarea1%nb
	  do j=1,gsarea1%nb
	    gsbnd1=>gsarea1%bnd(j)
	    write(1,"('-- bound ', i0)") j
	    write(1,"(i0, ' -- boundLineType')") gsbnd1%boundLineType
	    write(1,"(i0, ' -- boundlines count')") gsbnd1%nline
	    if (gsbnd1%boundLineType==1) then
	      do k=1,gsbnd1%nline
	        gsbndl1=>gsbnd1%line(k)
	        write(1,"('-- boundline ', i0)") k
	        write(1,"(i0, ' -- boundline format')") 1 !формат - по координатам
	        write(1,"(i0, ' -- boundline npanel')") gsbndl1%npanel
	        do j1=1,gsbndl1%npanel+1
	          write(1,*) gsbndl1%x(j1),' ',gsbndl1%y(j1)
	        enddo
		    enddo
      elseif (gsbnd1%boundLineType==2) then
	      do k=1,gsbnd1%nline
	        gsbndl_2=>gsbnd1%line2(k)
	        write(1,"('-- boundline ', i0)") k
	        write(1,*) dreal(gsbndl_2%zc),' ',dimag(gsbndl_2%zc),' zc'
		      write(1,*) gsbndl_2%r,' r'
		      write(1,*) gsbndl_2%g0,' g0'
		      write(1,*) gsbndl_2%dir,' dir'
		    enddo
	    endif
	  enddo
	  write(1,"('-- SOLVE RESULT')")
	  do j=1,gsarea1%nb
	    gsbnd1=>gsarea1%bnd(j)
	    write(1,"('-- bound ', i0)") j
	    select case (gsbnd1%boundLineType)
	    case(1)
        do k=1,gsbnd1%npanel
	        write(1,"(20E24.15)") (gsbnd1%psiom(k,j1),j1=1,gsarea1%umax)
	      enddo
	    case(2)
	      do k=1,gsbnd1%nline
		      write(1,"('-- boundline ', i0)") k
	        gsbndl_2=>gsbnd1%line2(k)
	        do k1=1,gsarea1%umax
		        ga=>gsbndl_2%ga(k1)
			      write(1,*) ga%typea, ga%bnda, ga%n
		        write(1,"(20E24.15)") (ga%valkk(j1),j1=1,ga%n)
		      enddo
		    enddo
	    end select
	  enddo
  enddo
  write(1,"('-- SOLVE RESULT constvals')")
  write(1,"(i0, ' -- constvaln')") gs%constvaln
  do k=1,gs%constvaln
    write(1,*) gs%constvals(k)
  enddo
enddo
!write(1,"('-- SOLVE RESULT constvals')")
!write(1,"(i0, ' -- constvaln')") gs%constvaln
!do k=1,gsMain%constvaln
!  write(1,*) gsMain%constvals(k)
!enddo
close(1)
end

subroutine pg_read_solve(fname)
!dec$ attributes dllexport:: pg_read_solve
!чтение расчета
use pgmod
integer(4) vers !версия выгрузки
integer(4) ival,ival2,i0,i,j,j1,k,blcount,k1
real(8), allocatable :: x(:),y(:)
real(8) x0,y0,r,g0
character(*) fname
type(TBounLineFuncApprox), pointer :: ga
call pg_deallocate_mem
OPEN (1,FILE=fname)
read(1,*) !--pg save file ...
read(1,*) vers !-- version
read(1,*) !--: save_comment
read(1,*) !----------------------------
read(1,*) ival !-- problems count
call pg_allocate_problems(ival)
do i0=1,gsMain%ngs
  read(1,*) !'-- PROBLEM i0
  call pg_bind_problem(i0)
  read(1,*) ival !-- domains count
  call pg_allocate_domains(ival)
  do i=1,gs%na
    read(1,*) !-- domain i
    call pg_bind_domain(i)
	  read(1,*) ival !-- domain equation
	  if (ival==8.or.ival==11) then
	    read(1,*) ival2 !-- domain equation2
	    call set_domain_equation(ival,ival2)
	  else
	    call pg_set_domain_equation(ival)
	  endif
	  read(1,*) !-- GEOMETRY
	  read(1,*) ival !-- bounds count
	  call pg_allocate_bounds(ival)
	  do j=1,gsarea%nb
	    read(1,*) !--bound j
	    call pg_bind_bound(j)
	    if (vers>=2) then
	      read(1,*) gsbnd%boundLineType
	    else
	      gsbnd%boundLineType=1
	    endif
	    read(1,*) blcount !-- boundlines count
	    select case (gsbnd%boundLineType)
	    case(1)
	      call pg_allocate_boundlines(blcount)
	      do k=1,gsbnd%nline
	        read(1,*) !--boundline k
		    read(1,*) ival !-- boundline format
	        select case (ival)
	        case(1) !формат - по координатам
	          read(1,*) ival !-- boundline npanel
		      allocate(x(ival+1),y(ival+1))
	          do j1=1,ival+1
	            read(1,*) x(j1),y(j1)
	          enddo
		      call pg_init_boundline_geom(k,ival,x,y)
		      deallocate(x,y)
	        end select
	      enddo
	    case(2)
	      call pg_allocate_boundlines2(blcount)
	      do k=1,gsbnd%nline
	        read(1,*) !--boundline k
		      read(1,*) x0,y0
		      read(1,*) r
		      read(1,*) g0
		      read(1,*) ival
		      call pg_init_cylider_particle(k,x0,y0,r,g0,ival,0)
		    enddo
	    end select
	  enddo
	  call pg_geom_postprocessor
	  read(1,*) !-- SOLVE RESULT
	  do j=1,gsarea%nb
	    call pg_bind_bound(j)
	    read(1,*) !-- bound j
	    select case (gsbnd%boundLineType)
	    case(1)
	      do k=1,gsbnd%npanel
	        read(1,*) (gsbnd%psiom(k,j1),j1=1,gsarea%umax)
	      enddo
	    case(2)
	      call pg_allocate_bound_gu2
	      read(1,*) !--boundline k
	  	  do k1=1,gsbnd%nline
	  	    call pg_bind_boundline(k1)
	  	    do k=1,gsarea%umax
	  	      ga=>gsbndl2%ga(k)
	  	      read(1,*) ga%typea, ga%bnda, ga%n
	  	      call pg_allocate_and_set_ga2(k,ga%n)
	  	      read(1,*) (ga%valkk(j1),j1=1,ga%n)
	  	    enddo
	  	  enddo
	    end select
	  enddo
  enddo
  if (vers>=3) then
    read(1,*) !-- SOLVE RESULT constvals
    read(1,*) ival !-- constvaln
    call pg_allocate_constvalind(ival)
    do k=1,gs%constvaln
      read(1,*) gs%constvals(k)
    enddo
  endif
enddo
close(1)
call pg_get_psioml
end

subroutine dns_drw_one_cell(nfile,mcell,ic,jc,nng,nnr,bndg,bndrv,rr,nchar,draw_uv)
!вывод в файл решения во внешней ячейке (круглая и квадратная ячейка)
use dnsmod_base
type(mcelltype) mcell
integer(4) nng,nnr,ic,jc,nfile
real(8) rr(nnr+1)
real(8) bndg(nng+1),bndrv(nng+1),bndrv_partcl(nng+1)  
character(1) nchar
logical draw_uv
bndrv_partcl=mcell%r
call dns_drw_one_cell_partcl(nfile,mcell,ic,jc,nng,nnr,bndg,bndrv_partcl,bndrv,rr,nchar,draw_uv)
end

subroutine dns_drw_one_cell_partcl(nfile,mcell,ic,jc,nng,nnr,bndg,bndrv_partcl,bndrv,rr,nchar,draw_uv)
!вывод в файл решения во внешней ячейке (круглая и квадратная ячейка)
!для произвольной формы частицы
use dnsmod_base
type(mcelltype) mcell
integer(4) nng,nnr,ic,jc,nfile
real(8) rr(nnr+1)
real(8) bndg(nng+1),bndrv(nng+1),bndrv_partcl(nng+1)  
character(1) nchar
logical draw_uv
external dns_write_point_psi_om_uv,dns_write_point_psi_om
if (draw_uv) then
  call dns_drw_one_cell2_partcl(nfile,mcell,ic,jc,nng,nnr,bndg,bndrv_partcl,bndrv,rr,nchar,dns_write_point_psi_om_uv)
else
  call dns_drw_one_cell2_partcl(nfile,mcell,ic,jc,nng,nnr,bndg,bndrv_partcl,bndrv,rr,nchar,dns_write_point_psi_om)
endif
end

subroutine dns_drw_one_cell2(nfile,mcell,ic,jc,nng,nnr,bndg,bndrv,rr,nchar,write_point)
!вывод в файл решения во внешней ячейке (круглая и квадратная ячейка)
use dnsmod_base
integer(4) nfile      !номер файла
type(mcelltype) mcell !информация об ЭЯ
integer(4) ic,jc      !индекс ячейки по x и y
integer(4) nng        !количество отрезков по углу
integer(4) nnr        !количество отрезков по радиусу
real(8) bndg(nng+1)   !разбиение по углу (вычисляется в dns_prepare_drw_one_cell_g)
real(8) bndrv(nng+1)  !расстояние от внешней границы до центра по улгу (вычисляется в dns_prepare_drw_one_cell_g)
real(8) bndrv_partcl(nng+1)  !расстояние от внутренней границы до центра по улгу
real(8) rr(nnr+1)     !разбиение по радиусу (вычисляется в dns_prepare_drw_one_cell_r)
character(1) nchar    !буква (индекс) в названии зоны
!write_point(nfile,xx,yy,r,g,mode) - процедура для вывода значения в точке
external write_point
bndrv_partcl=mcell%r
call dns_drw_one_cell2_partcl(nfile,mcell,ic,jc,nng,nnr,bndg,bndrv_partcl,bndrv,rr,nchar,write_point)
end

subroutine dns_drw_one_cell2_partcl(nfile,mcell,ic,jc,nng,nnr,bndg,bndrv_partcl,bndrv,rr,nchar,write_point)
!вывод в файл решения во внешней ячейке (круглая и квадратная ячейка)
!для произвольной формы частицы
use dnsmod_base
use pgmod
integer(4) nfile      !номер файла
type(mcelltype) mcell !информация об ЭЯ
integer(4) ic,jc      !индекс ячейки по x и y
integer(4) nng        !количество отрезков по углу
integer(4) nnr        !количество отрезков по радиусу
real(8) bndg(nng+1)   !разбиение по углу (вычисляется в dns_prepare_drw_one_cell_g)
real(8) bndrv(nng+1)  !расстояние от внешней границы до центра по улгу (вычисляется в dns_prepare_drw_one_cell_g)
real(8) bndrv_partcl(nng+1)  !расстояние от внутренней границы до центра по улгу
real(8) rr(nnr+1)     !разбиение по радиусу (вычисляется в dns_prepare_drw_one_cell_r)
character(1) nchar    !буква (индекс) в названии зоны
!write_point(nfile,xx,yy,r,g,mode) - процедура для вывода значения в точке
integer(4) mode ! см get_fun_xy_prepare_bi
integer(4) i,j
real(8) g,xx,yy,xx1,yy1
real(8) rrr
real(8) sred,cos45,eps,epsg
gs%const%epsMinDistBound=(rr(2)-rr(1))*d5
write(nfile,"('ZONE T=""area', A1, '_', i0, '_', i0, '"", I=', i0, ', J=', i0, ', F=POINT')") nchar,ic,jc,nnr+1,nng+1
cos45=dcos(pi*0.25d0)
eps=1.0d-4
epsg=1.0d-3
do i=1,nng+1
  g=bndg(i)
  do j=1,nnr+1
    rrr=sred(rr(1),rr(nnr+1),bndrv_partcl(i),bndrv(i),rr(j))
    xx=mcell%x(ic)+rrr*dcos(g)
 	  yy=mcell%y(jc)+rrr*dsin(g)
    xx1=xx
    yy1=yy
    mode=1
    if (j==1) mode=2
    if (j==nnr+1) then
      mode=0
      !чуть расширяем, чтобы ячейки перекрывались
		  if (dcos(g)>cos45-epsg) xx1=xx1+eps
		  if (dcos(g)<-cos45+epsg) xx1=xx1-eps
		  if (dsin(g)>cos45-epsg) yy1=yy1+eps
		  if (dsin(g)<-cos45+epsg) yy1=yy1-eps
    endif
    WRITE(nfile,"(F13.9, ' ', F13.9)",advance="no") xx1,yy1
    call write_point(nfile,xx,yy,rrr,g,mode)
  enddo
enddo
end

subroutine dns_drw_one_fictcell(nfile,mcell,ic,jc,nx,ny,nchar,write_point)
!вывод в файл решения во ячейке без твердой частицы (квадратная ячейка)
use dnsmod_base
use pgmod
integer(4) nfile      !номер файла
type(mcelltype) mcell !информация об ЭЯ
integer(4) ic,jc      !индекс ячейки по x и y
integer(4) nx        !количество отрезков по x
integer(4) ny        !количество отрезков по y
character(1) nchar    !буква (индекс) в названии зоны
!write_point(nfile,xx,yy,r,g,mode) - процедура для вывода значения в точке !!!r,g тут фиктивные
integer(4) i,j,mode,modex
real(8) xx,yy,x0,y0,dx,dy,eps,xx1,yy1
eps=1.0d-4
write(nfile,"('ZONE T=""area', A1, '_', i0, '_', i0, '"", I=', i0, ', J=', i0, ', F=POINT')") nchar,ic,jc,ny+1,nx+1
x0=mcell%x0+(ic-1)*mcell%dx
y0=mcell%y0+(jc-1)*mcell%dy
dx=mcell%dx/nx
dy=mcell%dy/ny
gs%const%epsMinDistBound=min(dx,dy)*d5
do i=0,nx
  xx=x0+i*dx
  xx1=xx
  mode=1
  if (i==0) then
    xx1=xx-eps
    mode=3
  endif
  if (i==nx) then
    xx1=xx+eps
    mode=3
  endif
  modex=mode
  do j=0,ny
    mode=modex
    yy=y0+j*dy
    yy1=yy
    if (j==0) then
      yy1=yy-eps
      mode=3
    endif
    if (j==ny) then
      yy1=yy+eps
      mode=3
    endif
    WRITE(nfile,"(F13.9, ' ', F13.9)",advance="no") xx1,yy1
    call write_point(nfile,xx,yy,d0,d0,mode)
  enddo
enddo
end

subroutine dns_prepare_drw_one_cell_r(mcell,ds0,minnr,maxnr,nnr,rr,nmax)
use dnsmod_base
type(mcelltype) mcell
integer(4) nnr,nmax,minnr,maxnr
real(8) rr(nmax),ds0
call ga_get_rr(mcell%r,mcell%dx*d5,rr,ds0,nnr,nmax,minnr,maxnr)  
nnr=nnr-1 !число отрезков
end

subroutine dns_prepare_drw_one_cell_g(mcell,nng_pi,nng,bndg,bndrv,nmax)
use dnsmod_base
type(mcelltype) mcell
integer(4) nng,nmax,nng_pi,i
real(8) bndg(nmax),bndrv(nmax)
nng=nng_pi
call ga_init_vneshg(nng,bndg,bndrv,mcell%dy*d5,mcell%dx*d5,3)
do i=2,nng+1
  bndg(i+nng)=pi+bndg(i)
  bndrv(i+nng)=bndrv(i)
enddo
nng=2*nng
end

subroutine dns_write_point_psi_om_(nfile,xx,yy,g,mode,bi)
use pgmod
integer(4) nfile
real(8) xx,yy,g,get_fun_xy_solve_bi,f,fom
integer(4) mode
type(TBound_info) bi
call get_fun_xy_prepare_bi(xx,yy,mode,bi)
f=get_fun_xy_solve_bi(xx,yy,1,d0,d0,bi)
fom=-get_fun_xy_solve_bi(xx,yy,3,d0,d0,bi)
WRITE(nfile,"(2(E14.5))",advance="no") f,fom
return
g=g
end

subroutine dns_write_point_psi_om(nfile,xx,yy,g,mode)
use pgmod
integer(4) nfile
real(8) xx,yy,g
integer(4) mode
type(TBound_info) bi
call dns_write_point_psi_om_(nfile,xx,yy,g,mode,bi)
WRITE(nfile,*)
end

subroutine dns_write_point_psi_om_uv(nfile,xx,yy,g,mode)
use pgmod
integer(4) nfile
real(8) xx,yy,g,u,v,ff(2) !,get_fun_xy_solve_bi
integer(4) mode
type(TBound_info) bi
call dns_write_point_psi_om_(nfile,xx,yy,g,mode,bi)
!u=get_fun_xy_solve_bi(xx,yy,2,d0,d1,bi)
!v=-get_fun_xy_solve_bi(xx,yy,2,d1,d0,bi)
call get_fun_xy_solve_gradient_bi(xx,yy,5,d0,d0,bi,ff)
u=ff(2)
v=-ff(1)
WRITE(nfile,"(2(E14.5))") u,v
end

subroutine dns_write_mcell_in2(nfile,mcell,ii0,ii1,jj0,jj1,mode_in,write_point)
use dnsmod_base
integer(4) nfile      
type(mcelltype) mcell 
integer(4) ii0,ii1,jj0,jj1  
integer(4) mode_in
external dns_write_particle_empty,write_point
call dns_write_mcell_in(nfile,mcell,ii0,ii1,jj0,jj1,mode_in,write_point,dns_write_particle_empty,.false.,.false.)
end

subroutine dns_write_mcell_in(nfile,mcell,ii0,ii1,jj0,jj1,mode_in,write_point,write_particle,write_time,need_write_particle)
use dnsmod_base
integer(4) nfile      !номер файла
type(mcelltype) mcell !информация об ЭЯ
integer(4) ii0,ii1,jj0,jj1  !индекс левой нижней и правой верхней ячейки
                            !ii0 может быть <1 - если моделируется только половина пористой области, а нужно выгрузить всю
integer(4) mode_in    !0 - сетка только по угловым узлам между ячейками
                      !1 - сетка по угловым узлам между ячейками, по серединам сторон ячеек и по частицам
                      !2 - сетка по угловым узлам между ячейками и по частицам (требует триангуляции в техплоте)
!write_point(nfile,x,y,i1,j1,mode) - процедура для вывода значения в точке
                        !mode = 0 - обычная точка
                        !     = 1 - частица (i1,j1 содержат индексы частицы)
!write_particle(nfile,x,y,i,j) - процедура для вывода частиц
logical need_write_particle   !выводить зону POINTS
logical write_time    !выводить время на экран
character(75) str
real(8) x,y,x1
integer(4) i,j,nx,ny,nc,mode,i1,j1
nx=ii1-ii0+1
ny=jj1-jj0+1
if (allocated(mcell%g)) then
  nc=0
  do i=ii0,ii1
    i1=i
    if (i<1) i1=1-i
    do j=jj0,jj1
      if (mcell%g(i1,j)>0) nc=nc+1
    enddo
  enddo
else
  nc=nx*ny
endif
if (write_time) call init_gs_time
mode=0
if (mode_in==0) then
	!сетка только по угловым узлам между ячейками
  write(nfile,"('ZONE T=""in"", I=', i4, ', J=', i4, ', F=POINT')") ny+1,nx+1
  do i=ii0-1,ii1
    x=mcell%x0+i*mcell%dx
    i1=i
    if (i1<ii0) i1=ii0
    if (i1>ii1) i1=ii1
    if (i1<1) i1=1-i1
    if (write_time) then
      write(str,"('Write file area_i. ', I0, ' (from ', I0, ')')") i-ii0+2, nx+1
      call write_by_time(str)
    endif
    do j=jj0-1,jj1
      y=mcell%y0+j*mcell%dy
      j1=j
      if (j1<1) j1=1
      if (j1>mcell%ny) j1=mcell%ny
      call write_point(nfile,x,y,i1,j1,mode)
    enddo
  enddo
elseif (mode_in==1) then
	!сетка по угловым узлам между ячейками, по серединам сторон ячеек и по частицам
  write(nfile,"('ZONE T=""in"", I=', i4, ', J=', i4, ', F=POINT')") ny*2+1,nx*2+1
  do i=0,nx*2
    x=mcell%x0+(i+(ii0-1)*2)*d5*mcell%dx
    i1=i/2+ii0
    if (i1<ii0) i1=ii0
    if (i1>ii1) i1=ii1
    if (i1<1) i1=1-i1
    if (write_time) then
      write(str,"('Write file area_i. ', I0, ' (from ', I0, ')')") i+1, nx*2+1
      call write_by_time(str)
    endif
    do j=0,ny*2
      y=mcell%y0+(j+(jj0-1)*2)*d5*mcell%dy
      j1=j/2+jj0
      if (j1<1) j1=1
      if (j1>mcell%ny) j1=mcell%ny
      mode=0
      if (mod(i,2)==1.and.mod(j,2)==1) mode=1 !частица
      call write_point(nfile,x,y,i1,j1,mode)
    enddo
	enddo
else
	!сетка по угловым узлам между ячейками и по частицам
	!требует триангуляции в техплоте
	write(nfile,"('ZONE T=""in"", I=', i4, ', F=POINT')") (ny+1)*(nx+1)+nc
  mode=0
  do i=ii0-1,ii1
    x=mcell%x0+i*mcell%dx
    i1=i
    if (i1<ii0) i1=ii0
    if (i1>ii1) i1=ii1
    if (i1<1) i1=1-i1
    if (write_time) then
      write(str,"('Write file area_i. ', I0, ' (from ', I0, ')')") i-ii0+2, nx+1
      call write_by_time(str)
    endif
    do j=jj0-1,jj1
      y=mcell%y0+j*mcell%dy
      j1=j
      if (j1<1) j1=1
      if (j1>mcell%ny) j1=mcell%ny
      call write_point(nfile,x,y,i1,j1,mode)
    enddo
  enddo
  mode=1
	do i=ii0,ii1
    i1=i
    if (i<1) i1=1-i
    do j=jj0,jj1
      if (allocated(mcell%g)) then
        if (mcell%g(i1,j)<=0) cycle
      endif
      x=mcell%x(i1)
      if (i<1) x=d2*mcell%x0-x
      call write_point(nfile,x1,mcell%y(j),i1,j,mode)
    enddo
  enddo
endif
if (need_write_particle) then
  write(nfile,"('ZONE T=""points"", I=', i4, ', F=POINT')") nc
  do i=ii0,ii1
    i1=i
    if (i<1) i1=1-i
    do j=jj0,jj1
      if (allocated(mcell%g)) then
        if (mcell%g(i1,j)<=0) cycle
      endif
      x=mcell%x(i1)
      if (i<1) x=d2*mcell%x0-x
      call write_particle(nfile,x,mcell%y(j),i1,j)
    enddo
  enddo
endif
if (write_time) call print_total_time
end

subroutine dns_write_particle_empty(nfile,x,y,i,j)
integer(4) nfile,i,j
real(8) x,y
return
nfile=nfile
i=i
j=j
x=x
y=y
end
program main
  use mod
  h=5.0d0
  call main_1(3)
  pause

end
subroutine main_1(g)
  use mod
  integer(4) g
  real(8) fg,psi_numerical
  select case (g)
  case (1)
    call pg_start
    call pg_allocate_problems(1)
    call pg_bind_problem(1)
    call pg_allocate_domains(1)
    call pg_bind_domain(1)
    call pg_set_domain_equation(3)
    call pg_allocate_bounds(1)
    call pg_bind_bound(1)
    call init_geom
    call pg_geom_postprocessor
    call init_gu(1)
    call ga_drw_trmesh(1)
    call pg_get_matrix
    call closing(1)
    call pg_solve
    call pg_get_psioml
    print *, "psi=",psi_numerical(d1,d1)
    !print *, "psi=",psi_numerical(dcos(pi*d5),dsin(pi*d5))
    call draw
    call print_psi_omega
    call pg_finish
    
  case (2)
    call pg_start
    call pg_allocate_problems(1)
    call pg_bind_problem(1)
    call pg_allocate_domains(1)
    call pg_bind_domain(1)
    call pg_set_domain_equation(21)
    call pg_allocate_bounds(1)
    call pg_bind_bound(1)
    call init_geom
    call pg_geom_postprocessor
    call init_gu(2)
    
    call pg_allocate_area(1)
    call pg_bind_areapart(1)
    call init_Mesh
    call pg_areageom_postprocessor
    
    call ga_drw_trmesh(1)
    
    call pg_allocate_area_gu
    call init_Meshval(1)
    call pg_get_matrix
    call closing(2)
    call pg_solve
    call pg_get_psioml
    print *, "fg=",fg(d2,d2)
    print *, "psi=",psi_numerical(d2,d2)
    call draw
    call print_psi_omega
    call print_fi_const(pi*0.5,1)
    call pg_finish
  case (3)
    call pg_start
    call pg_allocate_problems(1)
    call pg_bind_problem(1)
    call pg_allocate_domains(1)
    call pg_bind_domain(1)
    call pg_set_domain_equation(21)
    call pg_allocate_bounds(1)
    call pg_bind_bound(1)
    call init_geom
    call pg_geom_postprocessor
    call init_gu(2)
    
    call pg_allocate_area(1)
    call pg_bind_areapart(1)
    call init_Mesh_2
    call pg_areageom_postprocessor
    call ga_drw_trmesh(1)
    
    call pg_allocate_area_gu
    call init_Meshval(g)
    call pg_get_matrix
    call closing(g)
    call pg_solve
    call pg_get_psioml
    call draw
    call print_psi_omega
    call pg_finish
  end select 
end 
subroutine print_fi_const(fi,g1)
  use mod
  integer(4) ng,nr,i,g1
  real(8) x,y,r,g,pg_get_fun_xy,om,psi,dom,pg_f_psiomlxy,dpsi,fi
  open(1,File='data_fi_const.dat')
  WRITE(1,*) 'x y psi dpsi om dom'
  nr=30
  g=fi
  do i = 1, nr+1
    r=d1+(i-d1)*(h-d1)/nr
    select case (g1)
    
    case (2)
      if ( r>=R_1 .and. r<=R_2 ) then
        call pg_bind_domain(1)
        call pg_bind_bound(1)
      else if ( r>R_2 .and. r<=h ) then
        call pg_bind_domain(2)
        call pg_bind_bound(1)
      end if
    end select
    
    x=r*dcos(g)
    y=r*dsin(g)
    psi=pg_get_fun_xy(x,y,1,d0,d0,2)
    om=-pg_get_fun_xy(x,y,3,d0,d0,2)
    dom=-pg_get_fun_xy(x,y,4,d0,d1,2)
    dpsi=pg_get_fun_xy(x,y,2,d0,d1,2)
    WRITE(1,"(E13.5, ' ', E13.5, ' ',E15.5, ' ', E15.5, ' ', E13.5, ' ', E13.5)") x,y,psi,dpsi,om,dom
  end do
  close(1)
end
subroutine init_geom

  use mod
  integer(4) nj
  real(8) ds,ds2
  nj=100
  ds=pi/nj
  ds2=ds*d2
  call pg_allocate_boundlines(4)
  call pg_init_boundline_geomlineds2(1,ds,ds2,d1,d0,h,d0)
  call pg_init_boundline_geomcircleds(2,ds2,d0,d0,h,d0,pi)
  call pg_init_boundline_geomlineds2(3,ds2,ds,-h,d0,-d1,d0)
  call pg_init_boundline_geomcircleds(4,ds,d0,d0,d1,pi,d0)
  end

subroutine init_gu(g)
  use mod
  integer(4) g
  select case (g)
  case(1)
    !integer(4) n,pg_get_int
    !real(8), allocatable :: T(:)
    !real(8) T1(200)
    call pg_allocate_bound_gu
    call pg_allocate_constvalind(1)
    call pg_set_constvala(1,1,1)

    call pg_bind_boundline(1)
    call pg_init_boundline_gu_val_const(1,1,d0)
    call pg_init_boundline_gu_val_const(2,3,d0)

    call pg_bind_boundline(2)
    !n=pg_get_int(2,1)
    !allocate(T(n))
    !call pg_get_array_real(2,4,T,n)
    !call pg_init_boundline_gu_val(1,1,T)
    !call pg_get_array_real(2,8,T,n)
    !call pg_init_boundline_gu_val(1,4,-dsin(T-pi5)/10)
    !T=d0
    !call pg_init_boundline_gu_val(1,4,T)
    call pg_init_boundline_gu_val_const(1,3,d0)

    call pg_bind_boundline(3)
    call pg_init_boundline_gu_val_const(1,1,d0)
    call pg_init_boundline_gu_val_const(2,3,d0)

    call pg_bind_boundline(4)
    call pg_init_boundline_gu_val_const(1,1,d0)
    call pg_init_boundline_gu_val_const(2,2,d0)
    !deallocate(T)
  case(2)
    !integer(4) n,pg_get_int
    !real(8), allocatable :: T(:)
    !real(8) T1(200)
    call pg_allocate_bound_gu
    call pg_allocate_constvalind(1)
    call pg_set_constvala(1,1,1)

    call pg_bind_boundline(1)
    call pg_init_boundline_gu_val_const(1,1,d0)
    call pg_init_boundline_gu_val_const(2,3,d0)

    call pg_bind_boundline(2)
    !n=pg_get_int(2,1)
    !allocate(T(n))
    !call pg_get_array_real(2,4,T,n)
    !call pg_init_boundline_gu_val(1,1,T)
    !call pg_get_array_real(2,8,T,n)
    !call pg_init_boundline_gu_val(1,4,-dsin(T-pi5)/10)
    !T=d0
    !call pg_init_boundline_gu_val(1,4,T)
    call pg_init_boundline_gu_val_const(1,3,d0)

    call pg_bind_boundline(3)
    call pg_init_boundline_gu_val_const(1,1,d0)
    call pg_init_boundline_gu_val_const(2,3,d0)

    call pg_bind_boundline(4)
    call pg_init_boundline_gu_val_const(1,1,d0)
    call pg_init_boundline_gu_val_const(2,2,d0)
    !deallocate(T)
  end select
  end

subroutine print_psi_omega
  use mod
  integer(4) ng,i
  real(8) x,y,r,pg_get_fun_xy,om,psi,dom,g,pg_f_psiomlxy,dpsi
  open (2, FILE='om.dat')
  ng=30
  do i=1,ng+1
      g=(i-d1)*pi/ng
      !r=d1+(i-d1)*(h-d1)/ng
      r=d1
      x=r*dcos(g)
      y=r*dsin(g)
      psi=pg_get_fun_xy(x,y,1,d0,d0,2)
      if(r==d1) then
          dpsi=-pg_f_psiomlxy(x,y,2,0)
      else
          dpsi=pg_f_psiomlxy(x,y,2,0)
      end if
      !dpsi=pg_get_fun_xy(x+eps,y+eps,2,d0,d0,mode)
      om=-pg_get_fun_xy(x,y,3,d0,d0,2)
      if(r==d1) then
        dom=pg_f_psiomlxy(x,y,4,0)
      else 
        dom=-pg_f_psiomlxy(x,y,4,0)
      end if
      WRITE(2,"(E13.5, ' ', E13.5, ' ',E15.5, ' ', E15.5, ' ', E13.5, ' ', E13.5)") x,y,psi,dpsi,om,dom
  end do
  close(2)

  open (2, FILE='om_h.dat')
  do i=1,ng+1
      g=(i-d1)*pi/ng
      !r=d1+(i-d1)*(h-d1)/ng
      r=h
      x=r*dcos(g)
      y=r*dsin(g)
      psi=pg_get_fun_xy(x,y,1,d0,d0,2)
      if(r==d1) then
          dpsi=-pg_f_psiomlxy(x,y,2,0)
      else
          dpsi=pg_f_psiomlxy(x,y,2,0)
      end if
      !dpsi=pg_get_fun_xy(x+eps,y+eps,2,d0,d0,mode)
      om=-pg_get_fun_xy(x,y,3,d0,d0,2)
      if(r==d1) then
        dom=pg_f_psiomlxy(x,y,4,0)
      else 
        dom=-pg_f_psiomlxy(x,y,4,0)
      end if
      WRITE(2,"(E13.5, ' ', E13.5, ' ',E15.5, ' ', E15.5, ' ', E13.5, ' ', E13.5)") x,y,psi,dpsi,om,dom
  end do
  close(2)
  end subroutine print_psi_omega

subroutine draw
  use mod
  integer(4) nr,ng,i,j,mode
  real(8) g,r,x,y,psi,Vx,Vy,om,pg_get_fun_xy
  OPEN (1,FILE='data.dat')
  nr=20 !����� ����� �� r
  ng=30 !����� ����� �� gamma
  write(1,*) 'TITLE = "velosity"'
  write(1,*) 'VARIABLES = "X", "Y", "psi", "Vx", "Vy", "om"'
  write(1,"('ZONE T=""area"", I=', i0, ', J=', i0, ', F=POINT')") nr+1,ng+1
  do i=1,ng+1
    g=(i-d1)*pi/ng
    do j=1,nr+1
      r=d1+(j-d1)*(h-d1)/nr
      x=r*dcos(g)
      y=r*dsin(g)
      if (i==1.or.i==ng+1.or.j==1.or.j==nr+1) then
        mode=2
      else
        mode=1
      endif
      psi=pg_get_fun_xy(x,y,1,d0,d0,mode)
      Vx=pg_get_fun_xy(x,y,2,d0,d1,mode)
      Vy=-pg_get_fun_xy(x,y,2,d1,d0,mode)
      om=-pg_get_fun_xy(x,y,3,d0,d0,mode)
      if (i==1.and.j==nr+1) then
        Vx=pg_get_fun_xy(x-eps,y,2,d0,d1,mode)
        Vy=-pg_get_fun_xy(x-eps,y,2,d1,d0,mode)
      endif
      if (i==ng+1.and.j==nr+1) then
        Vx=pg_get_fun_xy(x+eps,y,2,d0,d1,mode)
        Vy=-pg_get_fun_xy(x+eps,y,2,d1,d0,mode)
      endif
      if (i==ng+1.and.j==1) then
        Vx=pg_get_fun_xy(x-eps,y,2,d0,d1,mode)
        Vy=-pg_get_fun_xy(x-eps,y,2,d1,d0,mode)
      endif
      if (i==1.and.j==1) then
        Vx=pg_get_fun_xy(x+eps,y,2,d0,d1,mode)
        Vy=-pg_get_fun_xy(x+eps,y,2,d1,d0,mode)
      endif
      !WRITE(1,"(F9.5, ' ', F9.5, ' ', F9.5, ' ', F9.5, ' ', F9.5, ' ', F9.5)") x,y,psi,Vx,Vy,om
      WRITE(1,"(E13.5, ' ', E13.5, ' ', E13.5, ' ', E13.5, ' ', E13.5, ' ', E13.5)") x,y,psi,Vx,Vy,om
    enddo
  enddo
  close(1)
  end

subroutine closing(g)
  use mod
  integer(4) i
  real(8), allocatable :: eq(:)
  real(8) b,dtt,tt
  integer(4) g
  select case (g)
  case(1,3)
    allocate(eq(gs%m%nx))
    b=d0
    eq=d0
    call pg_bind_boundline(2)
    dtt=pi/gsbndl%npanel
    do i=gsbndl%i_begin,gsbndl%i_end 
        tt=dtt*(i-gsbndl%i_begin+d5)
        eq(gsbnd%psiind(i,1))=d1
        eq(gs%constvalinds(1))=-dsin(tt)
        call pg_add_closing_eq_to_area(1,eq,b)
        eq(gsbnd%psiind(i,1))=d0
        eq(gs%constvalinds(1))=d0
    enddo
    do i=gsbndl%i_begin,gsbndl%i_end 
        eq(gsbnd%psiind(i,4))=d1
    enddo
    b=-d2/dtt
    call pg_add_closing_eq_to_area(1,eq,b)
    deallocate(eq)
  case(2)
    allocate(eq(gs%m%nx))
    b=d0
    eq=d0
    call pg_bind_boundline(2)
    dtt=pi/gsbndl%npanel
    do i=gsbndl%i_begin,gsbndl%i_end 
        tt=dtt*(i-gsbndl%i_begin+d5)
        eq(gsbnd%psiind(i,1))=d1
        eq(gs%constvalinds(1))=-dsin(tt)
        call pg_add_closing_eq_to_area(1,eq,b)
        eq(gsbnd%psiind(i,1))=d0
        eq(gs%constvalinds(1))=d0
    enddo
    do i=gsbndl%i_begin,gsbndl%i_end 
        eq(gsbnd%psiind(i,4))=d1
    enddo
    b=-d2/(dtt*h**3)
    call pg_add_closing_eq_to_area(1,eq,b)
    deallocate(eq)
  
  end select
  end
subroutine init_Mesh
  use mod
  real(8), allocatable :: rr1(:), gg1(:)
  integer(4) n1,ng1,pg_get_int,i
  call pg_bind_boundline(1)
  n1=pg_get_int(2,1)+1
  allocate (rr1(n1))
  call pg_get_array_real(2,1,rr1,n1)
  call pg_bind_boundline(2)
  ng1=pg_get_int(2,1)+1
  allocate (gg1(ng1))
  !call pg_get_array_real(2,8,gg1,ng1)
  do i=1, ng1
    gg1(i)=(i-d1)*pi/(ng1-1)
  enddo
  call ga_init_mesh_rcell_quads(rr1,gg1,n1,ng1)
  deallocate(rr1,gg1)
  end
subroutine init_Mesh_2
  use mod
  real(8), allocatable :: rr1(:), gg1(:)
  integer(4) ng1,pg_get_int,i
  ng1=40
  allocate (rr1(2))
  allocate (gg1(ng1))
  do i=1, ng1
    gg1(i)=(i-d1)*pi/(ng1-1)
  enddo
  rr1=[d1+((h-d1)/2)-eps,1+((h-d1)/2)+eps]
  call ga_init_mesh_rcell_quads(rr1,gg1,2,ng1)
  deallocate(rr1,gg1)
  end
subroutine init_Meshval(g1)
  use mod
  real(8), allocatable:: g(:),x(:),y(:)
  integer(4) ntr, pg_get_int,i
  integer(4) g1
  real(8) fg
  !real(8) fgsin
  !real(8) f_3
  ntr=pg_get_int(4,1)
  allocate (g(ntr),x(ntr),y(ntr))
  call pg_get_array_real(4,1,x,ntr)
  call pg_get_array_real(4,2,y,ntr)
  select case (g1)
  case (1)
    do i=1,ntr
      g(i)=fg(x(i),y(i))
    enddo
  case (2)
    do i=1,ntr
      g(i)=fgsin(x(i),y(i))
    enddo
  case (3)
    do i=1,ntr
      g(i)=f_3(x(i),y(i))
    enddo 
  end select 
  call pg_init_area_gu(g,1)
  deallocate (g,x,y)
  end
function fg(x,y)
  use mod
  real(8) fg,x,y
  fg=(3*y)/((x**d2+y**d2)**(2.5d0))
  end
function fgsin(x,y)
  use mod
  real(8) fg,x,y
  fg=-(pi*y*dcos((pi*(1-(x**2+y**2)**0.5))/(2*(h-1)))/(2*(h-1)*(x**2+y**2)**0.5))
  end
function psi_numerical(x,y)
  use mod
  real(8) psi_numerical,x,y,pg_get_fun_xy
  psi_numerical=(pg_get_fun_xy(x-eps,y,3,d0,d0,0)-2*pg_get_fun_xy(x,y,3,d0,d0,0)+pg_get_fun_xy(x+eps,y,3,d0,d0,0))/eps**2+(pg_get_fun_xy(x,y-eps,3,d0,d0,0)-2*pg_get_fun_xy(x,y,3,d0,d0,0)+pg_get_fun_xy(x,y+eps,3,d0,d0,0))/eps**2
  end
function f_3(x,y)
  use mod
  real(8) a,x,y
  a=-1/(2*eps)
  !b=d1+(3d0-eps)/(2*eps)
  f_3=-a*y/sqrt(x**2+y**2)
  end
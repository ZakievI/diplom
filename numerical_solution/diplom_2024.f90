subroutine main_1(g)
  use mod
  integer(4) g
  call pg_start
  call pg_allocate_problems(1)
  call pg_bind_problem(1)
  call pg_allocate_domains(1)
  call pg_bind_domain(1)
  select case (g)
  case (1)
    call pg_set_domain_equation(3)
  case (2,3,4)
    call pg_set_domain_equation(21)
  end select
  call pg_allocate_bounds(1)
  call pg_bind_bound(1)
  call init_geom
  call pg_geom_postprocessor
  call init_gu
  
  select case (g)
  case (2,3,4)
    call pg_allocate_area(1)
    call pg_bind_areapart(1)
    call init_Mesh(g)
    call pg_areageom_postprocessor
  end select
  
  
  select case (g)
  case (2,3,4)
    call pg_allocate_area_gu
    call init_Meshval(g)
  end select
  call ga_drw_trmesh(1)
  call pg_get_matrix
  call closing(g)
  call pg_solve
  call pg_get_psioml
  !call compute_pressure(2)
  !call draw(1,2)
  call error(g,1)
  !call print_fi_const(pi*0.5,1)
  call pg_finish
end 
subroutine init_geom
  use mod
  call pg_allocate_boundlines(4)
  call pg_init_boundline_geomlineds2(1,ds,ds2,d1,d0,h,d0)
  call pg_init_boundline_geomcircleds(2,ds2,d0,d0,h,d0,pi)
  call pg_init_boundline_geomlineds2(3,ds2,ds,-h,d0,-d1,d0)
  call pg_init_boundline_geomcircleds(4,ds,d0,d0,d1,pi,d0)
  end
subroutine init_gu
  use mod
    call pg_allocate_bound_gu
    call pg_allocate_constvalind(1)
    call pg_set_constvala(1,1,1)

    call pg_bind_boundline(1)
    call pg_init_boundline_gu_val_const(1,1,d0)
    call pg_init_boundline_gu_val_const(2,3,d0)

    call pg_bind_boundline(2)
    call pg_init_boundline_gu_val_const(1,3,d0)

    call pg_bind_boundline(3)
    call pg_init_boundline_gu_val_const(1,1,d0)
    call pg_init_boundline_gu_val_const(2,3,d0)

    call pg_bind_boundline(4)
    call pg_init_boundline_gu_val_const(1,1,d0)
    call pg_init_boundline_gu_val_const(2,2,d0)
    !deallocate(T)
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
      call pg_bind_domain(1)
      call pg_bind_bound(1)
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
      call pg_bind_domain(2)
      call pg_bind_bound(1)
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
subroutine error(number_task, bound)
  use mod
  integer(4) nr,ng,i,j,number_task,bound
  real(8) g,r,x,y,psi,Vx,Vy,om,pg_get_fun_xy,E_psi,E_u,E_v,E_om,max_E_psi,max_E_u,max_E_v,max_E_om
  OPEN (1,FILE='error_.dat') 
  write(1,*) 'TITLE = "velosity"'
  write(1,*) 'VARIABLES = "X", "Y", "E_psi", "E_u", "E_v", "E_om"'
  nr=30 
  ng=50
  write(1,"('ZONE T=""area"", I=', i0, ', J=', i0, ', F=POINT')") nr+1,ng+1
  E_psi=d0
  E_u=d0
  E_v=d0
  E_om=d0
  max_E_psi=d0
  max_E_u=d0
  max_E_v=d0
  max_E_om=d0
  do i=1,ng+1
    g=(i-d1)*pi/ng
    do j=1,nr+1
      r=d1+(j-d1)*(h-R_1)/nr
      select case (bound)
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
      psi=pg_get_fun_xy(x,y,1,d0,d0,0)
      Vx=pg_get_fun_xy(x,y,2,d0,d1,0)
      Vy=-pg_get_fun_xy(x,y,2,d1,d0,0)
      om=-pg_get_fun_xy(x,y,3,d0,d0,0)
      if (i==1.and.j==nr+1) then
        Vx=pg_get_fun_xy(x-eps,y,2,d0,d1,0)
        Vy=-pg_get_fun_xy(x-eps,y,2,d1,d0,0)
      endif
      if (i==ng+1.and.j==nr+1) then
        Vx=pg_get_fun_xy(x+eps,y,2,d0,d1,0)
        Vy=-pg_get_fun_xy(x+eps,y,2,d1,d0,0)
      endif
      if (i==ng+1.and.j==1) then
        Vx=pg_get_fun_xy(x-eps,y,2,d0,d1,0)
        Vy=-pg_get_fun_xy(x-eps,y,2,d1,d0,0)
      endif
      if (i==1.and.j==1) then
        Vx=pg_get_fun_xy(x+eps,y,2,d0,d1,0)
        Vy=-pg_get_fun_xy(x+eps,y,2,d1,d0,0)
      endif
      select case (number_task)     
    case(1)
        E_psi=((r**2-R_1**2)*(2*h**2+r**2-R_1**2)+(dlog(R_1/r))*(2*h*r)**2)*dsin(g)/(16*r)+psi
        E_u=-((x)**6+49*y**2-52*y**4+3*y**6+(48+5*y**2)*(x)**4+(-49-4*y**2+7*y**4)*(x)**2-50*((x)**2+y**2)**2*dlog((x)**2+y**2))/(16*((x)**2+y**2)**2)-Vx
        E_v=((x)*y*(49+(x)**4-50*y**2+y**4+2*(x)**2*(-25+y**2)))/(8*((x)**2+y**2)**2)-Vy
        if (i==1.and.j==nr+1) then
            !E_u=-((x-eps)**6+49*y**2-52*y**4+3*y**6+(48+5*y**2)*(x-eps)**4+(-49-4*y**2+7*y**4)*(x-eps)**2-50*((x-eps)**2+y**2)**2*dlog((x-eps)**2+y**2))/(16*((x-eps)**2+y**2)**2)-Vx
            !E_v=((x-eps)*y*(49+(x-eps)**4-50*y**2+y**4+2*(x-eps)**2*(-25+y**2)))/(8*((x-eps)**2+y**2)**2)-Vy
            E_u=d0
            E_v=d0
        else if (i==ng+1.and.j==nr+1) then
            !E_u=-((x+eps)**6+49*y**2-52*y**4+3*y**6+(48+5*y**2)*(x+eps)**4+(-49-4*y**2+7*y**4)*(x+eps)**2-50*((x+eps)**2+y**2)**2*dlog((x+eps)**2+y**2))/(16*((x+eps)**2+y**2)**2)-Vx
            !E_v=((x+eps)*y*(49+(x+eps)**4-50*y**2+y**4+2*(x+eps)**2*(-25+y**2)))/(8*((x+eps)**2+y**2)**2)-Vy
            E_u=d0
            E_v=d0
        else if (i==ng+1.and.j==1) then
            E_u=-((x-eps)**6+49*y**2-52*y**4+3*y**6+(48+5*y**2)*(x-eps)**4+(-49-4*y**2+7*y**4)*(x-eps)**2-50*((x-eps)**2+y**2)**2*dlog((x-eps)**2+y**2))/(16*((x-eps)**2+y**2)**2)-Vx
            E_v=((x-eps)*y*(49+(x-eps)**4-50*y**2+y**4+2*(x-eps)**2*(-25+y**2)))/(8*((x-eps)**2+y**2)**2)-Vy
        else if (i==1.and.j==1) then
            E_u=-((x+eps)**6+49*y**2-52*y**4+3*y**6+(48+5*y**2)*(x+eps)**4+(-49-4*y**2+7*y**4)*(x+eps)**2-50*((x+eps)**2+y**2)**2*dlog((x+eps)**2+y**2))/(16*((x+eps)**2+y**2)**2)-Vx
            E_v=((x+eps)*y*(49+(x+eps)**4-50*y**2+y**4+2*(x+eps)**2*(-25+y**2)))/(8*((x+eps)**2+y**2)**2)-Vy
        endif
        E_om=(y*(-25+x**2+y**2))/(2*(x**2+y**2))-om
    case(2)
        !E_psi=(-1-(-2*h*R_1+R_1**2)/(4*h**r)-(r*dlog(r)/(2*h))-(r*(-2*h-R_1-2*R_1*dlog(R_1)))/(4*h*R_1))*dsin(g)-psi
        E_psi=(((r - R_1)*(2*h*(r - R_1) + R_1*(r + R_1)) - 2*(r**2)*R_1*dlog(r))*dsin(g)/(4*h*r*R_1))-psi
        E_u=(11*(x)**4+9*y**2*(-1+y**2)+x**2*(9+20*y**2-20*dsqrt(x**2+y**2))-(x**2+y**2)**2*dlog((x)**2+y**2))/(20*(x**2+y**2)**2)-Vx
        E_v=(x*y*(9+x**2+y**2-10*dsqrt(x**2+y**2)))/(10*(x**2+y**2)**2)-Vy
        if (i==1.and.j==nr+1) then
            !E_u=-((x-eps)**6+49*y**2-52*y**4+3*y**6+(48+5*y**2)*(x-eps)**4+(-49-4*y**2+7*y**4)*(x-eps)**2-50*((x-eps)**2+y**2)**2*dlog((x-eps)**2+y**2))/(16*((x-eps)**2+y**2)**2)-Vx
            !E_v=((x-eps)*y*(49+(x-eps)**4-50*y**2+y**4+2*(x-eps)**2*(-25+y**2)))/(8*((x-eps)**2+y**2)**2)-Vy
            E_u=d0
            E_v=d0
        else if (i==ng+1.and.j==nr+1) then
            !E_u=-((x+eps)**6+49*y**2-52*y**4+3*y**6+(48+5*y**2)*(x+eps)**4+(-49-4*y**2+7*y**4)*(x+eps)**2-50*((x+eps)**2+y**2)**2*dlog((x+eps)**2+y**2))/(16*((x+eps)**2+y**2)**2)-Vx
            !E_v=((x+eps)*y*(49+(x+eps)**4-50*y**2+y**4+2*(x+eps)**2*(-25+y**2)))/(8*((x+eps)**2+y**2)**2)-Vy
            E_u=d0
            E_v=d0
        else if (i==ng+1.and.j==1) then
            E_u=(11*(x-eps)**4+9*y**2*(-1+y**2)+(x-eps)**2*(9+20*y**2-20*dsqrt((x-eps)**2+y**2))-((x-eps)**2+y**2)**2*dlog((x-eps)**2+y**2))/(20*((x-eps)**2+y**2)**2)-Vx
            E_v=((x-eps)*y*(9+(x-eps)**2+y**2-10*dsqrt((x-eps)**2+y**2)))/(10*((x-eps)*2+y**2)**2)-Vy
        else if (i==1.and.j==1) then
            E_u=(11*(x+eps)**4+9*y**2*(-1+y**2)+(x+eps)**2*(9+20*y**2-20*dsqrt((x+eps)**2+y**2))-((x+eps)**2+y**2)**2*dlog((x+eps)**2+y**2))/(20*((x+eps)**2+y**2)**2)-Vx
            E_v=((x+eps)*y*(9+(x+eps)**2+y**2-10*dsqrt((x+eps)**2+y**2)))/(10*((x+eps)*2+y**2)**2)-Vy
        endif
        E_om=(y*(-5+dsqrt(x**2+y**2)))/(5*dsqrt((x**2+y**2)**3))-om
    case(3) 
        if(r<=R_2)then
            E_psi=(17/(16*r)-r-(r**3)/16+(9*r*dlog(r))/4)*dsin(g)-psi
            E_u=-(x**6 + 17*y**2 - 20*y**4 + 3*y**6 + x**4*(16 + 5*y**2) + x**2*(-17 - 4*y**2 + 7*y**4) - 18*(x**2 + y**2)**2*dlog(x**2 + y**2))/(16*(x**2 + y**2)**2)-Vx
            E_v=(x*y*(17 + x**4 - 18*y**2 + y**4 + 2*x**2*(-9 + y**2)))/(8*(x**2 + y**2)**2)-Vy
            !E_om=((17/(16*r) - r - (r**3)/16 + 9/4*r*dlog(r))*dsin(g))/(r**2) - ((17/(8*r**3) + 9/(4*r) - (3*r)/8)*r*dsin(g) + (5/4 - 17/(16*r**2) - (3*r**2)/16 + (9*dlog(r))/4)*dsin(g))/r-om
            E_om=((-9 + r**2)*dsin(g))/(2*r)-om
            if (i==1.and.r==R_2) then
                E_u=d0
                E_v=d0
            end if
        else
            E_psi=(-4/r+(-4+9*dlog(3*d1))*r/4)*dsin(g)-psi
            E_u=(x**4*(-4 + 9*dlog(3*d1)) + 2*x**2*(-8 + y**2*(-4 + 9*dlog(3*d1))) + y**2*(16 + y**2*(-4 + 9*dlog(3*d1))))/(4*(x**2 + y**2)**2)-Vx
            E_v=-((8*x*y)/(x**2 + y**2)**2)-Vy
            !E_om=((-(4/r) + 1/4*r*(-4 + 9*dlog(3*d1)))*dsin(g))/r**2 - (-((8*dsin(g))/r**2) + (4/r**2 + 1/4*(-4 + 9*dlog(3*d1)))*dsin(g))/r-om
            E_om=-om
        end if
    end select
    if (max_E_psi <dabs(E_psi)) then
    max_E_psi=dabs(E_psi)
    end if
    if (max_E_u <dabs(E_u)) then
    max_E_u=dabs(E_u)
    end if
    if (max_E_v <dabs(E_v)) then
    max_E_v=dabs(E_v)
    end if
    if (max_E_om <dabs(E_om)) then
    max_E_om=dabs(E_om)
    end if
    WRITE(1,"(E13.5, ' ', E13.5, ' ', E13.5, ' ', E13.5, ' ', E13.5, ' ', E13.5)") x,y,dabs(E_psi),dabs(E_u),dabs(E_v),dabs(E_om)
    enddo
  enddo
  close(1)
  print *,"max_E_psi=:", max_E_psi,"max_E_u=:", max_E_u,"max_E_v=:", max_E_v,"max_E_om=:", max_E_om
  end
subroutine draw(g1, number_task)
  use mod
  
  integer(4) nr,ng,i,j,mode,g1,number_task
  real(8) g,r,x,y,psi,Vx,Vy,om,pg_get_fun_xy,E,max_E
  OPEN (1,FILE='data.dat')
  nr=30 
  ng=50
  E=0
  max_e=0
  mode=0
  write(1,*) 'TITLE = "velosity"'
  write(1,*) 'VARIABLES = "X", "Y", "psi", "Vx", "Vy", "om"'
  write(1,"('ZONE T=""area"", I=', i0, ', J=', i0, ', F=POINT')") nr+1,ng+1
  do i=1,ng+1
    g=(i-d1)*pi/ng
    do j=1,nr+1
      r=d1+(j-d1)*(h-R_1)/nr
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
      psi=pg_get_fun_xy(x,y,1,d0,d0,mode)
      Vx=pg_get_fun_xy(x,y,2,d0,d1,mode)
      Vy=-pg_get_fun_xy(x,y,2,d1,d0,mode)
      om=-pg_get_fun_xy(x,y,3,d0,d0,mode)
      select case (number_task)     
    case(1)
        E=((r**2-R_1**2)*(2*h**2+r**2-R_1**2)+(dlog(R_1/r))*(2*h*r)**2)*dsin(g)/(16*r)+psi
    case(2)
        E=(-1-(-2*h*R_1+R_1**2)/(4*h**r)-(r*dlog(r)/(2*h))-(r*(-2*h-R_1-2*R_1*dlog(R_1)))/(4*h*R_1))*dsin(g)-psi
    case(3) 
        if ( r>=R_1 .and. r<=R_2 ) then
            E=(17/(16*r)-r-(r**3)/16+(9*r*dlog(r))/4)*dsin(g)-psi
        else
            E=(-4/r+(-4+9*dlog(3*d1))*r/4)*dsin(g)-psi
        end if
    end select
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
      if (max_E <dabs(E)) then
          max_E=dabs(E)
      end if
      !WRITE(1,"(F9.5, ' ', F9.5, ' ', F9.5, ' ', F9.5, ' ', F9.5, ' ', F9.5)") x,y,psi,Vx,Vy,om
      WRITE(1,"(E13.5, ' ', E13.5, ' ', E13.5, ' ', E13.5, ' ', E13.5, ' ', E13.5)") x,y,psi,Vx,Vy,om
    enddo
  enddo
  close(1)
  print *,"nj=",nj,"  max_E=:", max_E
  !pause
  end
subroutine draw_square(g1)
  use mod
  integer(4) nr,ng,i,j,mode,k1,g1
  real(8) g,r,x,y,psi,Vx,Vy,om,pg_get_fun_xy,bndg(200),bndrv(200)
  ng=60 !число ячеек по gamma
  call ga_init_vneshg(ng,bndg,bndrv,h,h,8)

  OPEN (1,FILE='data.dat')
  nr=20 !число ячеек по r

  write(1,*) 'TITLE = "velosity"'
  write(1,*) 'VARIABLES = "X", "Y", "psi", "Vx", "Vy", "om"'
  write(1,"('ZONE T=""area"", I=', i0, ', J=', i0, ', F=POINT')") (nr+1),2*(ng+1)
  do k1=1,2
    call pg_bind_domain(k1)
  do i=1,ng+1
    select case(g1)    
    case(1)
        g=bndg(i)+pi*(k1-d1)
    case(2)
        g=bndg(i)+pi*(k1-d1)+pi5
    end select
    do j=1,nr+1
      r=R_1+(j-d1)*(bndrv(i)-R_1)/nr
      x=h+r*dcos(g)
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
      if (i==1.and.j==nr+1.and.k1==1) then
        Vx=pg_get_fun_xy(x-eps,y,2,d0,d1,mode)
        Vy=-pg_get_fun_xy(x-eps,y,2,d1,d0,mode)
      else if (i==ng+1.and.j==nr+1.and.k1==1) then
        Vx=pg_get_fun_xy(x-eps,y,2,d0,d1,mode)
        Vy=-pg_get_fun_xy(x-eps,y,2,d1,d0,mode)
      else if (i==ng+1.and.j==1.and.k1==1) then
        Vx=pg_get_fun_xy(x-eps,y,2,d0,d1,mode)
        Vy=-pg_get_fun_xy(x-eps,y,2,d1,d0,mode)
      else if (i==1.and.j==1.and.k1==1) then
        Vx=pg_get_fun_xy(x-eps,y,2,d0,d1,mode)
        Vy=-pg_get_fun_xy(x-eps,y,2,d1,d0,mode)
        !!!!
      else if (i==1.and.j==nr+1.and.k1==2) then
        Vx=pg_get_fun_xy(x+eps,y,2,d0,d1,mode)
        Vy=-pg_get_fun_xy(x+eps,y,2,d1,d0,mode)
      else if (i==ng+1.and.j==nr+1.and.k1==2) then
        Vx=pg_get_fun_xy(x+eps,y,2,d0,d1,mode)
        Vy=-pg_get_fun_xy(x+eps,y,2,d1,d0,mode)
      else if (i==ng+1.and.j==1.and.k1==2) then
        Vx=pg_get_fun_xy(x+eps,y,2,d0,d1,mode)
        Vy=-pg_get_fun_xy(x+eps,y,2,d1,d0,mode)
      else if (i==1.and.j==1.and.k1==2) then
        Vx=pg_get_fun_xy(x+eps,y,2,d0,d1,mode)
        Vy=-pg_get_fun_xy(x+eps,y,2,d1,d0,mode)
      endif
      WRITE(1,"(F13.5, ' ', F13.5, ' ', F13.5, ' ', F13.5, ' ', F13.5, ' ', F13.5)") x,y,psi,Vx,Vy,om
    enddo
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
  
  allocate(eq(gs%m%nx))
  b=d0
  eq=d0
  call pg_bind_boundline(2)
  dtt=pi/gsbndl%npanel
  do i=gsbndl%i_begin,gsbndl%i_end 
      tt=dtt*(i-gsbndl%i_begin+d5)
      eq(gsbnd%psiind(i,1))=d1
      eq(gs%constvalinds(1))=-dsin(tt)
      call pg_add_closing_eq_to_area(gsarea%i,eq,b) 
      eq(gsbnd%psiind(i,1))=d0
      eq(gs%constvalinds(1))=d0
  enddo
  do i=gsbndl%i_begin,gsbndl%i_end 
      eq(gsbnd%psiind(i,4))=d1
  enddo
  select case (g)
  case(1)
      b=-d2/dtt
  case(2)
      b=-d2/(dtt*h**3)
  case (3,4)
      b=0
  end select
  call pg_add_closing_eq_to_area(gsarea%i,eq,b)  
  deallocate(eq)
  
 end
subroutine init_Mesh(g)
  use mod
  integer(4) g
  real(8), allocatable :: rr1(:), gg1(:)
  integer(4) n1,pg_get_int,i,ng1
  select case (g)
  case (2,4)
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
  case (3)
    allocate (rr1(2))
    allocate (gg1(nj))
    do i=1, nj
      gg1(i)=(i-d1)*pi/(nj-1)
    enddo
    rr1=[R_2-eps_zd3,R_2+eps_zd3]
    call ga_init_mesh_rcell_quads(rr1,gg1,2,nj)
  end select 
  deallocate(rr1,gg1)
  end
subroutine init_Meshval(g1)
  use mod
  real(8), allocatable:: g(:),x(:),y(:)
  integer(4) ntr, pg_get_int,i
  integer(4) g1
  ntr=pg_get_int(4,1)
  allocate (g(ntr),x(ntr),y(ntr))
  call pg_get_array_real(4,1,x,ntr)
  call pg_get_array_real(4,2,y,ntr)
  select case (g1)
  case (1)
    do i=1,ntr
      g(i)=-(pi*y(i)*dcos((pi*(1-(x(i)**2+y(i)**2)**0.5))/(2*(h-1)))/(2*(h-1)*(x(i)**2+y(i)**2)**0.5))
    enddo
  case (2)
    do i=1,ntr
      g(i)=(3*y(i))/((x(i)**d2+y(i)**d2)**(2.5d0))
    enddo
  case (3)
    cof_a=-d1/(2*eps_zd3)
    do i=1,ntr
      g(i)=-cof_a*y(i)/sqrt(x(i)**2+y(i)**2)
    enddo 
  case(4)
      g(i)=-2*(x(i)*y(i))/(x(i)**2+y(i)**2)
  end select 
  call pg_init_area_gu(g,1)
  deallocate (g,x,y)
 end
subroutine print_fi_const(fi,g1)
  use mod
  integer(4) nr,i,g1
  real(8) x,y,r,g,pg_get_fun_xy,om,psi,dom,dpsi,fi
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
subroutine main_2
  use mod
  call pg_start
  call pg_allocate_problems(1)
  call pg_bind_problem(1)
  call pg_allocate_domains(2)
  
  
  call pg_bind_domain(1)
  call pg_set_domain_equation(3)
  call pg_allocate_bounds(1)
  call pg_bind_bound(1)
  call init_geom_1(1)
  call pg_geom_postprocessor
  call init_gu_1(1)
  
  
  call pg_bind_domain(2)
  call pg_set_domain_equation(3)
  call pg_allocate_bounds(1)
  call pg_bind_bound(1)
  call init_geom_2(1)
  call pg_geom_postprocessor
  call init_gu_2(1)
  
  
  
  call ga_drw_trmesh(0)
  call pg_get_matrix
  call closing_2(1)
  call pg_bind_domain(2)
  call pg_bind_bound(1)
  call closing(3)
  call pg_solve
  call pg_get_psioml
  call draw(2,3)
  call error(3,2)
  !call print_fi_const(pi*0.5,2)
  !call print_psi_omega
  call pg_finish
  end
subroutine init_geom_1(g)
    !1-цилиндр в цилиндре 
    !2- цилиндр в кдвадратной ячейке с горизонтальным разрезом 
    !3-цилиндр в кдвадратной ячейке с вертиканым разрезом 
    !4-цилиндр в кдвадратной ячейке с разрезом под углом
    use mod
    integer(4) g
    select case (g)
        case(1)
          call pg_allocate_boundlines(4)
          call pg_init_boundline_geomlineds2(1,ds,(ds+ds2)/2,R_1,d0,R_2,d0)
          call pg_init_boundline_geomcircleds(2,(ds+ds2)/2,d0,d0,R_2,d0,pi)
          call pg_init_boundline_geomlineds2(3,(ds+ds2)/2,ds,-R_2,d0,-R_1,d0)
          call pg_init_boundline_geomcircleds(4,ds,d0,d0,R_1,pi,d0)
        case (2)
          call pg_allocate_boundlines(6)
          call pg_init_boundline_geomlineds(1,ds2,d0,d0,h-R_1,d0)
          call pg_init_boundline_geomcircleds(2,ds,h,d0,R_1,pi,d0) 
          call pg_init_boundline_geomlineds(3,ds2,h+R_1,d0,2*h,d0)
          call pg_init_boundline_geomlineds(4,ds2,2*h,d0,2*h,h)
          call pg_init_boundline_geomlineds(5,ds2,2*h,h,d0,h)
          call pg_init_boundline_geomlineds(6,ds2,d0,h,d0,d0)
       case (3)
          call pg_allocate_boundlines(6)
          call pg_init_boundline_geomlineds(1,ds2,d0,-h,h,-h)
          call pg_init_boundline_geomlineds2(2,ds2,ds,h,-h,h,-R_1)  
          call pg_init_boundline_geomcircleds(3,ds,h,d0,R_1,3*pi/2,pi/2) 
          call pg_init_boundline_geomlineds2(4,ds,ds2,h,R_1,h,h) 
          call pg_init_boundline_geomlineds(5,ds2,h,h,d0,h)
          call pg_init_boundline_geomlineds(6,ds2,d0,h,d0,-h)
      case (4)
          call pg_allocate_boundlines(6)
          call pg_init_boundline_geomlineds(1,ds2,d0,-h,d0,h)
          call pg_init_boundline_geomlineds(2,ds2,d0,h,h-sin(pi/4)*R_1,R_1*sin(pi/4))
          call pg_init_boundline_geomcircleds(3,ds,h,d0,R_1,3*pi/4,7*pi/4) 
         call pg_init_boundline_geomlineds(2,ds2,h+sin(pi/4)*R_1,-R_1*sin(pi/4),2*h,-h)
          call pg_init_boundline_geomlineds(5,ds2,2*h,-h,d0,-h)
     end select
  end
subroutine init_geom_2(g)
    !1-цилиндр в цилиндре 
    !2-цилиндр в кдвадратной ячейке с горизонтальным разрезом 
    !3-цилиндр в кдвадратной ячейке с вертиканым разрезом
    !4-цилиндр в кдвадратной ячейке с разрезом под углом
    use mod
    integer(4) g
    select case (g)
        case(1)
          call pg_allocate_boundlines(4)
          call pg_init_boundline_geomlineds2(1,(ds+ds2)/2,ds2,R_2,d0,h,d0)
          call pg_init_boundline_geomcircleds(2,ds2,d0,d0,h,d0,pi)
          call pg_init_boundline_geomlineds2(3,ds2,(ds+ds2)/2,-h,d0,-R_2,d0)
          call pg_copy_boundline(4,1,1,1,2)
        case (2)
          call pg_allocate_boundlines(6)
          call pg_init_boundline_geomlineds(1,ds2,h-R_1,d0,d0,d0)
          call pg_init_boundline_geomlineds(2,ds2,d0,d0,d0,-h)
          call pg_init_boundline_geomlineds(3,ds2,d0,-h,2*h,-h)
          call pg_init_boundline_geomlineds(4,ds2,2*h,-h,2*h,d0)
          call pg_init_boundline_geomlineds(5,ds2,2*h,d0,h+R_1,d0)
          call pg_init_boundline_geomcircleds(6,ds,h,d0,R_1,2*pi,pi)
       case (3)
          call pg_allocate_boundlines(6)
          call pg_init_boundline_geomlineds(1,ds2,h,-h,2*h,-h)
          call pg_init_boundline_geomlineds(2,ds2,2*h,-h,2*h,h)
          call pg_init_boundline_geomlineds(3,ds2,2*h,h,h,h)
          call pg_init_boundline_geomlineds2(4,ds2,ds,h,h,h,R_1)
          call pg_init_boundline_geomcircleds(5,ds,h,d0,R_1,5*pi/2,3*pi/2) 
          call pg_init_boundline_geomlineds2(6,ds,ds2,h,-R_1,h,-h)
       case (4)
          call pg_allocate_boundlines(6)
          call pg_init_boundline_geomlineds(1,ds2,d0,h,2*h,h)
          call pg_init_boundline_geomlineds(2,ds2,2*h,h,2*h,-h)
          call pg_init_boundline_geomlineds(3,ds2,2*h,-h,h+sin(pi/4)*R_1,-R_1*sin(pi/4))
          call pg_init_boundline_geomcircleds(4,ds,h,d0,R_1,7*pi/4,11*pi/4) 
          call pg_init_boundline_geomlineds(5,ds2,h-sin(pi/4)*R_1,R_1*sin(pi/4),d0,h)
        end select 
   end
subroutine init_gu_1(g)
    !1-цилиндр в цилиндре, нижняя область с действием массовых сил
    !2-цилиндр в кдвадратной ячейке с горизонтальным разрезом, верхняя область без действия массовых сил
    !3-цилиндр в кдвадратной ячейке с вертиканым разрезом
    !4-цилиндр в кдвадратной ячейке с разрезом под углом
    use mod
    integer(4) g
    select case (g)
        case(1)
      call pg_allocate_bound_gu
  
      call pg_bind_boundline(1)
      call pg_init_boundline_gu_val_const(1,1,d0)
      call pg_init_boundline_gu_val_const(2,3,d0)
  
      call pg_bind_boundline(3)
      call pg_init_boundline_gu_val_const(1,1,d0)
      call pg_init_boundline_gu_val_const(2,3,d0)
  
      call pg_bind_boundline(4)
      call pg_init_boundline_gu_val_const(1,1,d0)
      call pg_init_boundline_gu_val_const(2,2,d0)
    case (2)
      call pg_allocate_bound_gu
      call pg_allocate_constvalind(2)
      call pg_set_constvala(1,1,1) 
      call pg_set_constvala(2,2,2) 

      call pg_bind_boundline(2)
      call pg_init_boundline_gu_val_const(1,1,d0)
      call pg_init_boundline_gu_val_const(2,2,d0)
  
      call pg_bind_boundline(4)
      call pg_init_boundline_gu_val_const(1,2,d0)
      call pg_init_boundline_gu_val_const(2,4,d0)
         
      call pg_bind_boundline(5)
      call pg_init_boundline_gu_constval(1,1,1)
      call pg_init_boundline_gu_val_const(2,3,d0)

      call pg_bind_boundline(6)
      call pg_init_boundline_gu_val_const(1,2,d0)
      call pg_init_boundline_gu_val_const(2,4,d0)
    case (3)
      call pg_allocate_bound_gu
      call pg_allocate_constvalind(2)
      call pg_set_constvala(1,1,1) 
      call pg_set_constvala(2,2,2) 

      call pg_bind_boundline(1)
      call pg_init_boundline_gu_constval(1,1,1)
      call pg_init_boundline_gu_val_const(2,3,d0)
  
      call pg_bind_boundline(3)
      call pg_init_boundline_gu_val_const(1,1,d0)
      call pg_init_boundline_gu_val_const(2,2,d0)
         
      call pg_bind_boundline(5)
      call pg_init_boundline_gu_constval(1,1,2)
      call pg_init_boundline_gu_val_const(2,3,d0)

      call pg_bind_boundline(6)
      call pg_init_boundline_gu_val_const(1,2,d0)
      call pg_init_boundline_gu_val_const(2,4,d0)

      end select 
    end
subroutine init_gu_2(g)
    !1-цилиндр в цилиндре, верхняя область без действия массовых сил
    !2-цилиндр в кдвадратной ячейке с горизонтальным разрезом, верхняя область без действия массовых сил
    !3-цилиндр в кдвадратной ячейке с вертиканым разрезом
    !4-цилиндр в кдвадратной ячейке с разрезом под углом
      use mod
      integer g
      select case (g)
      case(1)
        call pg_allocate_bound_gu
        call pg_allocate_constvalind(1)
        call pg_set_constvala(2,1,1)
    
        call pg_bind_boundline(1)
        call pg_init_boundline_gu_val_const(1,1,d0)
        call pg_init_boundline_gu_val_const(2,3,d0)
    
        call pg_bind_boundline(2)
        call pg_init_boundline_gu_val_const(1,3,d0)
    
        call pg_bind_boundline(3)
        call pg_init_boundline_gu_val_const(1,1,d0)
        call pg_init_boundline_gu_val_const(2,3,d0)
    
        call pg_bind_boundline(4)
        call pg_init_boundline_gu_gen(1,1,.false.,1,1,2,1,d0,[.true.],[1],[d1])
        call pg_init_boundline_gu_gen(2,3,.false.,1,1,2,1,d0,[.true.],[3],[d1])
        
        call pg_bind_domain(1)
        call pg_bind_bound(1)
        call pg_bind_boundline(2)
        call pg_init_boundline_gu_gen(1,2,.false.,2,1,4,1,d0,[.true.],[2],[-d1])
       case(2)
        call pg_allocate_bound_gu
      
        call pg_bind_boundline(1)
        call pg_init_boundline_gu_gen(1,1,.false.,1,1,1,1,d0,[.true.],[1],[d1])
        call pg_init_boundline_gu_gen(2,3,.false.,1,1,1,1,d0,[.true.],[3],[d1])
        
        call pg_bind_boundline(2)
        call pg_init_boundline_gu_val_const(1,2,d0)
        call pg_init_boundline_gu_val_const(2,4,d0)
        
        call pg_bind_boundline(3)
        call pg_init_boundline_gu_constval(1,1,2)
        call pg_init_boundline_gu_val_const(2,3,d0)
      
        call pg_bind_boundline(4)
        call pg_init_boundline_gu_val_const(1,2,d0)
        call pg_init_boundline_gu_val_const(2,4,d0)
           
        call pg_bind_boundline(5)
        call pg_init_boundline_gu_gen(1,1,.false.,1,1,3,1,d0,[.true.],[1],[d1])
        call pg_init_boundline_gu_gen(2,3,.false.,1,1,3,1,d0,[.true.],[3],[d1])
    
        call pg_bind_boundline(6)
        call pg_init_boundline_gu_val_const(1,1,d0)
        call pg_init_boundline_gu_val_const(2,2,d0)
        
        call pg_bind_domain(1)
        call pg_bind_bound(1)
        call pg_bind_boundline(1)
        call pg_init_boundline_gu_gen(1,2,.false.,2,1,1,1,d0,[.true.],[2],[-d1])
        call pg_bind_boundline(3)
        call pg_init_boundline_gu_gen(1,2,.false.,2,1,5,1,d0,[.true.],[2],[-d1])
       case(3)
        call pg_allocate_bound_gu
      
        call pg_bind_boundline(1)
        call pg_init_boundline_gu_constval(1,1,1)
        call pg_init_boundline_gu_val_const(2,3,d0)
        
        call pg_bind_boundline(2)
        call pg_init_boundline_gu_val_const(1,2,d0)
        call pg_init_boundline_gu_val_const(2,4,d0)
        
        call pg_bind_boundline(3)
        call pg_init_boundline_gu_constval(1,1,2)
        call pg_init_boundline_gu_val_const(2,3,d0)
      
        call pg_bind_boundline(4)
        call pg_init_boundline_gu_gen(1,1,.false.,1,1,4,1,d0,[.true.],[1],[d1])
        call pg_init_boundline_gu_gen(2,3,.false.,1,1,4,1,d0,[.true.],[3],[d1])
           
        call pg_bind_boundline(5)
        call pg_init_boundline_gu_val_const(1,1,d0)
        call pg_init_boundline_gu_val_const(2,2,d0)
        
        call pg_bind_boundline(6)
        call pg_init_boundline_gu_gen(1,1,.false.,1,1,2,1,d0,[.true.],[1],[d1])
        call pg_init_boundline_gu_gen(2,3,.false.,1,1,2,1,d0,[.true.],[3],[d1])
        
        call pg_bind_domain(1)
        call pg_bind_bound(1)
        call pg_bind_boundline(4)
        call pg_init_boundline_gu_gen(1,2,.false.,2,1,4,1,d0,[.true.],[2],[-d1])
        call pg_init_boundline_gu_gen(2,4,.false.,2,1,4,1,d0,[.true.],[4],[-d1])
        call pg_bind_boundline(2)
        call pg_init_boundline_gu_gen(1,2,.false.,2,1,6,1,d0,[.true.],[2],[-d1])
        call pg_init_boundline_gu_gen(2,4,.false.,2,1,6,1,d0,[.true.],[4],[-d1])
         end  select 
     end
subroutine closing_2(g)
    !1-цилиндр в цилиндре, нижняя область с действием массовых сил
    !2-цилиндр в кдвадратной ячейке с горизонтальным разрезом, верхняя область без действия массовых сил
    !3-цилиндр в кдвадратной ячейке с вертиканым разрезом
    !4-цилиндр в кдвадратной ячейке с разрезом под углом
  use mod 
  integer(4) i,j,k2,k1,i1
  real(8), allocatable :: eq(:)
  real(8) b,dtt,tt
  integer(4) g
 
  select case (g)
  case(1)
    allocate(eq(gs%m%nx))
    dtt=pi/gsbndl%npanel
    eq=d0
    do i=gsbndl%i_begin,gsbndl%i_end
        call pg_bind_domain(1)
        call pg_bind_bound(1)
        call pg_bind_boundline(2)
        i1=i-gsbndl%i_begin+1
        tt=dtt*(i-gsbndl%i_begin+d5)
        k1=gsbnd%psiind(i,4)
        eq(k1)=d1
        call pg_bind_domain(2)
        call pg_bind_bound(1)
        call pg_bind_boundline(4)
        j=gsbndl%i_end-i1+1
        k2=gsbnd%psiind(j,4)
        eq(k2)=d1
        b=-F_m*dsin(tt)
        call pg_add_closing_eq_to_area(1,eq,b)
        eq(k1)=d0
        eq(k2)=d0
    enddo
    deallocate(eq)
  case(2)
    call pg_bind_domain(1)
    call pg_bind_bound(1)
    call pg_bind_boundline(1)
    allocate(eq(gs%m%nx))
    eq=d0
    do i=gsbndl%i_begin,gsbndl%i_end
        call pg_bind_domain(1)
        call pg_bind_bound(1)
        call pg_bind_boundline(1)
        i1=i-gsbndl%i_begin+1
        k1=gsbnd%psiind(i,4)
        eq(k1)=d1
        call pg_bind_domain(2)
        call pg_bind_bound(1)
        call pg_bind_boundline(1)
        j=gsbndl%i_end-i1+1
        k2=gsbnd%psiind(j,4)
        eq(k2)=d1
        b=-F_m
        call pg_add_closing_eq_to_area(1,eq,b)
        eq(k1)=d0
        eq(k2)=d0
    enddo
    eq=d0
    call pg_bind_domain(1)
    call pg_bind_bound(1)
    call pg_bind_boundline(3)
    do i=gsbndl%i_begin,gsbndl%i_end
        call pg_bind_domain(1)
        call pg_bind_bound(1)
        call pg_bind_boundline(3)
        i1=i-gsbndl%i_begin+1
        k1=gsbnd%psiind(i,4)
        eq(k1)=d1
        call pg_bind_domain(2)
        call pg_bind_bound(1)
        call pg_bind_boundline(5)
        j=gsbndl%i_end-i1+1
        k2=gsbnd%psiind(j,4)
        eq(k2)=d1
        b=-F_m
        call pg_add_closing_eq_to_area(1,eq,b)
        eq(k1)=d0
        eq(k2)=d0
    enddo
    deallocate(eq)
    end select
end 
subroutine closing_3(g)
  use mod
  integer(4) i,nnx,g
  real(8), allocatable :: eq(:)
  real(8) b,dtt,dxx1
  nnx=gs%m%nx  
  allocate(eq(gs%m%nx))
  b=d0
  eq=d0
  select case (g)
  case (1)
       !замыкание для границы на цилиндре   
       call pg_bind_domain(1)
       call pg_bind_bound(1)
       call pg_bind_boundline(2)
       dtt=pi/gsbndl%npanel
       do i=gsbndl%i_begin,gsbndl%i_end 
           eq(gsbnd%psiind(i,4))=d1
       enddo
       call pg_bind_domain(2)
       call pg_bind_bound(1)
       call pg_bind_boundline(6)
       do i=gsbndl%i_begin,gsbndl%i_end 
           eq(gsbnd%psiind(i,4))=d1
       enddo
        b=-F_m*d2/dtt
       call pg_add_closing_eq_to_area(2,eq,b)
       
       ! замыкания на верхней границе 
       eq=d0
       call pg_bind_domain(1)
       call pg_bind_bound(1)
       call pg_bind_boundline(5)
       do i=gsbndl%i_begin,gsbndl%i_end 
           eq(gsbnd%psiind(i,4))=d1
       enddo
       b=0
       call pg_add_closing_eq_to_area(1,eq,b) 
       deallocate(eq)
  case (2)
       !замыкание для границы на цилиндре   
       call pg_bind_domain(1)
       call pg_bind_bound(1)
       call pg_bind_boundline(3)
       dtt=pi/gsbndl%npanel
       do i=gsbndl%i_begin,gsbndl%i_end 
           eq(gsbnd%psiind(i,4))=d1
       enddo
       call pg_bind_domain(2)
       call pg_bind_bound(1)
       call pg_bind_boundline(5)
       do i=gsbndl%i_begin,gsbndl%i_end 
           eq(gsbnd%psiind(i,4))=d1
       enddo
        b=d0 !!!!
       call pg_add_closing_eq_to_area(2,eq,b)
       
       ! замыкания на верхней границе 
       eq=d0
       call pg_bind_domain(1)
       call pg_bind_bound(1)
       call pg_bind_boundline(5)
       dxx1=h/gsbndl%npanel
       do i=gsbndl%i_begin,gsbndl%i_end 
           eq(gsbnd%psiind(i,4))=d1
       enddo
       call pg_bind_domain(2)
       call pg_bind_bound(1)
       call pg_bind_boundline(3)
       do i=gsbndl%i_begin,gsbndl%i_end 
           eq(gsbnd%psiind(i,4))=d1
       enddo
       b=-F_m*h/dxx1 !!!!
       call pg_add_closing_eq_to_area(1,eq,b) 
       deallocate(eq)          
      end select
 end
subroutine main_3(g)
  use mod
  !g-1 -квадрат с кругом постановка №1,
  !  3 -квадрат с кругом постановка №2
  integer(4) g
  !real(8) fg,psi_numerical
  call pg_start
  call pg_allocate_problems(1)
  call pg_bind_problem(1)
  call pg_allocate_domains(2)
  
  
  call pg_bind_domain(1)
  call pg_set_domain_equation(3)
  call pg_allocate_bounds(1)
  call pg_bind_bound(1)
  call init_geom_1(g+1)!2-квадрат с кругом постановка №1, 3-квадрат с кругом постановка №2
  call pg_geom_postprocessor
  call init_gu_1(g+1)!2-квадрат с кругом постановка №1, 3-квадрат с кругом постановка №2
  
  
  call pg_bind_domain(2)
  call pg_set_domain_equation(3)
  call pg_allocate_bounds(1)
  call pg_bind_bound(1)
  call init_geom_2(g+1)!2-квадрат с кругом постановка №1, 3-квадрат с кругом постановка №2
  call pg_geom_postprocessor
  call init_gu_2(g+1)!2-квадрат с кругом постановка №1, 3-квадрат с кругом постановка №2
  
  call ga_drw_trmesh(0)
  call pg_get_matrix
  call closing_2(g+1)!1-цилиндр в цилиндре, нижняя область с действием массовых сил 2-квадрат с кругом постановка №1, 3-квадрат с кругом постановка №2
  call closing_3(g)!1-квадрат с кругом постановка №1, 2-квадрат с кругом постановка №2
  call pg_solve
  call pg_get_psioml
  call compute_pressure(1)
  call draw_square(g)!1-квадрат с кругом постановка №1, 2-квадрат с кругом постановка №2
  !call print_fi_const(pi*0.5,2)
  !call print_psi_omega
  call pg_finish
end
!subroutine compute_pressure
!use mod
!real(8) dfdl
!integer(4) i,n
!external dfdl
!gsareapart=>ap
!call ga_init_mesh_quadcell_ds(d0,d1,d0,d1,0.1d0,4)
!!call ga_init_geom_circle_symmetr(d0,d0,h,.false.,40)
!call ga_drw_trmesh(0)
!end
subroutine init_geom_1_cheking
  use mod
  call pg_allocate_boundlines(4)
  call pg_init_boundline_geomlineds2(1,ds,ds,d0,d0,h,d0)
  call pg_init_boundline_geomlineds2(2,ds,ds,h,d0,h,2*h)
  call pg_init_boundline_geomlineds2(3,ds,ds,h,2*h,d0,2*h)
  call pg_init_boundline_geomlineds2(4,ds,ds,d0,2*h,d0,d0)
end 
subroutine init_geom_2_cheking
  use mod
  call pg_allocate_boundlines(4)
  call pg_init_boundline_geomlineds2(1,ds,ds,h,d0,2*h,d0)
  call pg_init_boundline_geomlineds2(2,ds,ds,2*h,d0,2*h,2*h)
  call pg_init_boundline_geomlineds2(3,ds,ds,2*h,2*h,h,2*h)
  call pg_init_boundline_geomlineds2(4,ds,ds,h,2*h,h,d0)
end 
subroutine init_gu_1_cheking
 use mod
 call pg_allocate_bound_gu
 call pg_allocate_constvalind(1)
 call pg_set_constvala(1,1,1) 
 
 call pg_bind_boundline(1)
 call pg_init_boundline_gu_val_const(1,1,d0)
 call pg_init_boundline_gu_val_const(2,3,d0)
 
 call pg_bind_boundline(3)
 call pg_init_boundline_gu_constval(1,1,1)
 call pg_init_boundline_gu_val_const(2,3,d0)
 
 call pg_bind_boundline(4)
 call pg_init_boundline_gu_val_const(1,2,d0)
 call pg_init_boundline_gu_val_const(2,4,d0)
end
subroutine init_gu_2_cheking
 use mod
 call pg_allocate_bound_gu
 
 call pg_bind_boundline(1)
 call pg_init_boundline_gu_val_const(1,1,d0)
 call pg_init_boundline_gu_val_const(2,3,d0)
  
 call pg_bind_boundline(2)
 call pg_init_boundline_gu_val_const(1,2,d0)
 call pg_init_boundline_gu_val_const(2,4,d0)
 
 call pg_bind_boundline(3)
 call pg_init_boundline_gu_constval(1,1,1)
 call pg_init_boundline_gu_val_const(2,3,d0)
 
 call pg_bind_boundline(4)
 call pg_init_boundline_gu_gen(1,1,.false.,1,1,2,1,d0,[.true.],[1],[d1])
 call pg_init_boundline_gu_gen(2,3,.false.,1,1,2,1,d0,[.true.],[3],[d1])
        
 call pg_bind_domain(1)
 call pg_bind_bound(1)
 call pg_bind_boundline(2)
 call pg_init_boundline_gu_gen(1,2,.false.,2,1,4,1,d0,[.true.],[2],[-d1])
 call pg_init_boundline_gu_gen(2,4,.false.,2,1,4,1,d0,[.true.],[4],[-d1])
end
subroutine closing_cheking
 use mod
 integer(4) i
 real(8), allocatable :: eq(:)
 real(8) b,dxx
 b=d0
 call pg_bind_domain(1)
 call pg_bind_bound(1)
 call pg_bind_boundline(2)
 allocate(eq(gs%m%nx))
 eq=d0
     call pg_bind_domain(1)
     call pg_bind_bound(1)
     call pg_bind_boundline(3)
     dxx=h/gsbndl%npanel
     do i=gsbndl%i_begin,gsbndl%i_end 
         eq(gsbnd%psiind(i,4))=d1
     enddo
     call pg_bind_domain(2)
     call pg_bind_bound(1)
     call pg_bind_boundline(3)
     do i=gsbndl%i_begin,gsbndl%i_end 
         eq(gsbnd%psiind(i,4))=d1
     enddo
     b=-F_m*h/dxx !!!!
     call pg_add_closing_eq_to_area(1,eq,b)      
 deallocate(eq)
end
subroutine main_3_cheking
  use mod
  call pg_start
  call pg_allocate_problems(1)
  call pg_bind_problem(1)
  call pg_allocate_domains(2)
  
  
  call pg_bind_domain(1)
  call pg_set_domain_equation(3)
  call pg_allocate_bounds(1)
  call pg_bind_bound(1)
  call init_geom_1_cheking
  call pg_geom_postprocessor
  call init_gu_1_cheking
  
  call ga_drw_trmesh(0)
  
  call pg_bind_domain(2)
  call pg_set_domain_equation(3)
  call pg_allocate_bounds(1)
  call pg_bind_bound(1)
  call init_geom_2_cheking
  call pg_geom_postprocessor
  call init_gu_2_cheking
  
  call ga_drw_trmesh(0)
  call pg_get_matrix
  call closing_cheking
  call pg_solve
  call pg_get_psioml
  call draw_square_cheking
  call pg_finish
end
subroutine draw_square_cheking
  use mod
  integer(4) nx,ny,i,j,k1 !,mode
  real(8) dx,dy,x,y,psi,Vx,Vy,om,pg_get_fun_xy

  OPEN (1,FILE='data.dat')
  nx=10
  ny=20
  dx=h/nx
  dy=2*h/ny
  write(1,*) 'TITLE = "velosity"'
  write(1,*) 'VARIABLES = "X", "Y", "psi", "Vx", "Vy", "om"'
  write(1,"('ZONE T=""area"", I=', i0, ', J=', i0, ', F=POINT')") (ny+1),2*(nx+1)
  do k1=1,2
    call pg_bind_domain(k1)
  do i=1,nx+1 
     x=h*(k1-d1)+dx*(i-d1)
    do j=1,ny+1
      y=dy*(j-d1)
      psi=pg_get_fun_xy(x,y,1,d0,d0,0)
      Vx=pg_get_fun_xy(x,y,2,d0,d1,0)
      Vy=-pg_get_fun_xy(x,y,2,d1,d0,0)
      om=-pg_get_fun_xy(x,y,3,d0,d0,0)
      if (i==1.and.j==ny+1) then
        Vx=pg_get_fun_xy(x-eps,y,2,d0,d1,0)
        Vy=-pg_get_fun_xy(x-eps,y,2,d1,d0,0)
      endif
      if (i==nx+1.and.j==ny+1) then
        Vx=pg_get_fun_xy(x+eps,y,2,d0,d1,0)
        Vy=-pg_get_fun_xy(x+eps,y,2,d1,d0,0)
      endif
      if (i==nx+1.and.j==1) then
        Vx=pg_get_fun_xy(x-eps,y,2,d0,d1,0)
        Vy=-pg_get_fun_xy(x-eps,y,2,d1,d0,0)
      endif
      if (i==1.and.j==1) then
        Vx=pg_get_fun_xy(x+eps,y,2,d0,d1,0)
        Vy=-pg_get_fun_xy(x+eps,y,2,d1,d0,0)
      endif
      WRITE(1,"(E13.5, ' ', E13.5, ' ', E13.5, ' ', E13.5, ' ', E13.5, ' ', E13.5)") x,y,psi,Vx,Vy,om
    enddo
  enddo
  enddo
  close(1)
end
subroutine build_mesh(g)
 use mod
 real(8) bndg(200),bndrv(200),rr(200)
 integer(4) ng,nr,ngg,g
 ng=60
 nr=20
 bndg=d0
 bndrv=d0
 rr=d0
 call pg_allocate_area(1)
 call pg_bind_areapart(1)
 select case(g)
 case(1)
    call ga_init_vneshg(ng,bndg,bndrv,h,h,2)
    do ngg=ng+2,2*ng+1
     bndg(ngg)=bndg(ngg-ng)+pi
     bndrv(ngg)=bndrv(ngg-ng)
    end do
    call sred_aray(d1,h,rr,nr+1)
    call ga_init_mesh_rcell_quads(rr,bndg,nr+1,2*(ng+1))
    call ga_mesh_square(d1,h,2*ng,bndrv,bndg)
    gsareapart%zm=gsareapart%zm+h
 case(2)
    call ga_init_vneshg(ng,bndg,bndrv,h,h,2)
    call sred_aray(d1,h,rr,nr+1)
    call ga_init_mesh_rcell_quads(rr,bndg,nr+1,(ng+1))
 end select 

end
subroutine compute_pressure(g)
!g=
!   1-задача с квадратом с вырезанным кругом по центру, координаты квадрата (-h,0)*(h,2*h)
!   2-задача с полукругом, центр в точке (0,0)
!
 use mod
 integer(4) i
 real(8) dfdl,dfdl2
 integer(4) ng,nr,g
 external dfdl,dfdl2
 call build_mesh(g)
 call pg_areageom_postprocessor
 call ga_drw_trmesh(0)
 allocate(ff(gsareapart%n),ffknow(gsareapart%n),err(gsareapart%n))
 ff=d0
 ffknow=.false.
 select case (g)
 case(1) 
    do i=1,gsareapart%n
        if (dabs(dreal(gsareapart%zm(i)))<eps .or. dabs(dreal(gsareapart%zm(i))-2*h)<eps )then
         ffknow(i)=.true.
         ff(i)=d0
        endif
    end do
    call intf_dxdy(gsareapart%zm,gsareapart%trm,gsareapart%n,gsareapart%ntr,gsareapart%npe,ff,ffknow,dfdl,d0)
 case(2)
    do i=1,gsareapart%n
        if (dabs(dreal(gsareapart%zm(i))-h)<eps .or. dabs(dreal(gsareapart%zm(i))+h)<eps )then
         ffknow(i)=.true.
         ff(i)=d0
        endif
    end do
    call intf_dxdy(gsareapart%zm,gsareapart%trm,gsareapart%n,gsareapart%ntr,gsareapart%npe,ff,ffknow,dfdl2,d0)
 end select 
 OPEN (1,FILE='pressure.dat')
 write(1,*) 'TITLE = "velosity"'
 write(1,*) 'VARIABLES = "X", "Y", "pressure"'
 ng=60
 nr=20
 select case(g)
 case(1)
    write(1,"('ZONE T=""area"", I=', i0, ', J=', i0, ', F=POINT')") 2*(ng+1),nr+1
    do i=1,gsareapart%n
    WRITE(1,"(F13.5, ' ', F13.5, ' ', F13.5, ' ', F13.5, ' ', F13.5, ' ', F13.5)") dreal(gsareapart%zm(i)),dimag(gsareapart%zm(i)),ff(i)
    enddo
 case(2)
    write(1,"('ZONE T=""area"", I=', i0, ', J=', i0, ', F=POINT')") ng+1,nr+1
    do i=1,gsareapart%n
    WRITE(1,"(F13.5, ' ', F13.5, ' ', F13.5, ' ', F13.5, ' ', F13.5, ' ', F13.5)") dreal(gsareapart%zm(i)),dimag(gsareapart%zm(i)),ff(i)
    enddo
 end select 
end

function dfdl(x,y,xder,yder)
use mod
real(8) grad(2)
real(8) dfdl,x,y,xder,yder,der1,der2,f
complex(8) lder
lder=dcmplx(xder,yder)
lder=lder/cdabs(lder)
der1=dreal(lder)
der2=dimag(lder)
f=d0
if(x<=h) then
!if(x>=d0.and.x<=h)then
 call pg_bind_domain(1) 
 call pg_get_fun_xy_gradient(x,y,6,0,grad)
 if (dabs(der1)>eps) f=(grad(2)+F_m)*der1
 if (dabs(der2)>eps) f=f-grad(1)*der2
else !if(x>h.and.x<=2*h) then
 call pg_bind_domain(2) 
 call pg_get_fun_xy_gradient(x,y,6,0,grad)
 if (dabs(der1)>eps) f=grad(2)*der1
 if (dabs(der2)>eps) f=f-grad(1)*der2
end if 
dfdl=f
end

function dfdl2(x,y,xder,yder)
use mod
real(8) grad(2)
real(8) dfdl2,x,y,xder,yder,der1,der2,f
complex(8) lder
lder=dcmplx(xder,yder)
lder=lder/cdabs(lder)
der1=dreal(lder)
der2=dimag(lder)
call pg_bind_domain(1) 
call pg_get_fun_xy_gradient(x,y,6,0,grad)
f=d0
if (dabs(der1)>eps) f=(grad(2)+F_m)*der1
if (dabs(der2)>eps) f=f-grad(1)*der2
dfdl2=f
end
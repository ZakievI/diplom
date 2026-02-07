program main
  use mod
  call main_1()
  pause
end
subroutine main_1()
    use mod
    H1              =5d0
    L1              =10*d1
    nj              =100
    ds              =pi/nj
    ds2             =ds
    ds_pg           =ds
    dlt             = 0.01d0
    call pg_start
    call pg_allocate_problems(1)
    call pg_bind_problem(1)
    call pg_allocate_domains(1)
    call pg_bind_domain(1)
    call pg_set_domain_equation(3)
    call pg_allocate_bounds(1)
    call pg_bind_bound(1)
    call init_geom(1)
    call pg_geom_postprocessor
    call init_gu
    call build_mesh()
    call pg_areageom_postprocessor
    call ga_drw_trmesh(1)
    call pg_get_matrix
    call pg_solve
    call pg_get_psioml
    call draw_square
    gs_test_point_near_bound = .false.
    call init_zaplat
    !call build_bound()
    call build_curve()
    call find_derivative()
    !call checking_derivatives()
    !call find_concentration()
    
    call find_Jacobian()
    call find_concentration_by_Jacobian()
    call build_mesh_1()
    call filling_in_the_cells()
    !call draw_mesh(3)
    !call draw_Curves(2)
    call draw_mesh_with_cell(1)
    
    !call build_time_isolines()
    !call solve
    call pg_finish
    !call testing()
    
end

subroutine draw_square
  use mod
  integer(4) nr,ng,i,j,mode
  real(8) g,r,x,y,psi,Vx,Vy,om,pg_get_fun_xy,bndg(200),bndrv(200)
  ng=60 !число ячеек по gamma
  call ga_init_vneshg(ng,bndg,bndrv,H1,L1/2,8)

  OPEN (1,FILE='data.dat')
  nr=20 !число ячеек по r

  write(1,*) 'TITLE = "velosity"'
  write(1,*) 'VARIABLES = "X", "Y", "psi", "Vx", "Vy", "om"'
  write(1,"('ZONE T=""area"", I=', i0, ', J=', i0, ', F=POINT')") (nr+1),(ng+1)
  do i=1,ng+1
    g=bndg(i)
    do j=1,nr+1
      r=d1+(j-d1)*(bndrv(i)-d1)/nr
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
      else if (i==ng+1.and.j==nr+1) then
        Vx=pg_get_fun_xy(x-eps,y,2,d0,d1,mode)
        Vy=-pg_get_fun_xy(x-eps,y,2,d1,d0,mode)
      else if (i==ng+1.and.j==1) then
        Vx=pg_get_fun_xy(x-eps,y,2,d0,d1,mode)
        Vy=-pg_get_fun_xy(x-eps,y,2,d1,d0,mode)
      else if (i==1.and.j==1) then
        Vx=pg_get_fun_xy(x-eps,y,2,d0,d1,mode)
        Vy=-pg_get_fun_xy(x-eps,y,2,d1,d0,mode)
      endif
      WRITE(1,"(F13.5, ' ', F13.5, ' ', F13.5, ' ', F13.5, ' ', F13.5, ' ', F13.5)") x,y,psi,Vx,Vy,om
  enddo
  enddo
  close(1)
end


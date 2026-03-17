program main
  use mod
  call main_1()
  pause
end
subroutine main_1()
    use mod
    integer(4) :: max_iter
    real(8) :: max_delta, tol
    real(8), allocatable :: prev_force(:)
    H1              = 5d0
    L1              = 16*d1
    nj              = 100
    ds              = pi/nj
    ds2             = ds
    ds_pg           = ds
    dlt             = 0.01d0
    tol             = 1e-8
    max_iter        = 20
    max_delta       = tol + d1
    iteration       = 0
    call init_error_history()
    do while (iteration < max_iter .and. max_delta > tol)
        write(*,"('Iteration=',i0)") iteration
        call run_iteration(prev_force, max_delta)
        write(*,"('Error=',ES12.5)") max_delta
        call append_error_history(iteration, max_delta)
        iteration = iteration + 1
    end do
  end

subroutine init_error_history()
    integer(4), parameter :: unit_error = 97
    open(unit_error, file='data/error_history.dat', status='replace', action='write')
    write(unit_error, '(A)') '# iteration error'
    close(unit_error)
end

subroutine append_error_history(iteration_value, error_value)
    integer(4), intent(in) :: iteration_value
    real(8), intent(in) :: error_value
    integer(4), parameter :: unit_error = 97
    open(unit_error, file='data/error_history.dat', status='old', position='append', action='write')
    write(unit_error, '(I0,1X,ES24.16)') iteration_value, error_value
    close(unit_error)
end

subroutine run_iteration(prev_force, max_delta)
    use mod
    real(8), intent(out) :: max_delta
    real(8), allocatable :: prev_force(:)
    integer(4) :: n_
    real(8), allocatable :: x_(:), y_(:)
    logical :: curve_hits_bottom
    logical curve_reaches_bottom_after_cylinder
    ! if (allocated(mesh) .and. allocated(mesh%F_trm)) then
    !     allocate(prev_force(size(mesh%F_trm)))
    !     prev_force = mesh%F_trm
    ! end if
    call pg_start
    call pg_allocate_problems(1)
    call pg_bind_problem(1)

    if (iteration == 0) then
      call pg_allocate_domains(1)
      call pg_bind_domain(1)
      call pg_set_domain_equation(3)
      call pg_allocate_bounds(1)
      call pg_bind_bound(1)
      call init_geom(1)
      call pg_geom_postprocessor
      call init_gu
      !call pg_allocate_area(1)
      !call pg_bind_areapart(1)
      !call init_Mesh
      call build_mesh()
      call pg_areageom_postprocessor
    else
      call pg_allocate_domains(2)

      call pg_bind_domain(1)
      call pg_set_domain_equation(3)
      call pg_allocate_bounds(1)
      call pg_bind_bound(1)
      call extract_curve_positive_x_fast(x_, y_, n_, int(L1/50/ds_compute_curves))
      curve_hits_bottom = curve_reaches_bottom_after_cylinder(x_, y_, n_)
      if (curve_hits_bottom) then
        call init_bottom_geom_touchdown(x_, y_, n_)
      else
        call init_bottom_geom(x_,y_,n_)
      end if
      call pg_geom_postprocessor
      if (curve_hits_bottom) then
        call init_bottom_gu_touchdown()
      else
        call init_bottom_gu()
      end if
      call ga_drw_trmesh(1)

      call pg_bind_domain(2)
      call pg_set_domain_equation(21)
      call pg_allocate_bounds(1)
      call pg_bind_bound(1)
      if (curve_hits_bottom) then
        call init_top_geom_touchdown(x_, y_, n_)
      else
        call init_top_geom(x_,y_,n_)
      end if
      call pg_geom_postprocessor
      call ga_drw_trmesh(2)
      if (curve_hits_bottom) then
        call init_top_gu_touchdown()
      else
        call init_top_gu()
      end if
      
      call ga_drw_trmesh(2)
    
      !call build_mesh()
      call pg_bind_domain(2)
      call pg_bind_bound(1)
      call build_mesh_for_several_areas()
      call pg_areageom_postprocessor
      call pg_bind_domain(2)
      call pg_allocate_area_gu
      call init_Meshval()
      !call pg_bind_domain(2)
      !call pg_bind_bound(1)
    end if
    
    call ga_drw_trmesh(0)
    call pg_get_matrix
    call pg_solve
    call pg_get_psioml
    call draw_square
    gs_test_point_near_bound = .false.
    call init_zaplat
    
    if (allocated(Curves)) then
        deallocate(Curves)
    end if
    !call solve()
    call build_curve()
    call draw_Curves(2)
    call find_derivative()
    
    call find_Jacobian()
    call find_concentration_by_Jacobian()
    call build_mesh_1()
    call filling_in_the_cells()
    call draw_mesh(3)
    call draw_mesh_with_cell(1)
    
    !call build_time_isolines()
    !call solve
    call compute_force_delta(prev_force, max_delta)
    call pg_finish
    !call testing()
        
    
    if (allocated(x_)) deallocate(x_)
    if (allocated(y_)) deallocate(y_)

  end
subroutine init_Meshval()
  use mod
  real(8), allocatable:: g(:),x(:),y(:)
  integer(4) ntr, pg_get_int,i
  ntr=pg_get_int(4,1)
  allocate (g(ntr),x(ntr),y(ntr))
  call pg_get_array_real(4,1,x,ntr)
  call pg_get_array_real(4,2,y,ntr)
  do i=1,ntr
    g(i)=n0*body_force(x(i),y(i))
  enddo
  call pg_init_area_gu(g,1)
  deallocate (g,x,y)
 end

subroutine compute_force_delta(prev_force, delta)
    use mod
    real(8), allocatable :: prev_force(:)
    real(8), intent(out) :: delta
    integer(4) :: nr = 150
    integer(4) :: ng = 100
    integer(4) :: i, j, point_count
    real(8) bndg(200),bndrv(200), value_, x, y, g, r, max_value
    real(8) get_psi
    delta = 0d0
    max_value = 0d0
    point_count = (nr + 1) * (ng + 1)
    call ga_init_vneshg(ng,bndg,bndrv,H1,L1/2,8)    
    if (.not.allocated(prev_force)) then
        allocate(prev_force((nr + 1)*(ng + 1)))
        delta = 1d+2
    end if
    if (.not.allocated(mesh) .or. .not.allocated(mesh%F_trm)) then
        delta = 1d+2
        return
    end if
    do i=1,ng+1
        g=bndg(i)
        do j=1,nr+1
            r=d1+(j-d1)*(bndrv(i)-d1)/nr
            x=r*dcos(g)
            y=r*dsin(g)
            ! value_ = force_from_mesh(x,y)
            valye_ = get_psi(x,y)
            if (iteration /= 0) then
                max_value = max(max_value, abs(value_), abs(prev_force((i-1)*(nr+1)+j)))
                delta  = delta + abs(value_ - prev_force((i-1)*(nr+1)+j))
            end if
            prev_force((i-1)*(nr+1)+j) = value_
        end do
    end do
    if (iteration /= 0) then
        if (max_value > 0d0) then
            delta = delta / max_value / point_count
        else
            delta = 0d0
        end if
    end if
  end

subroutine draw_square
    use mod
    integer(4) nr,ng,i,j,mode
    real(8) g,r,x,y,psi,Vx,Vy,om,pg_get_fun_xy,bndg(200),bndrv(200), yy
    CHARACTER(LEN=30) :: filename
    ng=60 !число ячеек по gamma
    call ga_init_vneshg(ng,bndg,bndrv,H1,L1/2,8)

    WRITE(filename, '(A, I0, A)') 'data/data_fluid_', iteration, '.dat'
    OPEN (1,FILE=filename, STATUS='unknown')
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
            if (allocated(boundary_section)) then 
                call pg_bind_domain(2)
                call pg_bind_bound(1)
                if (x >= 0d0) then
                    call dcsiez(size(boundary_section), dreal(boundary_section), dimag(boundary_section), 1, x, yy)
                    if (y < yy) then
                        call pg_bind_domain(1)
                        call pg_bind_bound(1)
                    end if
                end if
            end if
            if (i==1.or.i==ng+1.or.j==1.or.j==nr+1) then
                mode=2
            else
                mode=1
            endif
            psi=pg_get_fun_xy(x,y,1,d0,d0,mode)
            !call get_uxuy(x, y, Vx, Vy)
            Vx=pg_get_fun_xy(x,y,2,d0,d1,mode)
            Vy=-pg_get_fun_xy(x,y,2,d1,d0,mode)
            om=-pg_get_fun_xy(x,y,3,d0,d0,mode)
            !if (i==1.and.j==nr+1) then
            !  Vx=pg_get_fun_xy(x-eps,y,2,d0,d1,mode)
            !  Vy=-pg_get_fun_xy(x-eps,y,2,d1,d0,mode)
            !else if (i==ng+1.and.j==nr+1) then
            !  Vx=pg_get_fun_xy(x-eps,y,2,d0,d1,mode)
            !  Vy=-pg_get_fun_xy(x-eps,y,2,d1,d0,mode)
            !else if (i==ng+1.and.j==1) then
            !  Vx=pg_get_fun_xy(x-eps,y,2,d0,d1,mode)
            !  Vy=-pg_get_fun_xy(x-eps,y,2,d1,d0,mode)
            !else if (i==1.and.j==1) then
            !  Vx=pg_get_fun_xy(x-eps,y,2,d0,d1,mode)
            !  Vy=-pg_get_fun_xy(x-eps,y,2,d1,d0,mode)
            !endif
            WRITE(1,"(F13.5, ' ', F13.5, ' ', F13.5, ' ', F13.5, ' ', F13.5, ' ', F13.5)") x,y,psi,Vx,Vy,om
        enddo
    enddo
    close(1)
    end subroutine
subroutine extract_curve_positive_x_fast(x_out, y_out, n_out, step) ! извлекает часть кривой для построения границы, где x >= 0, с шагом step
    use mod
    implicit none

    integer(4), intent(out) :: n_out
    real(8), allocatable, intent(out) :: x_out(:), y_out(:)

    integer(4) :: k
    integer(4) :: left, right, mid, start_index
    integer(4) :: n, cnt, i, step

    k = index_extreme_particles + 1
    n = Curves(k)%n

    ! Быстрые отсекающие проверки
    if (n <= 0) then
        n_out = 0
        return
    end if
    if (Curves(k)%x(n) < 0.d0) then
        n_out = 0
        return
    end if
    if (Curves(k)%x(1) >= 0.d0) then
        start_index = 1
    else
        ! Двоичный поиск первого индекса, где x >= 0
        left  = 1
        right = n
        do while (left + 1 < right)
            mid = (left + right) / 2
            if (Curves(k)%x(mid) >= 0.d0) then
                right = mid
            else
                left = mid
            end if
        end do
        start_index = right
    end if

    n_out = (n - start_index) / step + 1

    if (allocated(x_out)) deallocate(x_out)
    if (allocated(y_out)) deallocate(y_out)
    if (allocated(boundary_section)) deallocate(boundary_section)
    allocate(x_out(n_out), y_out(n_out),boundary_section(n_out))
    cnt = 0
    do i = start_index, n, step
        cnt = cnt + 1
        x_out(cnt) = Curves(k)%x(i)
        y_out(cnt) = Curves(k)%y(i)
        boundary_section(cnt) = dcmplx(x_out(cnt), y_out(cnt))
    end do
    
    x_out(1) = d0
    y_out(1) = 1d0
    boundary_section(1) = dcmplx(x_out(1), y_out(1))
    x_out(n_out) = Curves(k)%x(Curves(k)%n)
    y_out(n_out) = Curves(k)%y(Curves(k)%n)
    boundary_section(n_out) = dcmplx(x_out(n_out), y_out(n_out))
    end subroutine

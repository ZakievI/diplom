subroutine init_geom(g)
    use mod
    integer(4) g
    select case (g)
        case (1)
            call pg_allocate_boundlines(6)
            call pg_init_boundline_geomcircleds(1,ds,d0,d0,d1,pi,d0) 
            call pg_init_boundline_geomlineds(2,ds2,d1,d0,L1*d5,d0)
            call pg_init_boundline_geomlineds(3,ds2,L1*d5,d0,L1*d5,H1)
            call pg_init_boundline_geomlineds(4,ds2,L1*d5,H1,-L1*d5,H1)
            call pg_init_boundline_geomlineds(5,ds2,-L1*d5,H1,-L1*d5,d0)
            call pg_init_boundline_geomlineds(6,ds2,-L1*d5,d0,-d1,d0)
    end select
end

subroutine init_mesh
    use mod
    integer(4), parameter:: nr1 = 20 !����� ����� � ������� rr
    integer(4), parameter:: ng1 = 20 !����� ����� � ������� gg
    real(8) rr(nr1) !���������� x ����� �� �������
    real(8) gg(ng1) !���������� y ����� �� ����������
    call sred_aray(d1,h1,rr,nr1)
    call sred_aray(d0,pi,gg,ng1)
    call pg_allocate_area(1)
    call pg_bind_areapart(1)
    call ga_init_mesh_rcell_quads(rr,gg,nr1,ng1)
    call pg_areageom_postprocessor
    end subroutine

subroutine init_gu
    use mod
    call pg_allocate_bound_gu
    
    call pg_bind_boundline(1)
    call pg_init_boundline_gu_val_const(1,1,d0)
    call pg_init_boundline_gu_val_const(2,2,d0)

    call pg_bind_boundline(2)
    call pg_init_boundline_gu_val_const(1,1,d0)
    call pg_init_boundline_gu_val_const(2,3,d0)

    call pg_bind_boundline(3)
    call pg_init_boundline_gu_val_const(1,2,d0)
    call pg_init_boundline_gu_val_const(2,4,d0)

    call pg_bind_boundline(4)
    call pg_init_boundline_gu_val_const(1,1,H1)
    call pg_init_boundline_gu_val_const(2,3,d0)
    
    call pg_bind_boundline(5)
    call pg_init_boundline_gu_val_const(1,2,d0)
    call pg_init_boundline_gu_val_const(2,4,d0)
    
    call pg_bind_boundline(6)
    call pg_init_boundline_gu_val_const(1,1,d0)
    call pg_init_boundline_gu_val_const(2,3,d0)
    end

subroutine init_top_geom(x, y, n)
    use mod
    integer(4) :: n
    real(8) :: x(n), y(n)
    call pg_allocate_boundlines(6)
    call pg_init_boundline_geompolylineds(1, ds2, x, y, n-1)
    call pg_init_boundline_geomlineds(2, ds2, x(n), y(n), L1*d5, H1)
    call pg_init_boundline_geomlineds(3, ds2, L1*d5, H1, -L1*d5, H1)
    call pg_init_boundline_geomlineds(4, ds2, -L1*d5, H1, -L1*d5, d0)
    call pg_init_boundline_geomlineds(5, ds2, -L1*d5, d0, -d1, d0)
    call pg_init_boundline_geomcircleds(6, ds, d0, d0, d1, pi, pi/2) 
    end subroutine init_top_geom

subroutine init_bottom_geom(x, y, n)
    use mod
    integer(4) :: n
    real(8) :: x(n), y(n), xt(n), yt(n)
    xt = x(size(x):1:-1)
    yt = y(size(x):1:-1)
    call pg_allocate_boundlines(4)
    call pg_init_boundline_geomcircleds(1, ds, d0, d0, d1, pi/2, d0)
    call pg_init_boundline_geomlineds(2, ds2, d1, d0, L1/2, d0)
    call pg_init_boundline_geomlineds(3, ds2, L1/2, d0, xt(1), yt(1))
    call pg_init_boundline_geompolylineds(4, ds2, xt, yt, n-1)
    end subroutine init_bottom_geom

subroutine init_top_gu()
    use mod
    real(8), allocatable :: f_tau(:)
    call pg_allocate_bound_gu
    
    call pg_bind_boundline(1)
    call pg_init_boundline_gu_gen(1,1,.false.,1,1,4,1,d0,[.true.],[1],[d1])
    call pg_init_boundline_gu_gen(2,3,.false.,1,1,4,1,d0,[.true.],[3],[d1])

    call pg_bind_boundline(2)
    call pg_init_boundline_gu_val_const(1,2,d0)
    call pg_init_boundline_gu_val_const(2,4,d0)

    call pg_bind_boundline(3)
    call pg_init_boundline_gu_val_const(1,1,H1)
    call pg_init_boundline_gu_val_const(2,3,d0)

    call pg_bind_boundline(4)
    call pg_init_boundline_gu_val_const(1,2,d0)
    call pg_init_boundline_gu_val_const(2,4,d0)
    
    call pg_bind_boundline(5)
    call pg_init_boundline_gu_val_const(1,1,d0)
    call pg_init_boundline_gu_val_const(2,3,d0)
    
    call pg_bind_boundline(6)
    call pg_init_boundline_gu_val_const(1,1,d0)
    call pg_init_boundline_gu_val_const(2,2,d0)

    call pg_bind_domain(1)
    call pg_bind_bound(1)
    call pg_bind_boundline(4)
    call pg_init_boundline_gu_gen(1,2,.false.,2,1,1,1,d0,[.true.],[2],[-d1])

    call pg_init_boundline_gu_gen(2,4,.false.,2,1,1,1,d0,[.true.],[4],[-d1])
    call compute_f_tau(f_tau)
    call pg_init_boundline_gu_gen_var(2, 0, f_tau) ! f_tau со знаком "-"
    deallocate(f_tau)
    end subroutine init_top_gu

subroutine init_bottom_gu()
    use mod
    call pg_allocate_bound_gu
    
    !call pg_allocate_constvalind(2)
    !call pg_set_constvala(1,1,1) 
    !call pg_set_constvala(2,2,2)

    call pg_bind_boundline(1)
    call pg_init_boundline_gu_val_const(1,1,d0)
    call pg_init_boundline_gu_val_const(2,2,d0)

    call pg_bind_boundline(2)
    call pg_init_boundline_gu_val_const(1,1,d0)
    call pg_init_boundline_gu_val_const(2,3,d0)

    call pg_bind_boundline(3)
    call pg_init_boundline_gu_val_const(1,2,d0)
    call pg_init_boundline_gu_val_const(2,4,d0)

    end subroutine init_bottom_gu

subroutine compute_f_tau(f_tau)
    use mod 
    real(8), allocatable, intent(out) :: f_tau(:)
    complex(8), allocatable :: f_c(:), u__(:)
    real(8), allocatable :: x_c(:), y_c(:), consetr_c(:), v_x_c(:), v_y_c(:), u_x(:), u_y(:)
    complex(8) :: v__, compute_Stocks_F
    integer(4) :: n_panel, i
    external compute_Stocks_F 
    call pg_bind_domain(1)
    call pg_bind_bound(1)
    n_panel = gsbnd%line(4)%npanel
    allocate(f_tau(n_panel), u_x(n_panel), u_y(n_panel), u__(n_panel))
    allocate(x_c(n_panel), y_c(n_panel), f_c(n_panel), consetr_c(n_panel), v_x_c(n_panel), v_y_c(n_panel))
    x_c = gsbnd%xc(gsbnd%line(4)%i_begin:gsbnd%line(4)%i_end)
    y_c = gsbnd%yc(gsbnd%line(4)%i_begin:gsbnd%line(4)%i_end)

    call dcsiez(Curves(index_extreme_particles + 1)%n, Curves(index_extreme_particles + 1)%x,Curves(index_extreme_particles + 1)%c,&
             n_panel, x_c, consetr_c)
    call dcsiez(Curves(index_extreme_particles + 1)%n, Curves(index_extreme_particles + 1)%x,Curves(index_extreme_particles + 1)%V_x,&
             n_panel, x_c, v_x_c)
    call dcsiez(Curves(index_extreme_particles + 1)%n, Curves(index_extreme_particles + 1)%x,Curves(index_extreme_particles + 1)%V_y,&
             n_panel, x_c, v_y_c)
    call dcsiez(Curves(index_extreme_particles + 1)%n, Curves(index_extreme_particles + 1)%x,dreal(Curves(index_extreme_particles + 1)%u_m),&
             n_panel, x_c, u_x)
    call dcsiez(Curves(index_extreme_particles + 1)%n, Curves(index_extreme_particles + 1)%x,dimag(Curves(index_extreme_particles + 1)%u_m),&
             n_panel, x_c, u_y)
    deallocate(x_c, y_c)
    
    do i = 1, n_panel
        v__ = cmplx(v_x_c(i), v_y_c(i))
        u__ = cmplx(u_x(i), u_y(i))
        f_c(i) = compute_Stocks_F(consetr_c(i), u__(i), v__, dlt)
        f_tau(i) = (dreal(f_c(i)) * (gsbnd%x(gsbnd%line(4)%i_begin + i) - gsbnd%x(gsbnd%line(4)%i_begin + i - 1) + dimag(f_c(i)) * (gsbnd%y(gsbnd%line(4)%i_begin + i) - gsbnd%y(gsbnd%line(4)%i_begin + i - 1))))/&
        dsqrt((gsbnd%x(gsbnd%line(4)%i_begin + i) - gsbnd%x(gsbnd%line(4)%i_begin + i - 1))**2 + (gsbnd%y(gsbnd%line(4)%i_begin + i) - gsbnd%y(gsbnd%line(4)%i_begin + i - 1))**2)
    end do
    deallocate( f_c, consetr_c, v_x_c, v_y_c)
    end subroutine compute_f_tau 
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
subroutine init_geom_new(g, H1, L1)
    use mod
    real(8) :: H1,L1
    integer g
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
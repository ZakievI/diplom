subroutine init_geom_new(g, H1, L1)
    use mod
    real(8) :: H1,L1
    integer g
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
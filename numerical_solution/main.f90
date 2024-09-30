program main
  use mod
  !real(8):: H1=2*d1, L1=6*d1
  !nj=100
  !ds=pi/nj
  !ds2=ds
  !call pg_start
  !call pg_allocate_problems(1)
  !call pg_bind_problem(1)
  !call pg_allocate_domains(1)
  !call pg_bind_domain(1)
  !call pg_set_domain_equation(3)
  !call pg_allocate_bounds(1)
  !call pg_bind_bound(1)
  !call init_geom_new(1, H1, L1)
  !call pg_geom_postprocessor
  !call build_mesh_2(H1, L1)
  !call pg_areageom_postprocessor
  !call ga_drw_trmesh(1)
  !call pg_finish
  call testing()
  pause
end
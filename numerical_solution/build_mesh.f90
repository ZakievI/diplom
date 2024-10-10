subroutine build_mesh()
    use mod
    real(8) bndg(200),bndrv(200),rr(200)
    integer(4) ng,nr
    ng=60
    nr=20
    bndg=d0
    bndrv=d0
    rr=d0
    call pg_allocate_area(1)
    call pg_bind_areapart(1)
    call ga_init_vneshg(ng,bndg,bndrv,H1,L1*d5,2)
    call sred_aray(d1,H1,rr,nr+1)
    call ga_init_mesh_rcell_quads(rr,bndg,nr+1,(ng+1))
    call ga_mesh_square(d1,H1,ng,bndrv,bndg)
    !gsareapart%zm=gsareapart%zm+H1
end
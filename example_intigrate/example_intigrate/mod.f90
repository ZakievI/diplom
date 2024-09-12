module mod
use pgmod
use ifdxdy_mod
type(areapart), target :: ap
real(8), allocatable :: ff(:),err(:)
logical, allocatable :: ffknow(:)
integer(2), allocatable :: cell_domain(:)
integer(4) nf
real(8) maxErr,ds
logical emptyBoundEq
logical use_subdomain
end module
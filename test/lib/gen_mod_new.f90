module gen_mod
!include 'link_fnl_static.h'
include 'link_fnl_shared.h'
!use LSLRG_INT

IMPLICIT NONE

real(8),PARAMETER :: d0=0.0d0
real(8),PARAMETER :: d1=1.0d0
real(8),PARAMETER :: d2=2.0d0
real(8),PARAMETER :: d5=0.50d0
real(8),PARAMETER :: pi=3.14159265358979323d0
real(8),PARAMETER :: pi2=pi*d2
real(8),PARAMETER :: pi5=pi*d5
complex(8),PARAMETER :: ii=(0.0d0,1.0d0)
complex(8),PARAMETER :: c0=(0.0d0,0.0d0)
complex(8),PARAMETER :: c1=(1.0d0,0.0d0)

real(8) dbsi0,dbsi1,dbsk0,dbsk1,dcsval,dcsitg,dcsder,dbsy1,dbsy0,zarg
external dbsi0,dbsi1,dbsk0,dbsk1,dcsval,dcsitg,dcsder,dbsy1,dbsy0,zarg



end module gen_mod

module mod
use pgmod
real(8) h
real(8) k
real(8),PARAMETER :: eps=1d-4
real(8),PARAMETER :: eps_zd3=1d-3
real(8) cof_a
real(8) R_2,R_1
logical newgu
real(8) F_m
integer(4) nj ! ����������� ��������� ������� �� ������ 
real(8) ds !������ ������ ��� ���������� ������� 
real(8) ds2 !������ ������ ��� �������� �������
!type(areapart), target :: ap
real(8), allocatable :: ff(:),err(:)
logical, allocatable :: ffknow(:)

end

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
    real(8) H1, L1
    integer(4) nj ! ����������� ��������� ������� �� ������ 
    real(8) ds !������ ������ ��� ���������� ������� 
    real(8) ds2 !������ ������ ��� �������� �������
    real(8), allocatable :: ff(:),err(:)
    logical, allocatable :: ffknow(:)
    real(8) :: st                           = 100d0  
    integer(4), parameter :: N_arr          = 1500
    integer(4), parameter :: num_particle   = 50
    real(8) :: cord_extreme_particles
    integer(4) :: index_extreme_particles   = 0
    integer(4) :: N_part_1                  = 10
    integer(4) :: N_part_2                  
    integer(4) :: N_part_3                  = 10 ! = N_part_1 /(L1 - d1))
    integer(4) :: N_part_4                  = 10
    integer(4) :: gs_use_parallel_build_cerves = 1
    integer(4) :: gs_use_parallel_build_grafic = 0
    type :: Curve
        integer(4) :: n
        real(8), allocatable :: x(:), y(:), t(:)
        real(8), allocatable :: c(:)
        real(8), allocatable :: s(:)
    end type
    type :: Mesh_1
        integer(4) :: n_i, n_j
        real(8), allocatable :: x_y_(:,:), t(:,:), c(:,:), v_x(:,:), v_y(:,:), s(:,:)
    end type
    type(Curve), allocatable :: Curves(:)
    type(Mesh_1), allocatable :: mesh
contains
    subroutine finalize(obj)
        type(Curve), intent(inout) :: obj
        deallocate(obj%x,obj%y)
    end subroutine
end

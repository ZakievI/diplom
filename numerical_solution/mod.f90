module mod
    use pgmod
    real(8), PARAMETER :: eps=1d-4
    real(8), PARAMETER :: eps_zd3=1d-3
    !real(8) cof_a
    !real(8) R_2,R_1
    logical newgu
    logical zaplat
    real(8) ds_pg
    real(8) d_zapl, cc1_zapl, dd1_zapl, cc2_zapl, dd2_zapl
    real(8) F_m
    real(8), PARAMETER :: n0 = 1d0
    real(8) H1, L1
    integer(4) iteration
    integer(4) nj ! ����������� ��������� ������� �� ������ 
    real(8) ds !������ ������ ��� ���������� ������� 
    real(8) ds2 !������ ������ ��� �������� �������
    real(8), allocatable :: ff(:),err(:)
    logical, allocatable :: ffknow(:)
    real(8) :: st                           = 5d0
    real(8) :: mu                           = 1 
    integer(4), parameter :: N_arr          = 10000
    integer(4), parameter :: num_particle   = 30
    real(8) :: cord_extreme_particles
    integer(4) :: index_extreme_particles   = 1
    integer(4) :: N_part_1                  = 30
    integer(4) :: N_part_2                  
    integer(4) :: N_part_3                  = 70 ! = N_part_1 /(L1 - d1))
    integer(4) :: N_part_4                  = 20
    integer(4) :: use_parallel_build_cerves = 0
    integer(4) :: use_parallel_build_cerves_J = 1
    integer(4) :: gs_use_parallel_build_grafic = 0
    real(8)    :: dlt
    complex(8), allocatable :: boundary_section(:) 
    
    
    type :: Curve
        integer(4) :: n
        real(8), allocatable :: x(:), y(:), t(:), V_x(:), V_y(:)
        real(8), allocatable :: c(:)
        real(8), allocatable :: du1dx1(:), du2dx1(:), du1dx2(:), du2dx2(:)
        real(8), allocatable :: J_11(:), J_12(:), J_21(:), J_22(:), J(:)
        complex(8), allocatable :: u_m(:) !скорость жидкости, только для вычисления силы первой вылетевшей частицы
        real(8), allocatable :: s(:)
        real(8), allocatable :: discrepancy_du1dx1(:), discrepancy_du2dx1(:), discrepancy_du1dx2(:), discrepancy_du2dx2(:)
        real(8), allocatable :: aprox_du1dx1(:), aprox_du2dx1(:), aprox_du1dx2(:), aprox_du2dx2(:)
    end type
    
    real(8) :: ds_compute_curves            = 1d-2
    real(8) :: tol_compute_curves           = 1d-4
    
    type :: Mesh_1
        integer(4) :: n_i, n_j
        real(8), allocatable :: t(:), c(:), v_x(:), v_y(:), s(:)
        complex(8), allocatable :: z_m(:) !координаты узлов сетки в области
        complex(8), allocatable :: v_m(:) !скорость частицы
        complex(8), allocatable :: f_m(:) !плотность масс сил
        integer(4) :: npe = 4 !число углов в элементе (3,4) - для выделения памяти
        integer(4) n !число узлов
        integer(4) ntr !количество ячеек в сетке
        integer(4), allocatable :: trm(:,:) !(npe,ntr) !индексы вершин треугольной сетки
        real(4), allocatable :: F_trm(:) ! значение rot(f) в ячейке
    end type
    type(Curve), allocatable, target :: Curves(:)
    type(Curve), pointer :: current_Curve => null()
    !$omp threadprivate(current_Curve)
    type(Mesh_1), allocatable :: mesh
    real(8) :: top_coordinat                = 5d0
    real(8) :: bottom_coordinat             = 0d0
contains
function body_force(x,y) result(fm)
    real(8), intent(in) :: x,y
    real(8) :: fm
    if (allocated(mesh) .and. allocated(mesh%F_trm) .and. mesh%ntr > 0) then
        fm = force_from_mesh(x,y)
    else
        fm = 0d0
    end if
    end function body_force

function force_from_mesh(x,y) result(fm)
    real(8), intent(in) :: x,y
    real(8) :: fm
    complex(8) :: z1, z2, z3, z4
    integer(4) :: i
    fm = d0
    do i = 1, mesh%ntr
        z1 = mesh%z_m(mesh%trm(1, i))
        z2 = mesh%z_m(mesh%trm(2, i))
        z3 = mesh%z_m(mesh%trm(3, i))
        z4 = mesh%z_m(mesh%trm(4, i))
        if (point_in_triangle(x,y,z1,z2,z3) .or. point_in_triangle(x,y,z1,z3,z4)) then
            fm = mesh%F_trm(i)
            return
        end if
    end do
    end function force_from_mesh

logical function point_in_triangle(x,y,z1,z2,z3)
    real(8), intent(in) :: x,y
    complex(8), intent(in) :: z1, z2, z3
    real(8) :: x1, y1, x2, y2, x3, y3
    real(8) :: d1, d2, d3
    logical :: has_neg, has_pos
    x1 = dreal(z1)
    y1 = dimag(z1)
    x2 = dreal(z2)
    y2 = dimag(z2)
    x3 = dreal(z3)
    y3 = dimag(z3)
    d1 = (x - x2) * (y1 - y2) - (x1 - x2) * (y - y2)
    d2 = (x - x3) * (y2 - y3) - (x2 - x3) * (y - y3)
    d3 = (x - x1) * (y3 - y1) - (x3 - x1) * (y - y1)
    has_neg = (d1 < -eps) .or. (d2 < -eps) .or. (d3 < -eps)
    has_pos = (d1 > eps) .or. (d2 > eps) .or. (d3 > eps)
    point_in_triangle = .not. (has_neg .and. has_pos)
    end function point_in_triangle

subroutine finalize(obj)
        type(Curve), intent(inout) :: obj
        deallocate(obj%x,obj%y)
    end subroutine
end

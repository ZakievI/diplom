module func_mod
  
  
interface
  
  !!!�����, ����� ����� ���� ���������� allocatable ������ ��� ��������
  subroutine line_array_ds(x1,y1,x2,y2,ds,x,y,s,n,min_np)
  integer(4) n
  real(8) x1,y1,x2,y2 !���������� ������ �������
  real(8) ds !��� ��������������
  real(8), allocatable :: x(:),y(:),s(:)
  integer(4) min_np !����������� ����� �����
  end
  
  function itoa(i) result(res)
  character(:),allocatable :: res
  integer,intent(in) :: i
  end
  
  function ftoa(f) result(res)
  character(:),allocatable :: res
  real(8),intent(in) :: f
  end
  
  subroutine LS_approx3(n,m,x,y,cc,fcn_LS_approx)
  integer(4) n !���������� ������� �������������
  integer(4) m !����� ����� (x_j,y_j)
  real(8), target :: x(m),y(m) !���������� ���������������� �����
  !real(8) w(m)      !��� ������ ����� L=sum_(j=1)^m w(j)*(F(cc,x_j)-y_j)^2
  real(8) cc(n)   !������������ ����������������� �������
  real(8) fcn_LS_approx
  external fcn_LS_approx !���������������� ������� fcn_LS_approx(m,cc,x) x - ���������� ����� ����� x 
  end
  
  subroutine LS_approx4(n,m,x,y,w,cc,fcn_LS_approx)
  integer(4) n !���������� ������� �������������
  integer(4) m !����� ����� (x_j,y_j)
  real(8), target :: x(m),y(m) !���������� ���������������� �����
  real(8), target :: w(m)      !��� ������ ����� L=sum_(j=1)^m w(j)*(F(cc,x_j)-y_j)^2
  real(8) cc(n)   !������������ ����������������� �������
  real(8) fcn_LS_approx
  external fcn_LS_approx !���������������� ������� fcn_LS_approx(m,cc,x) x - ���������� ����� ����� x 
  end

end interface
  
  ABSTRACT INTERFACE
  FUNCTION fcn_LS_approx_template(n,cc,x) RESULT(f)
      integer(4) n
      real(8) cc(n)
      real(8) x,f
  END FUNCTION fcn_LS_approx_template
  END INTERFACE
  
  PROCEDURE(fcn_LS_approx_template), POINTER :: func_fcn_LS_approx
  real(8), pointer :: func_fcn_LS_approx_x(:),func_fcn_LS_approx_y(:),func_fcn_LS_approx_w(:)
  
endmodule

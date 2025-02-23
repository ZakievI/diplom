!для всех функций
!g - тип группы
             !1-gsbnd (граница)
			 !2-gsbndl (участок границы)
			 !3-gs%m (главная матрица)
			 !4-gsarea%a (область)
       !5-gs%constval   
       !6-gsareapart (участок области)

function pg_get_int(g,k)
!dec$ attributes dllexport:: pg_get_int
!получить целое число для текущей границы
use pgmod
integer(4) pg_get_int
integer(4) g !тип группы
integer(4) k !тип значения
             !ig==1
               !1-npanel
			 !ig==2
			   !1-npanel
			   !2-i_begin
			   !3-i_end
			 !ig==3
			   !1-nx
			   !2-nu
			 !ig==4
			   !1-ntr
			   !2-n
      !ig==5
         !1-constvaln
      !ig==6
			   !1-ntr
			   !2-n
         !3-ntr_begin
         !4-ntr_end
         !5-n_begin
         !6-n_end
select case (g)        
case (1) !g==1
  select case (k)
  case (1) 
    pg_get_int=gsbnd%npanel
  end select
case (2) !g==2
  select case (k)
  case (1) 
    pg_get_int=gsbndl%npanel
  case (2) 
    pg_get_int=gsbndl%i_begin
  case (3) 
    pg_get_int=gsbndl%i_end
  end select
case (3) !g==3
  select case (k)
  case (1) 
    pg_get_int=gs%m%nx
  case (2) 
    pg_get_int=gs%m%nu
  end select
case (4) !g==4
  select case (k)
  case (1) 
    pg_get_int=gsarea%a%ntr
  case (2) 
    pg_get_int=gsarea%a%n
  end select
case (5) !g==5
  select case (k)
  case (1) 
    pg_get_int=gs%constvaln
  end select
case (6) !g==6
  select case (k)
  case (1) 
    pg_get_int=gsareapart%ntr
  case (2) 
    pg_get_int=gsareapart%n
  case (3) 
    pg_get_int=gsareapart%ntr_begin
  case (4) 
    pg_get_int=gsareapart%ntr_end
  case (5) 
    pg_get_int=gsareapart%n_begin
  case (6) 
    pg_get_int=gsareapart%n_end
  end select
end select
end

subroutine pg_get_array_real(g, k, arr,n)
!dec$ attributes dllexport:: pg_get_array_real
!получить массив вещественных чисел
use pgmod
integer(4) n,i,j
integer(4) g !тип группы
integer(4) k !тип значения
             !g==1
               !1-x(npanel+1)
			         !2-y(npanel+1)
			         !3-xc(npanel)
			         !4-yc(npanel)
			         !5-s(npanel+1)
			         !6-sc(npanel)
			         !7-l(npanel)
			         !8-tt(npanel) 
			       !g==2
               !1-x(npanel+1)
			         !2-y(npanel+1)
			         !3-xc(npanel)
			         !4-yc(npanel)
			         !5-s(npanel+1)
			         !6-sc(npanel)
			         !7-l(npanel)
			         !8-tt(npanel) 
             !g==4
			         !1-xmc(ntr)
			         !2-ymc(ntr)
			         !3-x(n)
			         !4-y(n)
             !g=5
               !1-constvals(constvaln)
             !g==6
			         !1-xmc(ntr)
			         !2-ymc(ntr)
			         !3-x(n)
			         !4-y(n)
real(8) arr(n)
select case (g)        
case (1) !g==1
  select case (k)
  case (1) 
    arr(1:gsbnd%npanel+1)=gsbnd%x(1:gsbnd%npanel+1)
  case (2) 
    arr(1:gsbnd%npanel+1)=gsbnd%y(1:gsbnd%npanel+1)
  case (3) 
    arr(1:gsbnd%npanel)=gsbnd%xc(1:gsbnd%npanel)
  case (4) 
    arr(1:gsbnd%npanel)=gsbnd%yc(1:gsbnd%npanel)
  case (5) 
    arr(1:gsbnd%npanel+1)=gsbnd%s(1:gsbnd%npanel+1)
  case (6) 
    arr(1:gsbnd%npanel)=gsbnd%sc(1:gsbnd%npanel)
  case (7) 
    arr(1:gsbnd%npanel)=gsbnd%l(1:gsbnd%npanel)
  case (8) 
    do i=1,gsbnd%npanel
	    arr(i)=zarg(gsbnd%ett(i))
	  enddo
  end select
case (2) !g==2
  select case (k)
  case (1) 
    arr(1:gsbndl%npanel+1)=gsbndl%x
  case (2) 
    arr(1:gsbndl%npanel+1)=gsbndl%y
  case (3) 
    arr(1:gsbndl%npanel)=gsbnd%xc(gsbndl%i_begin:gsbndl%i_end)
  case (4) 
    arr(1:gsbndl%npanel)=gsbnd%yc(gsbndl%i_begin:gsbndl%i_end)
  case (5) 
    arr(1:gsbndl%npanel+1)=gsbnd%s(gsbndl%i_begin:gsbndl%i_end+1)
  case (6) 
    arr(1:gsbndl%npanel)=gsbnd%sc(gsbndl%i_begin:gsbndl%i_end)
  case (7) 
    arr(1:gsbndl%npanel)=gsbnd%l(gsbndl%i_begin:gsbndl%i_end)
  case (8) 
    j=gsbndl%i_begin
    do i=1,gsbndl%npanel
	    arr(i)=zarg(gsbnd%ett(j))
	    j=j+1
	  enddo
  end select
case (4) !g==4
  select case (k)
  case (1) 
    arr(1:gsarea%a%ntr)=dreal(gsarea%a%zmc(:))
  case (2) 
    arr(1:gsarea%a%ntr)=dimag(gsarea%a%zmc(:))
  case (3) 
    arr(1:gsarea%a%n)=dreal(gsarea%a%zm(:))
  case (4) 
    arr(1:gsarea%a%n)=dimag(gsarea%a%zm(:))
  end select
case (5) !g==5
  select case (k)
  case (1)
    arr(1:gs%constvaln)=gs%constvals
  end select
case (6) !g==6
  select case (k)
  case (1) 
    arr(1:gsareapart%ntr)=dreal(gsarea%a%zmc(gsareapart%ntr_begin:gsareapart%ntr_end))
  case (2) 
    arr(1:gsareapart%ntr)=dimag(gsarea%a%zmc(gsareapart%ntr_begin:gsareapart%ntr_end))
  case (3) 
    if (gsarea%a%geom_inited) then
      arr(1:gsareapart%n)=dreal(gsarea%a%zm(gsareapart%n_begin:gsareapart%n_end))
    else
      arr(1:gsareapart%n)=dreal(gsareapart%zm(:))
    endif
  case (4) 
    if (gsarea%a%geom_inited) then
      arr(1:gsareapart%n)=dimag(gsarea%a%zm(gsareapart%n_begin:gsareapart%n_end))
    else
      arr(1:gsareapart%n)=dimag(gsareapart%zm(:))
    endif
  end select
end select
end

subroutine pg_get_array_real2(g, k, j, arr, n)
!dec$ attributes dllexport:: pg_get_array_real2
!получить массив вещественных чисел из двумерного массива
use pgmod
integer(4) n
integer(4) g !тип группы
integer(4) j !второй индекс 
integer(4) k !тип значения
             !g==1
               !1-psiom(npanel)  j=1- psi, 2- dpsi/dn, 3- om, 4- dom/dn
             !g==2
               !1-psiom(npanel)  j=1- psi, 2- dpsi/dn, 3- om, 4- dom/dn
real(8) arr(n)
select case (g)        
case (1) !g==1
  select case (k)
  case (1) 
    arr(1:gsbnd%npanel)=gsbnd%psiom(1:gsbnd%npanel,j)
  end select
case (2) !g==2
  select case (k)
  case (1) 
    arr(1:gsbndl%npanel)=gsbnd%psiom(gsbndl%i_begin:gsbndl%i_end,j)
  end select
end select
end

subroutine pg_get_array_int(g, k, arr, n)
!dec$ attributes dllexport:: pg_get_array_int
!получить массив целых чисел из двумерного массива
use pgmod
integer(4) n
integer(4) g !тип группы
integer(4) k !тип значения
			       !g=5
               !1-constvalinds(constvaln)
integer(4) arr(n)
select case (g)        
case (5) !g==5
  select case (k)
  case (1) 
    arr(1:gs%constvaln)=gs%constvalinds
  end select
end select
end


subroutine pg_get_array_int2(g, k, j, arr, n)
!dec$ attributes dllexport:: pg_get_array_int2
!получить массив целых чисел из одномерного массива
use pgmod
integer(4) n
integer(4) g !тип группы
integer(4) j !второй индекс
integer(4) k !тип значения
            !g==1
              !1-psiind(npanel) j=1- psi, 2- dpsi/dn, 3- om, 4- dom/dn
            !g==4
              !1-trm(ntr) j-номер вершины (1..3,<4>)
            !g==6
              !1-trm(ntr) j-номер вершины (1..3,<4>)
integer(4) arr(n)
select case (g)        
case (1) !g==1
  select case (k)
  case (1) 
    arr(1:gsbnd%npanel)=gsbnd%psiind(:,j)
  end select
case (4) !g==4
  select case (k)
  case (1) 
    arr(1:gsarea%a%ntr)=gsarea%a%trm(j,:)
  end select
case (6) !g==6
  select case (k)
  case (1) 
    if (gsarea%a%geom_inited) then
      arr(1:gsareapart%ntr)=gsarea%a%trm(j,gsareapart%ntr_begin:gsareapart%ntr_end)
    else
      arr(1:gsareapart%ntr)=gsareapart%trm(j,:)
    endif
  end select
end select
end

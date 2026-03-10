function get_funa(r,tet)
!аналитическое решение уравнения лапласа в области
!потенциальное течение
use pgmod2
real(8) r,tet,ka,kb
real(8) get_funa
kb=h**2/(h**2-d1)
ka=-kb
get_funa=(ka/r+kb*r)*dsin(tet) !-12.5d0
end

function get_funa_oss(r,tet,vx,vy)
!аналитическое решение уравнения на функцию тока при обтекании сферы
!осесимметричное течение
use pgmod2
real(8) r,tet,vx,vy,r3
real(8) get_funa_oss
r3=d1/r**3
get_funa_oss=d5*(r*dsin(tet))**2*(d1-d1*r3)
vx=d1+(d1-3.0d0*dcos(tet)**2)*d5*r3
vy=-0.75d0*dsin(d2*tet)*r3
end

function get_func(zz)
!аналитическое решение обтекания пористого цилиндра находящегося в пористой среде
use pgmod2
complex(8) zz,ff
real(8) get_func,k
if (cdabs(zz)>d1) then
  k=(d1-kkr)/(d1+kkr)
  ff=zz+k/zz
else
  k=d2*kkr/(d1+kkr)
  ff=k*zz
endif
get_func=dimag(ff)
end

function get_funb(xx,yy)
!аналитическое решение уравнения лапласа, Пуассона или Гельмгольца в квадрате
use pgmod2
real(8) xx,yy
real(8) get_funb

!Лаплас
if (nfun==1) get_funb=xx
if (nfun==2) get_funb=d1-xx-yy+xx*yy

!Пуассон
if (nfun==3) get_funb=xx*xx
if (nfun==4) get_funb=dsin(xx)

!Гельмгольц
if (nfun==5) get_funb=dexp(xx) 
if (nfun==6) get_funb=dsin(k_helm*xx)
if (nfun==7) get_funb=yy*dsin(xx)
if (nfun==8) get_funb=yy*dsin(k_helm*xx)
if (nfun==9) get_funb=yy*dexp(k_helm*xx) 
if (nfun==10) get_funb=dsin(xx*yy) 
if (nfun==17) get_funb=(xx+d1)**3

!Гельмгольц неоднородный
if (nfun==11) get_funb=(xx*yy)**2
if (nfun==12) get_funb=dsin(xx*yy) 
if (nfun==13) get_funb=xx**2*dsin(yy)
if (nfun==14.or.nfun==15.or.nfun==16) get_funb=xx*yy*(d1-dexp((xx-d1)/e_helm))*(d1-dexp((yy-d1)/e_helm)) 
  !(часть коэффициентов при неизвестных в области = 0)
if (nfun==46) get_funb=dsin(xx)

!уравнение переноса (однородное)
if (nfun==26) get_funb=yy*yy+d2*xx
if (nfun==27) get_funb=yy*dexp(xx)

!уравнение переноса (неоднородное)
if (nfun==28) get_funb=xx*xx+yy*yy
if (nfun==29) get_funb=xx+yy
if (nfun==37.or.nfun==47) get_funb=xx*xx+yy*yy

!уравнение Лапласа осесимметричное
if (nfun==38) get_funb=xx

!уравнение осесимметричное для функции тока
if (nfun==39) get_funb=yy**2*d5

!уравнение Пуассона с неизвестными в правой части
if (nfun==49) get_funb=yy*dsin(xx)
end

function get_funb2(xx,yy)
!лапласиан аналитического решения уравнения лапласа, Пуассона или Гельмгольца в квадрате
use pgmod2
real(8) get_funb2,xx,yy
get_funb2=yy*(d1-dexp((yy-d1)/e_helm))*(-dexp((xx-d1)/e_helm)/e_helm)*(d2+xx/e_helm)+xx*(d1-dexp((xx-d1)/e_helm))*(-dexp((yy-d1)/e_helm)/e_helm)*(d2+yy/e_helm)
end

function get_funb1(xx,yy,ku)
!производная по нормали на границе от аналитического решения уравнения лапласа, Пуассона или Гельмгольца в квадрате
use pgmod2
integer(4) ku
real(8) xx,yy
real(8) get_funb1

!Лаплас

if (nfun==1) then
  !get_funb=xx
  if (ku==1) then
    get_funb1=d0
  elseif (ku==2) then
    get_funb1=d1
  elseif (ku==3) then
    get_funb1=d0
  elseif (ku==4) then
    get_funb1=-d1
  endif
endif

if (nfun==2) then
  !get_funb=d1-xx-yy+xx*yy
  if (ku==1) then
    get_funb1=d1-xx
  elseif (ku==2) then
    get_funb1=-d1+yy
  elseif (ku==3) then
    get_funb1=-d1+xx
  elseif (ku==4) then
    get_funb1=d1-yy
  endif
endif

!Пуассон

if (nfun==3) then
  !get_funb=xx*xx
  if (ku==1) then
    get_funb1=d0
  elseif (ku==2) then
    get_funb1=d2*xx
  elseif (ku==3) then
    get_funb1=d0
  elseif (ku==4) then
    get_funb1=-d2*xx
  endif
endif

if (nfun==4) then
  !get_funb=dsin(xx)
  if (ku==1) then
    get_funb1=d0
  elseif (ku==2) then
    get_funb1=dcos(xx)
  elseif (ku==3) then
    get_funb1=d0
  elseif (ku==4) then
    get_funb1=-dcos(xx)
  endif
endif

!Гельмгольц

if (nfun==5) then
  !get_funb=dexp(xx)
  if (ku==1) then
    get_funb1=d0
  elseif (ku==2) then
    get_funb1=dexp(xx)
  elseif (ku==3) then
    get_funb1=d0
  elseif (ku==4) then
    get_funb1=-dexp(xx)
  endif
endif

if (nfun==6) then
  !get_funb=dsin(k_helm*xx)
  if (ku==1) then
    get_funb1=d0
  elseif (ku==2) then
    get_funb1=k_helm*dcos(k_helm*xx)
  elseif (ku==3) then
    get_funb1=d0
  elseif (ku==4) then
    get_funb1=-k_helm*dcos(k_helm*xx)
  endif
endif

if (nfun==7) then
  !get_funb=yy*dsin(xx)
  if (ku==1) then
    get_funb1=-dsin(xx)
  elseif (ku==2) then
    get_funb1=yy*dcos(xx)
  elseif (ku==3) then
    get_funb1=dsin(xx)
  elseif (ku==4) then
    get_funb1=-yy*dcos(xx)
  endif
endif

if (nfun==8) then
  !get_funb=yy*dsin(k_helm*xx)
  if (ku==1) then
    get_funb1=-dsin(k_helm*xx)
  elseif (ku==2) then
    get_funb1=k_helm*yy*dcos(k_helm*xx)
  elseif (ku==3) then
    get_funb1=dsin(k_helm*xx)
  elseif (ku==4) then
    get_funb1=-k_helm*yy*dcos(k_helm*xx)
  endif
endif

if (nfun==9) then
  !get_funb=yy*dexp(k_helm*xx)
  if (ku==1) then
    get_funb1=-dexp(k_helm*xx)
  elseif (ku==2) then
    get_funb1=k_helm*yy*dexp(k_helm*xx)
  elseif (ku==3) then
    get_funb1=dexp(k_helm*xx)
  elseif (ku==4) then
    get_funb1=-k_helm*yy*dexp(k_helm*xx)
  endif
endif

if (nfun==10) then
  !get_funb=dsin(xx*yy)
  if (ku==1) then
    get_funb1=xx*dcos(xx*yy)
  elseif (ku==2) then
    get_funb1=yy*dcos(xx*yy)
  elseif (ku==3) then
    get_funb1=xx*dcos(xx*yy)
  elseif (ku==4) then
    get_funb1=-yy*dcos(xx*yy)
  endif
endif

if (nfun==17) then
  !get_funb=(xx+d1)**3
  if (ku==1) then
    get_funb1=d0
  elseif (ku==2) then
    get_funb1=3.0d0*(xx+d1)**2
  elseif (ku==3) then
    get_funb1=d0
  elseif (ku==4) then
    get_funb1=-3.0d0*(xx+d1)**2
  endif
endif

!Гельмгольц неоднородный

if (nfun==11) then
  !get_funb=(xx*yy)**2
  if (ku==1) then
    get_funb1=-d2*yy*xx**2
  elseif (ku==2) then
    get_funb1=d2*xx*yy**2
  elseif (ku==3) then
    get_funb1=d2*yy*xx**2
  elseif (ku==4) then
    get_funb1=-d2*xx*yy**2
  endif
endif

if (nfun==12) then
  !get_funb=dsin(xx*yy)
  if (ku==1) then
    get_funb1=-xx*dcos(xx*yy)
  elseif (ku==2) then
    get_funb1=yy*dcos(xx*yy)
  elseif (ku==3) then
    get_funb1=xx*dcos(xx*yy)
  elseif (ku==4) then
    get_funb1=-yy*dcos(xx*yy)
  endif
endif

if (nfun==13) then
  !get_funb=xx**2*dsin(yy)
  if (ku==1) then
    get_funb1=-xx**2*dcos(yy)
  elseif (ku==2) then
    get_funb1=d2*xx*dsin(yy)
  elseif (ku==3) then
    get_funb1=xx**2*dcos(yy)
  elseif (ku==4) then
    get_funb1=-d2*xx*dsin(yy)
  endif
endif

if (nfun==14.or.nfun==15.or.nfun==16) then
  !get_funb=xx*yy*(d1-dexp((xx-d1)/e_helm))*(d1-dexp((yy-d1)/e_helm))
  if (ku==1) then
    get_funb1=-xx*(d1-dexp((xx-d1)/e_helm))*((d1-dexp((yy-d1)/e_helm))+yy/e_helm*(-dexp((yy-d1)/e_helm)))
  elseif (ku==2) then
    get_funb1=yy*(d1-dexp((yy-d1)/e_helm))*((d1-dexp((xx-d1)/e_helm))+xx/e_helm*(-dexp((xx-d1)/e_helm)))
  elseif (ku==3) then
    get_funb1=xx*(d1-dexp((xx-d1)/e_helm))*((d1-dexp((yy-d1)/e_helm))+yy/e_helm*(-dexp((yy-d1)/e_helm)))
  elseif (ku==4) then
    get_funb1=-yy*(d1-dexp((yy-d1)/e_helm))*((d1-dexp((xx-d1)/e_helm))+xx/e_helm*(-dexp((xx-d1)/e_helm)))
  endif
endif

  !(часть коэффициентов при неизвестных в области = 0)
if (nfun==46) then
  !get_funb=dsin(xx)
  if (ku==1) then
    get_funb1=d0
  elseif (ku==2) then
    get_funb1=dcos(xx)
  elseif (ku==3) then
    get_funb1=d0
  elseif (ku==4) then
    get_funb1=-dcos(xx)
  endif
endif

!уравнение переноса (однородное)
if (nfun==26) then
  !get_funb=yy*yy+d2*xx
  if (ku==1) then
    get_funb1=-d2*yy
  elseif (ku==2) then
    get_funb1=d2
  elseif (ku==3) then
    get_funb1=d2*yy
  elseif (ku==4) then
    get_funb1=-d2
  endif
endif

if (nfun==27) then
  !get_funb=yy*dexp(xx)
  if (ku==1) then
    get_funb1=-dexp(xx)
  elseif (ku==2) then
    get_funb1=yy*dexp(xx)
  elseif (ku==3) then
    get_funb1=dexp(xx)
  elseif (ku==4) then
    get_funb1=-yy*dexp(xx)
  endif
endif

!уравнение переноса (неоднородное)
if (nfun==28) then
  !get_funb=xx*xx+yy*yy
  if (ku==1) then
    get_funb1=-d2*yy
  elseif (ku==2) then
    get_funb1=d2*xx
  elseif (ku==3) then
    get_funb1=d2*yy
  elseif (ku==4) then
    get_funb1=-d2*xx
  endif
endif

if (nfun==29) then
  !get_funb=xx+yy
  if (ku==1) then
    get_funb1=-d1
  elseif (ku==2) then
    get_funb1=d1
  elseif (ku==3) then
    get_funb1=d1
  elseif (ku==4) then
    get_funb1=-d1
  endif
endif

if (nfun==37.or.nfun==47) then
  !get_funb=xx*xx+yy*yy
  if (ku==1) then
    get_funb1=-d2*yy
  elseif (ku==2) then
    get_funb1=d2*xx
  elseif (ku==3) then
    get_funb1=d2*yy
  elseif (ku==4) then
    get_funb1=-d2*xx
  endif
endif

!уравнение Лапласа осесимметричное
if (nfun==38) then
  !get_funb=xx
  if (ku==1) then
    get_funb1=d0
  elseif (ku==2) then
    get_funb1=d1
  elseif (ku==3) then
    get_funb1=d0
  elseif (ku==4) then
    get_funb1=-d1
  endif
endif

!уравнение осесимметричное для функции тока
if (nfun==39) then
  !get_funb=yy**2*d5
  if (ku==1) then
    get_funb1=-yy
  elseif (ku==2) then
    get_funb1=d0
  elseif (ku==3) then
    get_funb1=yy
  elseif (ku==4) then
    get_funb1=d0
  endif
endif

end

function get_funbxy(xx,yy,ku)
!производная по нормали на границе от аналитического решения уравнения лапласа, Пуассона или Гельмгольца в квадрате
use pgmod2
integer(4) ku !1-d/dx, 2-d/dy
real(8) xx,yy
real(8) get_funbxy

get_funbxy=d0

!Гельмгольц
if (nfun==5) then
  !get_funb=dexp(xx)
  if (ku==1) then
    get_funbxy=dexp(xx)
  elseif (ku==2) then
    get_funbxy=d0
  endif
endif

if (nfun==9) then
  !get_funb=yy*dexp(k_helm*xx) 
  if (ku==1) then
    get_funbxy=k_helm*yy*dexp(k_helm*xx) 
  elseif (ku==2) then
    get_funbxy=dexp(k_helm*xx) 
  endif
endif
end

function get_fund(xx,yy)
!аналитическое решение бигармонического уравнения
use pgmod2
real(8) xx,yy
real(8) get_fund

!однородное бигармоническое
if (nfun==18) get_fund=xx
if (nfun==19) get_fund=xx**3
if (nfun==20) get_fund=yy*dsin(k_helm*xx)
if (nfun==21) get_fund=yy*dexp(k_helm*xx)
if (nfun==23) get_fund=xx**4

!неоднородное бигармоническое
if (nfun==22) get_fund=xx**4

!неоднородная система
if (nfun==24) get_fund=dsin(xx*yy)

!уравнение Бринкмана через единую функцию Грина
if (nfun==25) get_fund=yy*dexp(k_helm*xx)

!неоднородное бигармоническое уравнение с неизвестными в правой части
if (nfun==42) get_fund=xx**4*yy
if (nfun==43.or.nfun==44) get_fund=xx**2+yy**2
if (nfun==45.or.nfun==48) get_fund=yy*sin(xx)


!get_fund=-dsin(xx)
!get_fund=dsin(k_helm*xx)
end

function get_fund2(xx,yy)
!лапласиан аналитического решения бигармонического уравнения
use pgmod2
real(8) xx,yy
real(8) get_fund2

!однородное бигармоническое
if (nfun==18) get_fund2=d0
if (nfun==19) get_fund2=6.0d0*xx
if (nfun==20) get_fund2=yy*dsin(k_helm*xx)
if (nfun==21) get_fund2=yy*dexp(k_helm*xx)
if (nfun==23) get_fund2=6.0d0*xx

!неоднородное бигармоническое
if (nfun==22) get_fund2=12.0d0*xx**2

!неоднородная система
if (nfun==24) get_fund2=dsin(xx*yy)

!уравнение Бринкмана через единую функцию Грина
if (nfun==25) get_fund2=yy*dexp(k_helm*xx)*k_helm**2

!неоднородное бигармоническое уравнение с неизвестными в правой части
if (nfun==42) get_fund2=12.0d0*xx**2*yy
if (nfun==43.or.nfun==44) get_fund2=4.0d0
if (nfun==45.or.nfun==48) get_fund2=-yy*sin(xx)

!get_fund2=dsin(xx)
!get_fund2=dsin(k_helm*xx)
end

function fz_puasson(z)
real(8) fz_puasson,f_puasson
complex(8) z
fz_puasson=f_puasson(dreal(z),dimag(z))
end

function f_puasson(xx,yy)
use pgmod2
real(8) f_puasson,xx,yy
real(8) get_funb2,get_funb,f_helmgolz

f_puasson=d0

!Пуассон
if (nfun==3) f_puasson=d2
if (nfun==4) f_puasson=-dsin(xx)

!Гельмгольц неоднородный
if (nfun==11) f_puasson=d2*(xx**2+yy**2+(xx*yy)**2)
if (nfun==12) f_puasson=(d1-xx**2-yy**2)*dsin(xx*yy)
if (nfun==13) f_puasson=(d2+xx**2*(yy-d1))*dsin(yy)
if (nfun==14.or.nfun==15.or.nfun==16) f_puasson=get_funb2(xx,yy)+f_helmgolz(xx,yy)*get_funb(xx,yy) 
if (nfun==46) then
  if (xx<d5) then
    f_puasson=-dsin(xx)
  else
    f_puasson=d0
  endif
endif

!неоднородное бигармоническое
if (nfun==22) f_puasson=24.0d0

!неоднородное уравнение переноса
if (nfun==28) f_puasson=d2*(d2-xx**2-yy**3)
if (nfun==29) f_puasson=-xx-yy
if (nfun==37) f_puasson=4.0d0-d2*xx*yy**2
if (nfun==47) then
  if (xx<d5) then 
    f_puasson=4.0d0-d2*yy**2
  else
    f_puasson=4.0d0-d2*xx**2
  endif
endif
!f_puasson=d0
end

function fz_helmgolz(z)
real(8) fz_helmgolz,f_helmgolz
complex(8) z
fz_helmgolz=f_helmgolz(dreal(z),dimag(z))
end

function f_helmgolz(xx,yy)
!множитель перед f в уравнении гельмгольца
use pgmod2
real(8) f_helmgolz,xx,yy
if (nfun==4) f_helmgolz=d1
if (nfun==5) f_helmgolz=-d1
if (nfun==6) f_helmgolz=k_helm**2 
if (nfun==7) f_helmgolz=d1
if (nfun==8) f_helmgolz=k_helm**2 
if (nfun==9) f_helmgolz=-k_helm**2 
if (nfun==10) f_helmgolz=xx**2+yy**2
if (nfun==11) f_helmgolz=d2
if (nfun==12) f_helmgolz=d1
if (nfun==13) f_helmgolz=yy
if (nfun==14) f_helmgolz=d1/e_helm
if (nfun==15) f_helmgolz=(xx+yy)/e_helm
if (nfun==16) f_helmgolz=dcos(xx*yy)/e_helm
if (nfun==17) f_helmgolz=-6.0d0/(xx+d1)**2
if (nfun==24) f_helmgolz=xx**2+yy**2
if (nfun==25) f_helmgolz=-k_helm**2 
if (nfun==46) then
  if (xx<d5) then
    f_helmgolz=d0
  else
    f_helmgolz=d1
  endif
endif
end

function fz_syst(z)
real(8) fz_syst,f_syst
complex(8) z
fz_syst=f_syst(dreal(z),dimag(z))
end

function f_syst(xx,yy)
!множитель перед omega в первом уравнении системы
use pgmod2
real(8) f_syst,xx,yy
if (nfun==18.or.nfun==19.or.nfun==22.or.nfun==25) f_syst=d1
if (nfun==20) f_syst=-k_helm**2
if (nfun==21) f_syst=k_helm**2
if (nfun==23) f_syst=d2*xx
if (nfun==24) f_syst=-(xx**2+yy**2)

!f_syst=-d1
!f_syst=gs%const%k_helm**2
end

function fz_v(z,k)
real(8) fz_v,f_v
complex(8) z
integer(4) k
fz_v=f_v(dreal(z),dimag(z),k)
end

function f_v(xx,yy,k)
!компоненты скорости в уравнении переноса
use pgmod2
integer(4) k
real(8) f_v,xx,yy,rrb
if (nfun==8) then
  !(f(1)=0; f(2)=-k_helm**2; f(3)=0; f(4)=0)
  if (k==2) then 
    f_v=-k_helm**2
  endif
endif
if (nfun==13) then
  !(f(1)=(d2+xx**2*(yy-d1))*dsin(yy); f(2)=-y; f(3)=0; f(4)=0)
  if (k==1) then 
    f_v=(d2+xx**2*(yy-d1))*dsin(yy)
  elseif (k==2) then 
    f_v=-yy
  endif
endif
if (nfun==26.or.nfun==27) then
  if (nmain==8) then
    if (k==1) then 
      f_v=d1
    else 
      f_v=d0 
    endif
  endif
  if (nmain==25) then
    !(f(1)=0; f(2)=0; f(3)=d1; f(4)=d0)
    if (k==3) then 
      f_v=d1
    endif
  endif
endif
if (nfun==28) then
  if (nmain==8) then
    if (k==1) then 
      f_v=xx
    else 
      f_v=yy**2 
    endif
  endif
  if (nmain==25) then
    !(f(1)=d2*(d2-xx**2-yy**3); f(2)=0; f(3)=x; f(4)=y**2)
    if (k==1) then 
      f_v=d2*(d2-xx**2-yy**3)
    elseif (k==3) then 
      f_v=xx
    elseif (k==4) then 
      f_v=yy**2 
    endif
  endif
endif
if (nfun==29) then
  if (k==1) then 
    f_v=xx
  else 
    f_v=yy 
  endif
endif
if (nfun==30) then
  if (k==1) then 
    f_v=Pe
  else 
    f_v=d0
  endif
endif
if (nfun==31) then
  if (k==1) then 
    f_v=Pe*4.0d0*yy*(d1-yy)
  else 
    f_v=d0
  endif
endif
if (nfun==32) then
  rrb=s_darci/(xx**2+yy**2)
  if (k==1) then 
    f_v=xx*rrb
  else 
    f_v=yy*rrb
  endif
endif
if (nfun==37) then
  if (k==1) then 
    f_v=d0
  else 
    f_v=xx*yy
  endif
endif
if (nfun==38.or.nfun==39) then
  if (k==1) then 
    f_v=d0
  else 
    f_v=-k_oss/yy
  endif
endif
if (nfun==42) then
  if (k==1) then 
    f_v=4.0d0*yy*(6.0d0-xx**4)
  elseif (k==2) then 
    f_v=d1
  elseif (k==3) then 
    f_v=-xx**2
  elseif (k==4) then 
    f_v=xx
  elseif (k==5) then 
    f_v=-yy
  elseif (k==6) then 
    f_v=xx**3
  elseif (k==7) then 
    f_v=-xx**2*yy
  endif
endif
if (nfun==43) then
  if (k==1) then 
    f_v=-4.0d0*xx*yy
  elseif (k==4) then 
    f_v=yy
  elseif (k==5) then 
    f_v=xx
  endif
endif
if (nfun==44) then
  if (k==4) then 
    f_v=yy
  elseif (k==5) then 
    f_v=-xx
  endif
endif
if (nfun==45) then
  !(f(1)=x**2*y**2*cos(x); f(2)=x; f(3)=cos(x); f(4)=sin(x); f(5)=y; f(6)=x**2*y; f(7)=x*y)
  if (k==1) then 
    f_v=xx**2*yy**2*cos(xx)
  elseif (k==2) then 
    f_v=xx
  elseif (k==3) then 
    f_v=dcos(xx)
  elseif (k==4) then 
    f_v=dsin(xx)
  elseif (k==5) then 
    f_v=yy
  elseif (k==6) then 
    f_v=xx**2*yy
  elseif (k==7) then 
    f_v=xx*yy
  endif
endif
if (nfun==47) then
  f_v=d0
  if (k==1) then 
    if (xx>d5) f_v=xx
  else 
    if (xx<d5) f_v=yy
  endif
endif
if (nfun==48) then
  f_v=d0
  if (k==1) then 
    if (xx>d5) f_v=yy*(dsin(xx)-xx*dcos(xx))
  elseif (k==4) then 
    if (xx>d5) f_v=xx
  elseif (k==5) then
    if (xx<d5) f_v=yy
  endif
endif
if (nfun==49) then
  !(f(1)=y**2*cos(x)**2-y*sin(x); f(2)=x; f(3)=-y*cos(x); f(4)=-x*y)
  if (k==1) then 
    f_v=(yy*dcos(xx))**2-yy*dsin(xx)
  elseif (k==2) then 
    f_v=xx
  elseif (k==3) then 
    f_v=-yy*dcos(xx)
  elseif (k==4) then 
    f_v=-xx*yy
  endif
endif
end

function psi_inf(xx,yy,nf,nder,ia)
use pgmod2
real(8) psi_inf,xx,yy
complex(8) z,w
integer(4) nf,nder,ia
psi_inf=d0
if (nmain==22) then
  if (nfun==40) then
    if (nf==1) then
      if (nder==0) then
        psi_inf=yy
      elseif (nder==1) then
        psi_inf=d1
      endif
	endif
  elseif (nfun==41) then
    if (nf==1) then
	  !попытка сделать обтекание с вихрем
	  !!!не получилось
      z=dcmplx(xx,yy)
      w=c1/ii*cdlog(z-ii*d2)
	  psi_inf=dimag(w)
	endif
  endif
endif
return 
ia=ia
end

subroutine get_fune(xx,yy,psi1,u,v)
!аналитическое решение уравнения лапласа в области
!обтекание окружности безграничным потоком
use pgmod2
real(8) xx,yy,psi1,u,v
complex(8) z,w
z=dcmplx(xx,yy)
w=z+c1/z
psi1=dimag(w)
w=c1-c1/z**2
u=dreal(w)
v=-dimag(w)
end

        !COMPILER-GENERATED INTERFACE MODULE: Fri Apr 05 15:15:14 2024
        MODULE SPLINE_RAZRYV__genmod
          INTERFACE 
            SUBROUTINE SPLINE_RAZRYV(XX,YY,BB,CC,NR,NMAX,RAZRYV)
              INTEGER(KIND=4) :: NMAX
              REAL(KIND=8) :: XX(NMAX)
              REAL(KIND=8) :: YY(NMAX)
              REAL(KIND=8) :: BB(NMAX)
              REAL(KIND=8) :: CC(4,NMAX)
              INTEGER(KIND=4) :: NR
              INTEGER(KIND=4) :: RAZRYV(NMAX)
            END SUBROUTINE SPLINE_RAZRYV
          END INTERFACE 
        END MODULE SPLINE_RAZRYV__genmod

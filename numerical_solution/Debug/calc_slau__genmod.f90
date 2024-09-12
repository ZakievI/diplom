        !COMPILER-GENERATED INTERFACE MODULE: Mon Sep 09 10:15:05 2024
        MODULE CALC_SLAU__genmod
          INTERFACE 
            SUBROUTINE CALC_SLAU(A,B,N,NMAX,X,MODE,FIND_NORM,NORM)
              INTEGER(KIND=4) :: NMAX
              INTEGER(KIND=4) :: N
              REAL(KIND=8) :: A(NMAX,NMAX)
              REAL(KIND=8) :: B(NMAX)
              REAL(KIND=8) :: X(NMAX)
              INTEGER(KIND=4) :: MODE
              LOGICAL(KIND=4) :: FIND_NORM
              REAL(KIND=8) :: NORM
            END SUBROUTINE CALC_SLAU
          END INTERFACE 
        END MODULE CALC_SLAU__genmod

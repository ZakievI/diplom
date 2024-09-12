        !COMPILER-GENERATED INTERFACE MODULE: Fri Apr 05 15:15:09 2024
        MODULE CALC_SLAU_LS__genmod
          INTERFACE 
            SUBROUTINE CALC_SLAU_LS(A,B,M,N,NMAX,X,MODE,FIND_NORM,NORM)
              INTEGER(KIND=4) :: NMAX
              INTEGER(KIND=4) :: N
              REAL(KIND=8) :: A(NMAX,N)
              REAL(KIND=8) :: B(NMAX)
              INTEGER(KIND=4) :: M
              REAL(KIND=8) :: X(N)
              INTEGER(KIND=4) :: MODE
              LOGICAL(KIND=4) :: FIND_NORM
              REAL(KIND=8) :: NORM
            END SUBROUTINE CALC_SLAU_LS
          END INTERFACE 
        END MODULE CALC_SLAU_LS__genmod

        !COMPILER-GENERATED INTERFACE MODULE: Fri Apr 05 15:15:09 2024
        MODULE SLAU_VCL__genmod
          INTERFACE 
            SUBROUTINE SLAU_VCL(A,B,N,NMAX,X,TOFLOAT,MODE)
              INTEGER(KIND=4) :: NMAX
              REAL(KIND=8) :: A(NMAX,NMAX)
              REAL(KIND=8) :: B(NMAX)
              INTEGER(KIND=4) :: N
              REAL(KIND=8) :: X(NMAX)
              INTEGER(KIND=4) :: TOFLOAT
              INTEGER(KIND=4) :: MODE
            END SUBROUTINE SLAU_VCL
          END INTERFACE 
        END MODULE SLAU_VCL__genmod

        !COMPILER-GENERATED INTERFACE MODULE: Sat Dec 07 19:21:23 2024
        MODULE LS_APPROX2__genmod
          INTERFACE 
            SUBROUTINE LS_APPROX2(N,M,X,Y,W,CC,FCN_LS_APPROX)
              INTEGER(KIND=4) :: M
              INTEGER(KIND=4) :: N
              REAL(KIND=8) :: X(M)
              REAL(KIND=8) :: Y(M)
              REAL(KIND=8) :: W(M)
              REAL(KIND=8) :: CC(N)
              REAL(KIND=8) :: FCN_LS_APPROX
              EXTERNAL FCN_LS_APPROX
            END SUBROUTINE LS_APPROX2
          END INTERFACE 
        END MODULE LS_APPROX2__genmod

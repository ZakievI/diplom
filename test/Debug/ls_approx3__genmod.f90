        !COMPILER-GENERATED INTERFACE MODULE: Sat Dec 07 19:21:23 2024
        MODULE LS_APPROX3__genmod
          INTERFACE 
            SUBROUTINE LS_APPROX3(N,M,X,Y,CC,FCN_LS_APPROX)
              INTEGER(KIND=4) :: M
              INTEGER(KIND=4) :: N
              REAL(KIND=8) ,TARGET :: X(M)
              REAL(KIND=8) ,TARGET :: Y(M)
              REAL(KIND=8) :: CC(N)
              REAL(KIND=8) :: FCN_LS_APPROX
              EXTERNAL FCN_LS_APPROX
            END SUBROUTINE LS_APPROX3
          END INTERFACE 
        END MODULE LS_APPROX3__genmod

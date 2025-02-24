        !COMPILER-GENERATED INTERFACE MODULE: Sat Dec 07 19:21:06 2024
        MODULE BM_CALC_SLAU_BLOCK_DIAGONAL2__genmod
          INTERFACE 
            SUBROUTINE BM_CALC_SLAU_BLOCK_DIAGONAL2(BM,RSH,X,LAM,IMAX,  &
     &EPS,MAXITER_CONVERGENCE,SYSTEM_I,WRITE_ITER)
              USE SLAU_BLOCK
              TYPE (MAIN_BLOCK_MATRIX) ,TARGET :: BM
              REAL(KIND=8) :: RSH(BM%NU_ALL)
              REAL(KIND=8) ,TARGET :: X(BM%NX_ALL)
              REAL(KIND=8) :: LAM
              INTEGER(KIND=4) :: IMAX
              REAL(KIND=8) :: EPS
              INTEGER(KIND=4) :: MAXITER_CONVERGENCE
              INTEGER(KIND=4) :: SYSTEM_I
              LOGICAL(KIND=4) :: WRITE_ITER
              EXTERNAL WRITE_ITER
            END SUBROUTINE BM_CALC_SLAU_BLOCK_DIAGONAL2
          END INTERFACE 
        END MODULE BM_CALC_SLAU_BLOCK_DIAGONAL2__genmod

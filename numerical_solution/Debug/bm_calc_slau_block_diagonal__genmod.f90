        !COMPILER-GENERATED INTERFACE MODULE: Sat Dec 07 19:21:05 2024
        MODULE BM_CALC_SLAU_BLOCK_DIAGONAL__genmod
          INTERFACE 
            SUBROUTINE BM_CALC_SLAU_BLOCK_DIAGONAL(BM,RSH,X,LAM,IMAX,EPS&
     &,MAXITER_CONVERGENCE)
              USE SLAU_BLOCK
              TYPE (MAIN_BLOCK_MATRIX) ,TARGET :: BM
              REAL(KIND=8) :: RSH(BM%NU_ALL)
              REAL(KIND=8) :: X(BM%NX_ALL)
              REAL(KIND=8) :: LAM
              INTEGER(KIND=4) :: IMAX
              REAL(KIND=8) :: EPS
              INTEGER(KIND=4) :: MAXITER_CONVERGENCE
            END SUBROUTINE BM_CALC_SLAU_BLOCK_DIAGONAL
          END INTERFACE 
        END MODULE BM_CALC_SLAU_BLOCK_DIAGONAL__genmod

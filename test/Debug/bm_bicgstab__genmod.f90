        !COMPILER-GENERATED INTERFACE MODULE: Sat Dec 07 19:21:06 2024
        MODULE BM_BICGSTAB__genmod
          INTERFACE 
            SUBROUTINE BM_BICGSTAB(BM,RSH,X,IMAX,EPS,MAXITER_CONVERGENCE&
     &,WRITE_ITER)
              USE SLAU_BLOCK
              TYPE (MAIN_BLOCK_MATRIX) :: BM
              REAL(KIND=8) :: RSH(BM%NU_ALL)
              REAL(KIND=8) ,TARGET :: X(BM%NX_ALL)
              INTEGER(KIND=4) :: IMAX
              REAL(KIND=8) :: EPS
              INTEGER(KIND=4) :: MAXITER_CONVERGENCE
              LOGICAL(KIND=4) :: WRITE_ITER
              EXTERNAL WRITE_ITER
            END SUBROUTINE BM_BICGSTAB
          END INTERFACE 
        END MODULE BM_BICGSTAB__genmod

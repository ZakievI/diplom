        !COMPILER-GENERATED INTERFACE MODULE: Fri Apr 05 15:15:02 2024
        MODULE BM_ILUP_SOLVE__genmod
          INTERFACE 
            SUBROUTINE BM_ILUP_SOLVE(BMILU,X)
              USE SLAU_BLOCK
              TYPE (MAIN_BLOCK_MATRIX) ,TARGET :: BMILU
              REAL(KIND=8) :: X(:)
            END SUBROUTINE BM_ILUP_SOLVE
          END INTERFACE 
        END MODULE BM_ILUP_SOLVE__genmod

        !COMPILER-GENERATED INTERFACE MODULE: Sat Dec 07 19:21:06 2024
        MODULE BM_INIT_ILU__genmod
          INTERFACE 
            SUBROUTINE BM_INIT_ILU(BM,BMILU,INIT_P)
              USE SLAU_BLOCK
              TYPE (MAIN_BLOCK_MATRIX) ,TARGET :: BM
              TYPE (MAIN_BLOCK_MATRIX) ,TARGET :: BMILU
              LOGICAL(KIND=4) :: INIT_P
            END SUBROUTINE BM_INIT_ILU
          END INTERFACE 
        END MODULE BM_INIT_ILU__genmod

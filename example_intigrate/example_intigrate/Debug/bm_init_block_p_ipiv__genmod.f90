        !COMPILER-GENERATED INTERFACE MODULE: Fri Apr 05 15:15:02 2024
        MODULE BM_INIT_BLOCK_P_IPIV__genmod
          INTERFACE 
            SUBROUTINE BM_INIT_BLOCK_P_IPIV(BM,IPIV)
              USE SLAU_BLOCK
              TYPE (MAIN_BLOCK_MATRIX) ,TARGET :: BM
              INTEGER(KIND=4) :: IPIV(:)
            END SUBROUTINE BM_INIT_BLOCK_P_IPIV
          END INTERFACE 
        END MODULE BM_INIT_BLOCK_P_IPIV__genmod

        !COMPILER-GENERATED INTERFACE MODULE: Fri Apr 05 15:15:02 2024
        MODULE BM_ILU_SUM__genmod
          INTERFACE 
            SUBROUTINE BM_ILU_SUM(BMILU,B,BILU,K1,J,IMAX,SUM)
              USE SLAU_BLOCK
              TYPE (MAIN_BLOCK_MATRIX) ,TARGET :: BMILU
              TYPE (BLOCK_MATRIX) ,POINTER :: B
              TYPE (BLOCK_MATRIX) ,POINTER :: BILU
              INTEGER(KIND=4) :: K1
              INTEGER(KIND=4) :: J
              INTEGER(KIND=4) :: IMAX
              REAL(KIND=8) :: SUM
            END SUBROUTINE BM_ILU_SUM
          END INTERFACE 
        END MODULE BM_ILU_SUM__genmod

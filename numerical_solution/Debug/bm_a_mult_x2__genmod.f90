        !COMPILER-GENERATED INTERFACE MODULE: Mon Sep 09 10:15:03 2024
        MODULE BM_A_MULT_X2__genmod
          INTERFACE 
            SUBROUTINE BM_A_MULT_X2(ALF,BM,X,BET,BB,RES,BUFF)
              USE SLAU_BLOCK
              REAL(KIND=8) :: ALF
              TYPE (MAIN_BLOCK_MATRIX) ,TARGET :: BM
              REAL(KIND=8) :: X(:)
              REAL(KIND=8) :: BET
              REAL(KIND=8) :: BB(:)
              REAL(KIND=8) :: RES(:)
              REAL(KIND=8) :: BUFF(:)
            END SUBROUTINE BM_A_MULT_X2
          END INTERFACE 
        END MODULE BM_A_MULT_X2__genmod

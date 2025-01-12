        !COMPILER-GENERATED INTERFACE MODULE: Sat Dec 07 19:21:05 2024
        MODULE BM_A_MULT_X__genmod
          INTERFACE 
            SUBROUTINE BM_A_MULT_X(ALF,BM,X,RES,BUFF)
              USE SLAU_BLOCK
              REAL(KIND=8) :: ALF
              TYPE (MAIN_BLOCK_MATRIX) :: BM
              REAL(KIND=8) :: X(:)
              REAL(KIND=8) :: RES(:)
              REAL(KIND=8) :: BUFF(:)
            END SUBROUTINE BM_A_MULT_X
          END INTERFACE 
        END MODULE BM_A_MULT_X__genmod

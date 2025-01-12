        !COMPILER-GENERATED INTERFACE MODULE: Sat Dec 07 19:21:05 2024
        MODULE BM_GET_U_IJ__genmod
          INTERFACE 
            FUNCTION BM_GET_U_IJ(BMILU,I,J,U_IJ) RESULT(RES)
              USE SLAU_BLOCK
              TYPE (MAIN_BLOCK_MATRIX) ,TARGET :: BMILU
              INTEGER(KIND=4) :: I
              INTEGER(KIND=4) :: J
              REAL(KIND=8) :: U_IJ
              LOGICAL(KIND=4) :: RES
            END FUNCTION BM_GET_U_IJ
          END INTERFACE 
        END MODULE BM_GET_U_IJ__genmod

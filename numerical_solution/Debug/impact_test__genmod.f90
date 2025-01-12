        !COMPILER-GENERATED INTERFACE MODULE: Sun Dec 22 15:56:17 2024
        MODULE IMPACT_TEST__genmod
          INTERFACE 
            SUBROUTINE IMPACT_TEST(N,Y,Y_OUT,NUM,DLT)
              INTEGER(KIND=4) :: NUM
              INTEGER(KIND=4) :: N
              REAL(KIND=8) :: Y(N)
              REAL(KIND=8) :: Y_OUT(NUM,2)
              REAL(KIND=8) :: DLT
            END SUBROUTINE IMPACT_TEST
          END INTERFACE 
        END MODULE IMPACT_TEST__genmod

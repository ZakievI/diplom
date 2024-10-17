        !COMPILER-GENERATED INTERFACE MODULE: Thu Oct 17 21:12:12 2024
        MODULE IMPACT_TEST__genmod
          INTERFACE 
            SUBROUTINE IMPACT_TEST(N,Y,Y_OUT,K1,NUM)
              INTEGER(KIND=4) :: NUM
              INTEGER(KIND=4) :: N
              REAL(KIND=8) :: Y(N)
              REAL(KIND=8) :: Y_OUT(NUM,2)
              INTEGER(KIND=4) :: K1
            END SUBROUTINE IMPACT_TEST
          END INTERFACE 
        END MODULE IMPACT_TEST__genmod

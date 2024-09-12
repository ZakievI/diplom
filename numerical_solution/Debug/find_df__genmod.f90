        !COMPILER-GENERATED INTERFACE MODULE: Mon Sep 09 10:15:06 2024
        MODULE FIND_DF__genmod
          INTERFACE 
            FUNCTION FIND_DF(Z1,Z2,DS,DFDL,ITR,IB,K1,K2,IS_BOUND,MODE)
              COMPLEX(KIND=8) :: Z1
              COMPLEX(KIND=8) :: Z2
              REAL(KIND=8) :: DS
              REAL(KIND=8) :: DFDL
              EXTERNAL DFDL
              INTEGER(KIND=4) :: ITR
              INTEGER(KIND=4) :: IB
              INTEGER(KIND=4) :: K1
              INTEGER(KIND=4) :: K2
              LOGICAL(KIND=4) :: IS_BOUND
              INTEGER(KIND=4) :: MODE
              REAL(KIND=8) :: FIND_DF
            END FUNCTION FIND_DF
          END INTERFACE 
        END MODULE FIND_DF__genmod

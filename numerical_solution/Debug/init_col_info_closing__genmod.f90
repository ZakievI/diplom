        !COMPILER-GENERATED INTERFACE MODULE: Sat Nov 16 21:00:31 2024
        MODULE INIT_COL_INFO_CLOSING__genmod
          INTERFACE 
            SUBROUTINE INIT_COL_INFO_CLOSING(A,CI,B1,BL1,K1,B2,BL2,K2)
              USE PGMOD
              TYPE (TAREA) :: A
              TYPE (COL_INFO) ,TARGET :: CI
              TYPE (TBOUND) :: B1
              TYPE (TBOUNDLINE) :: BL1
              INTEGER(KIND=4) :: K1(:)
              TYPE (TBOUND) :: B2
              TYPE (TBOUNDLINE) :: BL2
              INTEGER(KIND=4) :: K2(:)
            END SUBROUTINE INIT_COL_INFO_CLOSING
          END INTERFACE 
        END MODULE INIT_COL_INFO_CLOSING__genmod

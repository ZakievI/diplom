        !COMPILER-GENERATED INTERFACE MODULE: Sat Dec 07 19:21:40 2024
        MODULE INIT_IND_PSIOM__genmod
          INTERFACE 
            SUBROUTINE INIT_IND_PSIOM(B,BL,IND,K)
              USE PGMOD
              TYPE (TBOUND) :: B
              TYPE (TBOUNDLINE) :: BL
              LOGICAL(KIND=4) :: IND(GS%M%NX)
              INTEGER(KIND=4) :: K
            END SUBROUTINE INIT_IND_PSIOM
          END INTERFACE 
        END MODULE INIT_IND_PSIOM__genmod

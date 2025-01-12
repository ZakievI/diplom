        !COMPILER-GENERATED INTERFACE MODULE: Sat Dec 07 19:21:39 2024
        MODULE INIT_BOUNDLINE_GU_GEN__genmod
          INTERFACE 
            SUBROUTINE INIT_BOUNDLINE_GU_GEN(GU,BNDF,IS_DIRECT,IA,IBND, &
     &IBNDL,N,C0_,SECOND_BND,BNDF2,C)
              USE PGMOD
              INTEGER(KIND=4) :: N
              TYPE (TBOUNDLINE_GU) :: GU
              INTEGER(KIND=4) :: BNDF
              LOGICAL(KIND=4) :: IS_DIRECT
              INTEGER(KIND=4) :: IA
              INTEGER(KIND=4) :: IBND
              INTEGER(KIND=4) :: IBNDL
              REAL(KIND=8) :: C0_
              LOGICAL(KIND=4) :: SECOND_BND(N)
              INTEGER(KIND=4) :: BNDF2(N)
              REAL(KIND=8) :: C(N)
            END SUBROUTINE INIT_BOUNDLINE_GU_GEN
          END INTERFACE 
        END MODULE INIT_BOUNDLINE_GU_GEN__genmod

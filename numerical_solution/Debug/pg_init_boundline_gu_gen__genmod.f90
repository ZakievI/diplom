        !COMPILER-GENERATED INTERFACE MODULE: Sat Nov 16 21:00:30 2024
        MODULE PG_INIT_BOUNDLINE_GU_GEN__genmod
          INTERFACE 
            SUBROUTINE PG_INIT_BOUNDLINE_GU_GEN(IGU,BNDF,IS_DIRECT,IA,  &
     &IBND,IBNDL,N,C0_,SECOND_BND,BNDF2,C)
              INTEGER(KIND=4) :: N
              INTEGER(KIND=4) :: IGU
              INTEGER(KIND=4) :: BNDF
              LOGICAL(KIND=4) :: IS_DIRECT
              INTEGER(KIND=4) :: IA
              INTEGER(KIND=4) :: IBND
              INTEGER(KIND=4) :: IBNDL
              REAL(KIND=8) :: C0_
              LOGICAL(KIND=4) :: SECOND_BND(N)
              INTEGER(KIND=4) :: BNDF2(N)
              REAL(KIND=8) :: C(N)
            END SUBROUTINE PG_INIT_BOUNDLINE_GU_GEN
          END INTERFACE 
        END MODULE PG_INIT_BOUNDLINE_GU_GEN__genmod

        !COMPILER-GENERATED INTERFACE MODULE: Mon Sep 09 10:15:06 2024
        MODULE GET_TRBOUND__genmod
          INTERFACE 
            SUBROUTINE GET_TRBOUND(TR,N,NTR,NPE,TRBOUND)
              INTEGER(KIND=4) :: NPE
              INTEGER(KIND=4) :: NTR
              INTEGER(KIND=4) :: TR(NPE,NTR)
              INTEGER(KIND=4) :: N
              LOGICAL(KIND=4) :: TRBOUND(NPE,NTR)
            END SUBROUTINE GET_TRBOUND
          END INTERFACE 
        END MODULE GET_TRBOUND__genmod

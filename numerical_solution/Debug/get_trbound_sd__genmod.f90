        !COMPILER-GENERATED INTERFACE MODULE: Sat Dec 07 19:21:11 2024
        MODULE GET_TRBOUND_SD__genmod
          INTERFACE 
            SUBROUTINE GET_TRBOUND_SD(TR,SD,N,NTR,NPE,TRBOUND,USE_SD)
              INTEGER(KIND=4) :: NPE
              INTEGER(KIND=4) :: NTR
              INTEGER(KIND=4) :: TR(NPE,NTR)
              INTEGER(KIND=2) :: SD(NTR)
              INTEGER(KIND=4) :: N
              INTEGER(KIND=2) :: TRBOUND(NPE,NTR)
              LOGICAL(KIND=4) :: USE_SD
            END SUBROUTINE GET_TRBOUND_SD
          END INTERFACE 
        END MODULE GET_TRBOUND_SD__genmod

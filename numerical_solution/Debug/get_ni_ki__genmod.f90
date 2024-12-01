        !COMPILER-GENERATED INTERFACE MODULE: Sat Nov 16 21:00:40 2024
        MODULE GET_NI_KI__genmod
          INTERFACE 
            SUBROUTINE GET_NI_KI(NEED_CASH,NI,NF,ONLY_NUMERICAL_INT,C,  &
     &NF_DUAL)
              USE PGMOD
              LOGICAL(KIND=4) :: NEED_CASH(0:)
              INTEGER(KIND=4) :: NI
              INTEGER(KIND=4) :: NF(20)
              LOGICAL(KIND=4) :: ONLY_NUMERICAL_INT(20)
              TYPE (TCASHINTEGRAL) :: C
              LOGICAL(KIND=4) :: NF_DUAL(20)
            END SUBROUTINE GET_NI_KI
          END INTERFACE 
        END MODULE GET_NI_KI__genmod

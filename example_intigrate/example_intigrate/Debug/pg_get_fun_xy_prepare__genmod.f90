        !COMPILER-GENERATED INTERFACE MODULE: Fri Apr 05 21:57:13 2024
        MODULE PG_GET_FUN_XY_PREPARE__genmod
          INTERFACE 
            SUBROUTINE PG_GET_FUN_XY_PREPARE(XX,YY,MODE,IN_BOUND,IBND,S,&
     &TT,BN)
              USE PGMOD
              REAL(KIND=8) :: XX
              REAL(KIND=8) :: YY
              INTEGER(KIND=4) :: MODE
              INTEGER(KIND=4) :: IN_BOUND
              INTEGER(KIND=4) :: IBND
              REAL(KIND=8) :: S
              REAL(KIND=8) :: TT
              TYPE (TBOUND_NEAR) :: BN
            END SUBROUTINE PG_GET_FUN_XY_PREPARE
          END INTERFACE 
        END MODULE PG_GET_FUN_XY_PREPARE__genmod

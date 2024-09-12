        !COMPILER-GENERATED INTERFACE MODULE: Fri Apr 05 21:57:13 2024
        MODULE PG_GET_FUN_XY_SOLVE__genmod
          INTERFACE 
            FUNCTION PG_GET_FUN_XY_SOLVE(XX,YY,NF,XDER,YDER,IN_BOUND,   &
     &IBND,S,TT,BN)
              USE PGMOD
              REAL(KIND=8) :: XX
              REAL(KIND=8) :: YY
              INTEGER(KIND=4) :: NF
              REAL(KIND=8) :: XDER
              REAL(KIND=8) :: YDER
              INTEGER(KIND=4) :: IN_BOUND
              INTEGER(KIND=4) :: IBND
              REAL(KIND=8) :: S
              REAL(KIND=8) :: TT
              TYPE (TBOUND_NEAR) :: BN
              REAL(KIND=8) :: PG_GET_FUN_XY_SOLVE
            END FUNCTION PG_GET_FUN_XY_SOLVE
          END INTERFACE 
        END MODULE PG_GET_FUN_XY_SOLVE__genmod

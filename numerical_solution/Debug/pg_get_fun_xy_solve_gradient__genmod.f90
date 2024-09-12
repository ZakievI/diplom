        !COMPILER-GENERATED INTERFACE MODULE: Mon Sep 09 10:15:22 2024
        MODULE PG_GET_FUN_XY_SOLVE_GRADIENT__genmod
          INTERFACE 
            RECURSIVE SUBROUTINE PG_GET_FUN_XY_SOLVE_GRADIENT(XX,YY,NF, &
     &XDER,YDER,IN_BOUND,IBND,S,TT,FF,BN)
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
              REAL(KIND=8) :: FF(2)
              TYPE (TBOUND_NEAR) :: BN
            END SUBROUTINE PG_GET_FUN_XY_SOLVE_GRADIENT
          END INTERFACE 
        END MODULE PG_GET_FUN_XY_SOLVE_GRADIENT__genmod

        !COMPILER-GENERATED INTERFACE MODULE: Mon Sep 09 10:15:22 2024
        MODULE PG_GET_BOUND_S_BT__genmod
          INTERFACE 
            FUNCTION PG_GET_BOUND_S_BT(XX,YY,IBND,SMIN,TTMIN,ONLY_BT1,BN&
     &)
              USE PGMOD
              REAL(KIND=8) :: XX
              REAL(KIND=8) :: YY
              INTEGER(KIND=4) :: IBND
              REAL(KIND=8) :: SMIN
              REAL(KIND=8) :: TTMIN
              LOGICAL(KIND=4) :: ONLY_BT1
              TYPE (TBOUND_NEAR) :: BN
              INTEGER(KIND=4) :: PG_GET_BOUND_S_BT
            END FUNCTION PG_GET_BOUND_S_BT
          END INTERFACE 
        END MODULE PG_GET_BOUND_S_BT__genmod

        !COMPILER-GENERATED INTERFACE MODULE: Sat Dec 07 19:21:44 2024
        MODULE DISTANCETOBOUND__genmod
          INTERFACE 
            FUNCTION DISTANCETOBOUND(XX,YY,IB,SMIN,TTMIN,BN)
              USE PGMOD
              REAL(KIND=8) :: XX
              REAL(KIND=8) :: YY
              INTEGER(KIND=4) :: IB
              REAL(KIND=8) :: SMIN
              REAL(KIND=8) :: TTMIN
              TYPE (TBOUND_NEAR) :: BN
              REAL(KIND=8) :: DISTANCETOBOUND
            END FUNCTION DISTANCETOBOUND
          END INTERFACE 
        END MODULE DISTANCETOBOUND__genmod

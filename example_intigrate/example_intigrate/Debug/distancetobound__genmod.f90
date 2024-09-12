        !COMPILER-GENERATED INTERFACE MODULE: Fri Apr 05 21:57:12 2024
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

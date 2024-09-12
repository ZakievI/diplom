        !COMPILER-GENERATED INTERFACE MODULE: Fri Apr 05 15:15:14 2024
        MODULE DISTANCETOLINE2__genmod
          INTERFACE 
            FUNCTION DISTANCETOLINE2(XX,YY,X,Y,NPANEL,J,SMIN,PX_MIN,    &
     &PY_MIN,FIRST_I)
              INTEGER(KIND=4) :: NPANEL
              REAL(KIND=8) :: XX
              REAL(KIND=8) :: YY
              REAL(KIND=8) :: X(NPANEL+1)
              REAL(KIND=8) :: Y(NPANEL+1)
              INTEGER(KIND=4) :: J
              REAL(KIND=8) :: SMIN
              REAL(KIND=8) :: PX_MIN
              REAL(KIND=8) :: PY_MIN
              INTEGER(KIND=4) :: FIRST_I
              REAL(KIND=8) :: DISTANCETOLINE2
            END FUNCTION DISTANCETOLINE2
          END INTERFACE 
        END MODULE DISTANCETOLINE2__genmod

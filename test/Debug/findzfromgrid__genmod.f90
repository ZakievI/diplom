        !COMPILER-GENERATED INTERFACE MODULE: Sat Dec 07 19:21:23 2024
        MODULE FINDZFROMGRID__genmod
          INTERFACE 
            FUNCTION FINDZFROMGRID(XX,YY,X0,Y0,DX,DY,NX,NY,FF,EMPTYVAL) &
     & RESULT(RES)
              INTEGER(KIND=4) :: NY
              INTEGER(KIND=4) :: NX
              REAL(KIND=8) :: XX
              REAL(KIND=8) :: YY
              REAL(KIND=8) :: X0
              REAL(KIND=8) :: Y0
              REAL(KIND=8) :: DX
              REAL(KIND=8) :: DY
              REAL(KIND=8) :: FF(NX+1,NY+1)
              REAL(KIND=8) :: EMPTYVAL
              REAL(KIND=8) :: RES
            END FUNCTION FINDZFROMGRID
          END INTERFACE 
        END MODULE FINDZFROMGRID__genmod

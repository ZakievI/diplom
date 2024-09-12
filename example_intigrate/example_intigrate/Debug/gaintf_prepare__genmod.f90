        !COMPILER-GENERATED INTERFACE MODULE: Fri Apr 05 21:57:21 2024
        MODULE GAINTF_PREPARE__genmod
          INTERFACE 
            SUBROUTINE GAINTF_PREPARE(NN,X1,Y1,X2,Y2,DS,XX,YY,MODE,IFS)
              USE PGMOD
              INTEGER(KIND=4) :: NN
              REAL(KIND=8) :: X1
              REAL(KIND=8) :: Y1
              REAL(KIND=8) :: X2
              REAL(KIND=8) :: Y2
              REAL(KIND=8) :: DS
              REAL(KIND=8) ,TARGET :: XX(NN)
              REAL(KIND=8) ,TARGET :: YY(NN)
              INTEGER(KIND=4) :: MODE
              TYPE (INTF_STRUCT) ,TARGET :: IFS
            END SUBROUTINE GAINTF_PREPARE
          END INTERFACE 
        END MODULE GAINTF_PREPARE__genmod

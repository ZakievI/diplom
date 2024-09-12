        !COMPILER-GENERATED INTERFACE MODULE: Fri Apr 05 21:57:13 2024
        MODULE GET_FUN_XY_SOLVE_GRADIENT_BI__genmod
          INTERFACE 
            SUBROUTINE GET_FUN_XY_SOLVE_GRADIENT_BI(XX,YY,NF,XDER,YDER, &
     &BI,FF)
              USE PGMOD
              REAL(KIND=8) :: XX
              REAL(KIND=8) :: YY
              INTEGER(KIND=4) :: NF
              REAL(KIND=8) :: XDER
              REAL(KIND=8) :: YDER
              TYPE (TBOUND_INFO) :: BI
              REAL(KIND=8) :: FF(2)
            END SUBROUTINE GET_FUN_XY_SOLVE_GRADIENT_BI
          END INTERFACE 
        END MODULE GET_FUN_XY_SOLVE_GRADIENT_BI__genmod

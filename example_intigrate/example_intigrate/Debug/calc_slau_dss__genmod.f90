        !COMPILER-GENERATED INTERFACE MODULE: Fri Apr 05 15:15:08 2024
        MODULE CALC_SLAU_DSS__genmod
          INTERFACE 
            SUBROUTINE CALC_SLAU_DSS(VALUES,COLUMNS,ROWINDEX,RHS,NROWS, &
     &NCOLS,NNONZEROS,SOLUTION)
              INTEGER(KIND=4) :: NNONZEROS
              INTEGER(KIND=4) :: NCOLS
              INTEGER(KIND=4) :: NROWS
              REAL(KIND=8) :: VALUES(NNONZEROS)
              INTEGER(KIND=4) :: COLUMNS(NNONZEROS)
              INTEGER(KIND=4) :: ROWINDEX(NROWS+1)
              REAL(KIND=8) :: RHS(NROWS)
              REAL(KIND=8) :: SOLUTION(NCOLS)
            END SUBROUTINE CALC_SLAU_DSS
          END INTERFACE 
        END MODULE CALC_SLAU_DSS__genmod

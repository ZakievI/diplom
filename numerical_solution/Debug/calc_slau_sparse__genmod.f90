        !COMPILER-GENERATED INTERFACE MODULE: Sat Dec 07 19:21:09 2024
        MODULE CALC_SLAU_SPARSE__genmod
          INTERFACE 
            SUBROUTINE CALC_SLAU_SPARSE(VALUES,COLUMNS,ROWINDEX,RHS,    &
     &NROWS,NCOLS,NNONZEROS,SOLUTION,MODE,FIND_NORM,NORM,PARDISO_EPS)
              INTEGER(KIND=4) :: NNONZEROS
              INTEGER(KIND=4) :: NCOLS
              INTEGER(KIND=4) :: NROWS
              REAL(KIND=8) :: VALUES(NNONZEROS)
              INTEGER(KIND=4) :: COLUMNS(NNONZEROS)
              INTEGER(KIND=4) :: ROWINDEX(NROWS+1)
              REAL(KIND=8) :: RHS(NROWS)
              REAL(KIND=8) :: SOLUTION(NCOLS)
              INTEGER(KIND=4) :: MODE
              LOGICAL(KIND=4) :: FIND_NORM
              REAL(KIND=8) :: NORM
              INTEGER(KIND=4) :: PARDISO_EPS
            END SUBROUTINE CALC_SLAU_SPARSE
          END INTERFACE 
        END MODULE CALC_SLAU_SPARSE__genmod

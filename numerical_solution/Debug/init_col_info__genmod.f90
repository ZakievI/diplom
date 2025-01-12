        !COMPILER-GENERATED INTERFACE MODULE: Sat Dec 07 19:21:32 2024
        MODULE INIT_COL_INFO__genmod
          INTERFACE 
            SUBROUTINE INIT_COL_INFO(A,CI,IND,BUFF_C,N,NEED_FICT)
              USE PGMOD
              INTEGER(KIND=4) :: N
              TYPE (TAREA) :: A
              TYPE (COL_INFO) :: CI
              LOGICAL(KIND=4) :: IND(N)
              INTEGER(KIND=4) :: BUFF_C(N)
              LOGICAL(KIND=4) :: NEED_FICT
            END SUBROUTINE INIT_COL_INFO
          END INTERFACE 
        END MODULE INIT_COL_INFO__genmod

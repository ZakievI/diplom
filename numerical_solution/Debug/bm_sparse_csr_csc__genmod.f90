        !COMPILER-GENERATED INTERFACE MODULE: Sat Dec 07 19:21:05 2024
        MODULE BM_SPARSE_CSR_CSC__genmod
          INTERFACE 
            SUBROUTINE BM_SPARSE_CSR_CSC(SP1,SP2,IU,IU2,M)
              USE SLAU_BLOCK
              TYPE (SPARSE_MATRIX) :: SP1
              TYPE (SPARSE_MATRIX) :: SP2
              INTEGER(KIND=4) :: IU
              INTEGER(KIND=4) :: IU2
              INTEGER(KIND=4) :: M
            END SUBROUTINE BM_SPARSE_CSR_CSC
          END INTERFACE 
        END MODULE BM_SPARSE_CSR_CSC__genmod

        !COMPILER-GENERATED INTERFACE MODULE: Sat Dec 07 19:21:23 2024
        MODULE MATR_VECT_P__genmod
          INTERFACE 
            SUBROUTINE MATR_VECT_P(MM,V,R,N,M,LDA)
              INTEGER(KIND=4) :: M
              INTEGER(KIND=4) :: N
              REAL(KIND=8) :: MM(N,M)
              REAL(KIND=8) :: V(M)
              REAL(KIND=8) :: R(N)
              INTEGER(KIND=4) :: LDA
            END SUBROUTINE MATR_VECT_P
          END INTERFACE 
        END MODULE MATR_VECT_P__genmod

        !COMPILER-GENERATED INTERFACE MODULE: Sat Dec 07 19:21:23 2024
        MODULE GETFOURIERSERIESCOEFF__genmod
          INTERFACE 
            FUNCTION GETFOURIERSERIESCOEFF(N,S,F,SMIN,SMAX,K,MODE)      &
     & RESULT(RES)
              INTEGER(KIND=4) :: N
              REAL(KIND=8) :: S(N)
              REAL(KIND=8) :: F(N)
              REAL(KIND=8) :: SMIN
              REAL(KIND=8) :: SMAX
              INTEGER(KIND=4) :: K
              INTEGER(KIND=4) :: MODE
              REAL(KIND=8) :: RES
            END FUNCTION GETFOURIERSERIESCOEFF
          END INTERFACE 
        END MODULE GETFOURIERSERIESCOEFF__genmod

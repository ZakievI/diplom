        !COMPILER-GENERATED INTERFACE MODULE: Fri Apr 05 15:15:10 2024
        MODULE INTF_DXDY2__genmod
          INTERFACE 
            SUBROUTINE INTF_DXDY2(ZM,TR,N,NTR,NPE,FF,FFKNOW,DFDL,DS)
              INTEGER(KIND=4) :: NPE
              INTEGER(KIND=4) :: NTR
              INTEGER(KIND=4) :: N
              COMPLEX(KIND=8) :: ZM(N)
              INTEGER(KIND=4) :: TR(NPE,NTR)
              REAL(KIND=8) :: FF(N)
              LOGICAL(KIND=4) :: FFKNOW(N)
              REAL(KIND=8) :: DFDL
              EXTERNAL DFDL
              REAL(KIND=8) :: DS
            END SUBROUTINE INTF_DXDY2
          END INTERFACE 
        END MODULE INTF_DXDY2__genmod

        !COMPILER-GENERATED INTERFACE MODULE: Mon Sep 09 10:15:06 2024
        MODULE INTF_DXDY2_SD__genmod
          INTERFACE 
            SUBROUTINE INTF_DXDY2_SD(ZM,TR,N,NTR,NPE,FF,FFKNOW,DFDL,DS, &
     &SD)
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
              INTEGER(KIND=2) :: SD(NTR)
            END SUBROUTINE INTF_DXDY2_SD
          END INTERFACE 
        END MODULE INTF_DXDY2_SD__genmod

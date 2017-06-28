!    This file is part of SRO_EOS.
!
!    SRO_EOS is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    SRO_EOS is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with SRO_EOS.  If not, see <http://www.gnu.org/licenses/>.
!
MODULE Beta_Mod

  USE Kind_Types_Mod, ONLY : DP, I4B
  USE Physical_Constants_Mod, ONLY : ZERO, ONE, TWO, BETA_0, R_2_3

  IMPLICIT NONE

CONTAINS

  SUBROUTINE BETA_SUB (derivatives,y,n,H,DH_DT,DH_DY,D2H_DT2,D2H_DY2,D2H_DTDY,&
          SIGMA,DSIGMA_DT,DSIGMA_DY,D2SIGMA_DT2,D2SIGMA_DTDY,D2SIGMA_DY2,  &
          BETA,DBETA_DN,DBETA_DT,DBETA_DY,D2BETA_DN2,D2BETA_DNDT,D2BETA_DNDY,&
          D2BETA_DT2,D2BETA_DTDY,D2BETA_DY2)

    IMPLICIT NONE

    INTEGER(I4B), INTENT(IN) :: derivatives
    REAL(DP), INTENT(IN) :: y, n
    REAL(DP), INTENT(IN) :: H, DH_DT, DH_DY, D2H_DT2, D2H_DY2, D2H_DTDY
    REAL(DP), INTENT(IN) :: SIGMA, DSIGMA_DT, DSIGMA_DY
    REAL(DP), INTENT(IN) :: D2SIGMA_DT2, D2SIGMA_DY2, D2SIGMA_DTDY
    REAL(DP), INTENT(OUT) :: BETA, DBETA_DY, DBETA_DN, DBETA_DT
    REAL(DP), INTENT(OUT) :: D2BETA_DY2, D2BETA_DN2, D2BETA_DT2
    REAL(DP), INTENT(OUT) :: D2BETA_DNDY,D2BETA_DNDT,D2BETA_DTDY

    BETA = BETA_0*(y*n*SIGMA)**R_2_3

    DBETA_DY = ZERO
    DBETA_DN = ZERO
    DBETA_DT = ZERO

    D2BETA_DY2 = ZERO
    D2BETA_DN2 = ZERO
    D2BETA_DT2 = ZERO

    D2BETA_DNDY = ZERO
    D2BETA_DNDT = ZERO
    D2BETA_DTDY = ZERO

    IF (derivatives == 0) RETURN

    DBETA_DY = R_2_3*BETA*(ONE/y + DSIGMA_DY/SIGMA)
    DBETA_DN = R_2_3*BETA/n
    DBETA_DT = R_2_3*BETA*DH_DT/H

    IF (derivatives == 1) RETURN

    D2BETA_DY2 = DBETA_DY*DBETA_DY/BETA + R_2_3*BETA*( - ONE/y/y &
               + D2SIGMA_DY2/SIGMA - (DSIGMA_DY/SIGMA)**TWO )
    D2BETA_DN2 = DBETA_DN*DBETA_DN/BETA - R_2_3*BETA/n/n
    D2BETA_DT2 = DBETA_DT*DBETA_DT/BETA &
               + R_2_3*BETA*( D2H_DT2/H - (DH_DT/H)**TWO )

    D2BETA_DNDY = R_2_3*DBETA_DY/n
    D2BETA_DNDT = R_2_3*DBETA_DT/n
    D2BETA_DTDY = R_2_3*DBETA_DY*DH_DT/H &
                + R_2_3*BETA*(D2H_DTDY/H - DH_DY*DH_DT/H**TWO )

    IF (derivatives == 2) RETURN

  END SUBROUTINE BETA_SUB

END MODULE Beta_Mod

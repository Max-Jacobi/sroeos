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
MODULE Nuclear_Radius_Mod

  USE Kind_Types_Mod, ONLY : DP, I4B
  USE Physical_Constants_Mod, ONLY : ZERO, ONE, TWO
  USE Surface_Properties_Mod, ONLY : fix_heavy_nuclei_size

  IMPLICIT NONE

CONTAINS

  SUBROUTINE NUCLEAR_RADIUS_SUB (derivatives,u,DELU,DDELU_DU,D2DELU_DU2,   &
          SIGMA,DSIGMA_DT,DSIGMA_DY,D2SIGMA_DT2,D2SIGMA_DY2,D2SIGMA_DTDY,  &
          BETA,DBETA_DY,DBETA_DN,DBETA_DT,D2BETA_DY2,D2BETA_DN2,D2BETA_DT2,&
          D2BETA_DNDY,D2BETA_DNDT,D2BETA_DTDY,R,DR_DU,DR_DN,DR_DT,DR_DY,   &
          D2R_DU2,D2R_DUDN,D2R_DUDT,D2R_DUDY,D2R_DN2,D2R_DNDT,D2R_DNDY,    &
          D2R_DT2,D2R_DTDY,D2R_DY2)

    IMPLICIT NONE

    INTEGER(I4B), INTENT(IN) :: derivatives
    REAL(DP), INTENT(IN) :: u, DELU, DDELU_DU, D2DELU_DU2
    REAL(DP), INTENT(IN) :: SIGMA, DSIGMA_DT, DSIGMA_DY
    REAL(DP), INTENT(IN) :: D2SIGMA_DT2, D2SIGMA_DY2, D2SIGMA_DTDY
    REAL(DP), INTENT(IN) :: BETA, DBETA_DY, DBETA_DN, DBETA_DT
    REAL(DP), INTENT(IN) :: D2BETA_DY2, D2BETA_DN2, D2BETA_DT2
    REAL(DP), INTENT(IN) :: D2BETA_DNDY,D2BETA_DNDT,D2BETA_DTDY
    REAL(DP), INTENT(OUT) :: R, DR_DU, DR_DN, DR_DT, DR_DY
    REAL(DP), INTENT(OUT) :: D2R_DU2, D2R_DUDN, D2R_DUDT, D2R_DUDY
    REAL(DP), INTENT(OUT) :: D2R_DN2, D2R_DNDT, D2R_DNDY
    REAL(DP), INTENT(OUT) :: D2R_DT2, D2R_DTDY, D2R_DY2
    REAL(DP) :: odu

    odu = ONE - u

    R = 9.0D0/TWO*SIGMA/BETA/DELU*u*odu

    DR_DU = ZERO
    DR_DN = ZERO
    DR_DT = ZERO
    DR_DY = ZERO

    D2R_DU2  = ZERO
    D2R_DUDN = ZERO
    D2R_DUDT = ZERO
    D2R_DUDY = ZERO

    D2R_DN2  = ZERO
    D2R_DNDT = ZERO
    D2R_DNDY = ZERO

    D2R_DT2  = ZERO
    D2R_DTDY = ZERO
    D2R_DY2  = ZERO

    IF (derivatives == 0) RETURN

    DR_DU = R*( - DDELU_DU/DELU + (ONE-TWO*u)/u/odu )
    DR_DN = - R*DBETA_DN/BETA
    DR_DT = R*(DSIGMA_DT/SIGMA-DBETA_DT/BETA)
    DR_DY = R*(DSIGMA_DY/SIGMA-DBETA_DY/BETA)

    IF (derivatives == 1) RETURN

    D2R_DU2  = DR_DU*DR_DU/R + R*( (DDELU_DU/DELU)**TWO - D2DELU_DU2/DELU &
             - TWO/(u*odu) - ((ONE-TWO*u)/u/odu)**TWO )
    D2R_DUDN = DR_DN*DR_DU/R
    D2R_DUDT = DR_DT*DR_DU/R
    D2R_DUDY = DR_DY*DR_DU/R

    D2R_DN2  = DR_DN*DR_DN/R &
             - R*( D2BETA_DN2/BETA  - (DBETA_DN/BETA)**TWO )
    D2R_DNDT = DR_DN*DR_DT/R &
             - R*( D2BETA_DNDT/BETA - DBETA_DN*DBETA_DT/BETA**TWO )
    D2R_DNDY = DR_DN*DR_DY/R &
             - R*( D2BETA_DNDY/BETA - DBETA_DN*DBETA_DY/BETA**TWO )

    D2R_DT2  = DR_DT*DR_DT/R &
             + R*( D2SIGMA_DT2/SIGMA - (DSIGMA_DT/SIGMA)**TWO ) &
             - R*( D2BETA_DT2/BETA - (DBETA_DT/BETA)** TWO )
    D2R_DTDY = DR_DT*DR_DY/R &
             + R*( D2SIGMA_DTDY/SIGMA - DSIGMA_DT*DSIGMA_DY/SIGMA**TWO ) &
             - R*( D2BETA_DTDY/BETA - DBETA_DT*DBETA_DY/BETA** TWO )
    D2R_DY2  = DR_DY*DR_DY/R &
             + R*( D2SIGMA_DY2/SIGMA - (DSIGMA_DY/SIGMA)**TWO ) &
             - R*( D2BETA_DY2/BETA - (DBETA_DY/BETA)** TWO )

    IF (derivatives == 2) RETURN

  END SUBROUTINE NUCLEAR_RADIUS_SUB

END MODULE Nuclear_Radius_Mod

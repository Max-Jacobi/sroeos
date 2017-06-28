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
MODULE Sigma_Mod

  USE Kind_Types_Mod, ONLY : DP, I4B
  USE Physical_Constants_Mod, ONLY : ZERO, ONE, TWO
  USE Surface_Properties_Mod, ONLY : Q1 => Surface_q, L1 => Surface_lambda, &
                                    SIGMA_0 => Surface_sigma
  IMPLICIT NONE

CONTAINS

  SUBROUTINE SIGMA_SUB (derivatives,y,T,H,DH_DT,DH_DY,D2H_DT2,D2H_DY2,D2H_DTDY,&
                SIGMA,DSIGMA_DT,DSIGMA_DY,D2SIGMA_DT2,D2SIGMA_DTDY,D2SIGMA_DY2)

    IMPLICIT NONE

    INTEGER(I4B), INTENT(IN) :: derivatives
    REAL(DP), INTENT(IN) :: y, T
    REAL(DP), INTENT(IN) :: H, DH_DT, DH_DY, D2H_DT2, D2H_DY2, D2H_DTDY
    REAL(DP), INTENT(OUT) :: SIGMA, DSIGMA_DT, DSIGMA_DY
    REAL(DP), INTENT(OUT) :: D2SIGMA_DT2, D2SIGMA_DTDY, D2SIGMA_DY2
    REAL(DP) :: NUM, DEN, DDEN_DY, D2DEN_DY2

    NUM = TWO**(L1+ONE) + Q1
    DEN = y**(-L1)+Q1+(ONE-y)**(-L1)

    SIGMA        = SIGMA_0*H*NUM/DEN

    DSIGMA_DY    = ZERO
    DSIGMA_DT    = ZERO

    D2SIGMA_DY2  = ZERO
    D2SIGMA_DT2  = ZERO
    D2SIGMA_DTDY = ZERO

    IF (derivatives == 0) RETURN

    DDEN_DY = -L1*Y**(-L1-ONE) + L1*(ONE-y)**(-L1-ONE)

    DSIGMA_DY    = SIGMA*(DH_DY/H - DDEN_DY/DEN)
    DSIGMA_DT    = SIGMA*(DH_DT/H)

    IF (derivatives == 1) RETURN

    D2DEN_DY2 = L1*(L1+ONE)*Y**(-L1-TWO) + L1*(L1+ONE)*(ONE-y)**(-L1-TWO)

    D2SIGMA_DY2 = DSIGMA_DY*DSIGMA_DY/SIGMA &
                + SIGMA*( D2H_DY2/H - (DH_DY/H)**TWO  &
                - D2DEN_DY2/DEN + (DDEN_DY/DEN)**TWO )

    D2SIGMA_DT2 = DSIGMA_DT*DSIGMA_DT/SIGMA &
                + SIGMA*( D2H_DT2/H - (DH_DT/H)**TWO )

    D2SIGMA_DTDY = DSIGMA_DY*DSIGMA_DT/SIGMA &
                 + SIGMA*( D2H_DTDY/H - DH_DT*DH_DY/H/H )

    IF (derivatives == 2) RETURN

  END SUBROUTINE SIGMA_SUB

END MODULE Sigma_Mod

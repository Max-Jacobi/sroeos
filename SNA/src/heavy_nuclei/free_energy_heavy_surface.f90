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
MODULE Free_Energy_Heavy_Surface_Mod

  USE Kind_Types_Mod, ONLY : DP, I4B
  USE Physical_Constants_Mod, ONLY : ZERO, ONE, TWO, R_3_2, R_5_2

  IMPLICIT NONE

CONTAINS

  SUBROUTINE FREE_ENERGY_SURFACE (derivatives,u,DELU,DDELU_DU,D2DELU_DU2,&
            BETA,DBETA_Dn,DBETA_DT,DBETA_Dy,D2BETA_Dn2,D2BETA_DnDT,D2BETA_DnDy,&
            D2BETA_DT2,D2BETA_DTDy,D2BETA_Dy2,F_SC,DF_SC_DU,DF_SC_DN,DF_SC_DT, &
            DF_SC_DY,D2F_SC_DU2,D2F_SC_DUDN,D2F_SC_DUDT,D2F_SC_DUDY,D2F_SC_DN2,&
            D2F_SC_DNDT,D2F_SC_DNDY,D2F_SC_DT2,D2F_SC_DTDY,D2F_SC_DY2)
    IMPLICIT NONE

    INTEGER(I4B), INTENT(IN) :: derivatives
    REAL(DP), INTENT(IN) :: u, DELU, DDELU_DU, D2DELU_DU2
    REAL(DP), INTENT(IN) :: BETA, DBETA_DY, DBETA_DN, DBETA_DT
    REAL(DP), INTENT(IN) :: D2BETA_DY2, D2BETA_DN2, D2BETA_DT2
    REAL(DP), INTENT(IN) :: D2BETA_DNDY,D2BETA_DNDT,D2BETA_DTDY
    REAL(DP), INTENT(OUT) :: F_SC, DF_SC_DU, DF_SC_DN, DF_SC_DT, DF_SC_DY
    REAL(DP), INTENT(OUT) :: D2F_SC_DU2, D2F_SC_DUDN, D2F_SC_DUDT, D2F_SC_DUDY
    REAL(DP), INTENT(OUT) :: D2F_SC_DN2, D2F_SC_DNDT, D2F_SC_DNDY
    REAL(DP), INTENT(OUT) :: D2F_SC_DT2, D2F_SC_DTDY, D2F_SC_DY2

    REAL(DP) :: odu

    F_SC = BETA*DELU

    DF_SC_DU = ZERO
    DF_SC_DN = ZERO
    DF_SC_DT = ZERO
    DF_SC_DY = ZERO

    D2F_SC_DU2  = ZERO
    D2F_SC_DUDN = ZERO
    D2F_SC_DUDT = ZERO
    D2F_SC_DUDY = ZERO

    D2F_SC_DN2  = ZERO
    D2F_SC_DNDT = ZERO
    D2F_SC_DNDY = ZERO

    D2F_SC_DT2  = ZERO
    D2F_SC_DTDY = ZERO
    D2F_SC_DY2  = ZERO

!   first derivatives of surface free energy w.r.t. u, n_i, y_I
    IF (derivatives == 0) RETURN

    DF_SC_DU = BETA*DDELU_DU
    DF_SC_DN = DBETA_DN*DELU
    DF_SC_DT = DBETA_DT*DELU
    DF_SC_DY = DBETA_DY*DELU

    IF (derivatives == 1) RETURN

!   second derivatives of surface free energy w.r.t. u, n_i, y_I

    D2F_SC_DU2  = BETA*D2DELU_DU2
    D2F_SC_DUDN = DBETA_DN*DDELU_DU
    D2F_SC_DUDT = DBETA_DT*DDELU_DU
    D2F_SC_DUDY = DBETA_DY*DDELU_DU

    D2F_SC_DN2  = D2BETA_DN2 *DELU
    D2F_SC_DNDT = D2BETA_DNDT*DELU
    D2F_SC_DNDY = D2BETA_DNDY*DELU

    D2F_SC_DT2  = D2BETA_DT2 *DELU
    D2F_SC_DTDY = D2BETA_DTDY*DELU
    D2F_SC_DY2  = D2BETA_DY2 *DELU

    IF (derivatives == 2) RETURN

  END SUBROUTINE FREE_ENERGY_SURFACE

END MODULE Free_Energy_Heavy_Surface_Mod

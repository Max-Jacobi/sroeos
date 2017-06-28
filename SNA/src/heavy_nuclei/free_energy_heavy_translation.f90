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
MODULE Free_Energy_Heavy_Translation_Mod

  USE Kind_Types_Mod, ONLY : DP, I4B
  USE Physical_Constants_Mod, ONLY : ZERO, ONE, TWO, R_3_2, R_5_2
  USE Surface_Properties_Mod, ONLY : fix_heavy_nuclei_size, A00

  IMPLICIT NONE

CONTAINS

  SUBROUTINE FREE_ENERGY_TRANSLATION ( &
      derivatives,n,T,u,DELU,DDELU_DU,D2DELU_DU2,                           &
      H,DH_DT,DH_Dy,D2H_DT2,D2H_Dy2,D2H_DTDy,A,DA_DU,DA_DN,DA_DT,DA_DY,     &
      D2A_DU2,D2A_DUDN,D2A_DUDT,D2A_DUDY,D2A_DN2,D2A_DNDT,D2A_DNDY,D2A_DT2, &
      D2A_DTDY,D2A_DY2,MU,DMU_DU,DMU_DN,DMU_DT,DMU_DY,D2MU_DU2,D2MU_DUDN,   &
      D2MU_DUDT,D2MU_DUDY,D2MU_DN2,D2MU_DNDT,D2MU_DNDY,D2MU_DT2,D2MU_DTDY,  &
      D2MU_DY2, F_TR,DF_TR_DU,DF_TR_DN,DF_TR_DT,DF_TR_DY,D2F_TR_DU2,        &
      D2F_TR_DUDN,D2F_TR_DUDT,D2F_TR_DUDY,D2F_TR_DN2,D2F_TR_DNDT,           &
      D2F_TR_DNDY,D2F_TR_DT2,D2F_TR_DTDY,D2F_TR_DY2 )


    IMPLICIT NONE

    INTEGER(I4B), INTENT(IN) :: derivatives
    REAL(DP), INTENT(IN) :: n, T
    REAL(DP), INTENT(IN) :: u, DELU, DDELU_DU, D2DELU_DU2
    REAL(DP), INTENT(IN) :: H, DH_DT, DH_Dy, D2H_DT2, D2H_Dy2, D2H_DTDy
    REAL(DP), INTENT(IN) :: A, DA_DU, DA_DN, DA_DT, DA_DY
    REAL(DP), INTENT(IN) :: D2A_DU2, D2A_DUDN, D2A_DUDT, D2A_DUDY
    REAL(DP), INTENT(IN) :: D2A_DN2, D2A_DNDT, D2A_DNDY
    REAL(DP), INTENT(IN) :: D2A_DT2, D2A_DTDY, D2A_DY2
    REAL(DP), INTENT(IN) :: MU, DMU_DU, DMU_DN, DMU_DT, DMU_DY
    REAL(DP), INTENT(IN) :: D2MU_DU2, D2MU_DUDN, D2MU_DUDT, D2MU_DUDY
    REAL(DP), INTENT(IN) :: D2MU_DN2, D2MU_DNDT, D2MU_DNDY
    REAL(DP), INTENT(IN) :: D2MU_DT2, D2MU_DTDY, D2MU_DY2

    REAL(DP), INTENT(OUT) :: F_TR, DF_TR_DU, DF_TR_DN, DF_TR_DT, DF_TR_DY
    REAL(DP), INTENT(OUT) :: D2F_TR_DU2, D2F_TR_DUDN, D2F_TR_DUDT, D2F_TR_DUDY
    REAL(DP), INTENT(OUT) :: D2F_TR_DN2, D2F_TR_DNDT, D2F_TR_DNDY
    REAL(DP), INTENT(OUT) :: D2F_TR_DT2, D2F_TR_DTDY, D2F_TR_DY2

    REAL(DP) :: odu

    odu = ONE -u

    IF (fix_heavy_nuclei_size) THEN
      F_TR =  u*odu*n*H*(MU-T)/A00
    ELSE
      F_TR =  u*odu*n*H*(MU-T)/A
    ENDIF

    DF_TR_DU = ZERO
    DF_TR_DN = ZERO
    DF_TR_DT = ZERO
    DF_TR_DY = ZERO

    D2F_TR_DU2  = ZERO
    D2F_TR_DUDN = ZERO
    D2F_TR_DUDT = ZERO
    D2F_TR_DUDY = ZERO

    D2F_TR_DN2  = ZERO
    D2F_TR_DNDT = ZERO
    D2F_TR_DNDY = ZERO

    D2F_TR_DT2  = ZERO
    D2F_TR_DTDY = ZERO
    D2F_TR_DY2  = ZERO

    IF (derivatives == 0) RETURN
!   first derivatives of translational free energy  w.r.t. u, n_i, T, y_i

    DF_TR_DU = (ONE-TWO*u)/u/odu-DA_DU/A+DMU_DU/(MU-T)
    DF_TR_DU = F_TR*DF_TR_DU

    DF_TR_DN = ONE/N - DA_DN/A + DMU_DN/(MU-T)
    DF_TR_DN = F_TR*DF_TR_DN

    DF_TR_DT = DH_DT/H + (DMU_DT-ONE)/(MU-T) - DA_DT/A
    DF_TR_DT = F_TR*DF_TR_DT

    DF_TR_DY = DH_DY/H - DA_DY/A + DMU_DY/(MU-T)
    DF_TR_DY = F_TR*DF_TR_DY

    IF (derivatives == 1) RETURN
!   second derivatives of translational free energy w.r.t. u, n_i, T, y_i
    D2F_TR_DU2  = - TWO/u/odu - ((ONE-TWO*u)/u/odu)**TWO  &
                - D2A_DU2/A + (DA_DU/A)**TWO &
                + D2MU_DU2/(MU-T) - (DMU_DU/(MU-T))**TWO
    D2F_TR_DU2  = DF_TR_DU*DF_TR_DU/F_TR + F_TR*D2F_TR_DU2

    D2F_TR_DUDN = - D2A_DUDN/A + DA_DU*DA_DN/A/A &
                + D2MU_DUDN/(MU-T) - DMU_DU*DMU_DN/(MU-T)**TWO
    D2F_TR_DUDN = DF_TR_DU*DF_TR_DN/F_TR + F_TR*D2F_TR_DUDN

    D2F_TR_DUDT = D2MU_DUDT/(MU-T) - (DMU_DT-ONE)*DMU_DU/(MU-T)**TWO &
                - D2A_DUDT/A + DA_DT*DA_DU/A/A
    D2F_TR_DUDT = DF_TR_DU*DF_TR_DT/F_TR + F_TR*D2F_TR_DUDT

    D2F_TR_DUDY = - D2A_DUDY/A + DA_DU*DA_DY/A/A &
                + D2MU_DUDY/(MU-T) - DMU_DY*DMU_DU/(MU-T)**TWO
    D2F_TR_DUDY = DF_TR_DU*DF_TR_DY/F_TR + F_TR*D2F_TR_DUDY

    D2F_TR_DN2  = - ONE/N/N - D2A_DN2/A + (DA_DN/A)**TWO &
                + D2MU_DN2/(MU-T) - (DMU_DN/(MU-T))**TWO
    D2F_TR_DN2  = DF_TR_DN*DF_TR_DN/F_TR + F_TR*D2F_TR_DN2

    D2F_TR_DNDT = - D2A_DNDT/A + DA_DN*DA_DT/A/A &
                + D2MU_DNDT/(MU-T) - DMU_DN*(DMU_DT-ONE)/(MU-T)**TWO
    D2F_TR_DNDT  = DF_TR_DN*DF_TR_DT/F_TR + F_TR*D2F_TR_DNDT

    D2F_TR_DNDY = - D2A_DNDY/A + DA_DN*DA_DY/A/A &
                + D2MU_DNDY/(MU-T) - DMU_DN*DMU_DY/(MU-T)**TWO
    D2F_TR_DNDY = DF_TR_DN*DF_TR_DY/F_TR + F_TR*D2F_TR_DNDY

    D2F_TR_DT2  = D2H_DT2/H - (DH_DT/H)**TWO &
                + D2MU_DT2/(MU-T) - ((DMU_DT-ONE)/(MU-T))**TWO &
                - D2A_DT2/A + (DA_DT/A)**TWO
    D2F_TR_DT2  = DF_TR_DT*DF_TR_DT/F_TR + F_TR*D2F_TR_DT2

    D2F_TR_DTDY = D2H_DTDY/H -DH_DT*DH_DY/H/H &
                - D2A_DTDY/A + DA_DT*DA_DY/A/A &
                + D2MU_DTDY/(MU-T) - (DMU_DT-ONE)*DMU_DY/(MU-T)**TWO
    D2F_TR_DTDY = DF_TR_DT*DF_TR_DY/F_TR + F_TR*D2F_TR_DTDY

    D2F_TR_DY2  = D2H_DY2/H - (DH_DY/H)**TWO &
                + D2MU_DY2/(MU-T) - (DMU_DY/(MU-T))**TWO &
                - D2A_DY2/A + (DA_DY/A)**TWO
    D2F_TR_DY2  = DF_TR_DY*DF_TR_DY/F_TR + F_TR*D2F_TR_DY2

    IF (derivatives == 2) RETURN

  END SUBROUTINE FREE_ENERGY_TRANSLATION

END MODULE Free_Energy_Heavy_Translation_Mod

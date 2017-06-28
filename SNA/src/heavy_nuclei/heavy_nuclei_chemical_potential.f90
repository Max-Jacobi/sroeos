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
MODULE Mu_Mod

  USE Kind_Types_Mod, ONLY : DP, I4B
  USE Physical_Constants_Mod, ONLY : ZERO, ONE, TWO, R_3_2, R_5_2, PI, &
                                    Hbarc_Square, Mass_n
  USE Surface_Properties_Mod, ONLY : fix_heavy_nuclei_size, A00

  IMPLICIT NONE

CONTAINS

  SUBROUTINE MU_SUB (derivatives,u,n,T,A,DA_DU,DA_DN,DA_DT,DA_DY,           &
      D2A_DU2,D2A_DUDN,D2A_DUDT,D2A_DUDY,D2A_DN2,D2A_DNDT,D2A_DNDY,D2A_DT2, &
      D2A_DTDY,D2A_DY2,MU,DMU_DU,DMU_DN,DMU_DT,DMU_DY,D2MU_DU2,D2MU_DUDN,   &
      D2MU_DUDT,D2MU_DUDY,D2MU_DN2,D2MU_DNDT,D2MU_DNDY,D2MU_DT2,D2MU_DTDY,  &
      D2MU_DY2)

    IMPLICIT NONE

    INTEGER(I4B), INTENT(IN) :: derivatives
    REAL(DP), INTENT(IN) :: u, n, T, A, DA_DU, DA_DN, DA_DT, DA_DY
    REAL(DP), INTENT(IN) :: D2A_DU2, D2A_DUDN, D2A_DUDT, D2A_DUDY
    REAL(DP), INTENT(IN) :: D2A_DN2, D2A_DNDT, D2A_DNDY
    REAL(DP), INTENT(IN) :: D2A_DT2, D2A_DTDY, D2A_DY2

    REAL(DP), INTENT(OUT) :: MU, DMU_DU, DMU_DN, DMU_DT, DMU_DY
    REAL(DP), INTENT(OUT) :: D2MU_DU2, D2MU_DUDN, D2MU_DUDT, D2MU_DUDY
    REAL(DP), INTENT(OUT) :: D2MU_DN2, D2MU_DNDT, D2MU_DNDY
    REAL(DP), INTENT(OUT) :: D2MU_DT2, D2MU_DTDY, D2MU_DY2

    REAL(DP) :: odu, n_Q, Dn_Q_DT

    odu = ONE - u
    n_Q = ((Mass_n*T)/(TWO*PI*Hbarc_Square))**R_3_2

    IF (fix_heavy_nuclei_size) THEN
      MU = T*LOG(u*odu*n/n_Q/A00**R_5_2)
    ELSE
      MU = T*LOG(u*odu*n/n_Q/A**R_5_2)
    ENDIF


    DMU_DU = ZERO
    DMU_DN = ZERO
    DMU_DT = ZERO
    DMU_DY = ZERO

    D2MU_DU2  = ZERO
    D2MU_DUDN = ZERO
    D2MU_DUDT = ZERO
    D2MU_DUDY = ZERO

    D2MU_DN2  = ZERO
    D2MU_DNDT = ZERO
    D2MU_DNDY = ZERO

    D2MU_DT2  = ZERO
    D2MU_DTDY = ZERO
    D2MU_DY2  = ZERO

    IF (derivatives == 0) RETURN

    Dn_Q_DT  = R_3_2*n_Q/T

    DMU_DU = T*( ONE/u - ONE/odu ) - R_5_2*T/A*DA_DU
    DMU_DN = T/n  - R_5_2*T/A*DA_DN
    DMU_DT = MU/T - R_3_2 - R_5_2*T/A*DA_DT
    DMU_DY =      - R_5_2*T/A*DA_DY

    IF (derivatives == 1) RETURN

    D2MU_DU2  = - T*( ONE/u/u + ONE/odu/odu ) &
                - R_5_2*T*( D2A_DU2/A - (DA_DU/A)**TWO )
    D2MU_DUDN = - R_5_2*T/A*( D2A_DUDN - DA_DU*DA_DN/A )
    D2MU_DUDT = ( ONE/u - ONE/odu ) &
                - R_5_2*( DA_DU/A + T*D2A_DUDT/A - T*DA_DU*DA_DT/A/A)
    D2MU_DUDY = - R_5_2*T*( D2A_DUDY/A - DA_DU*DA_DY/A/A )

    D2MU_DN2  = - T/n/n - R_5_2*T*( D2A_DN2/A - (DA_DN/A)**TWO )
    D2MU_DNDT = DMU_DN/T - R_5_2*T*( D2A_DNDT/A - DA_DN*DA_DT/A/A )
    D2MU_DNDY = - R_5_2*T*( D2A_DNDY/A - DA_DN*DA_DY/A/A )

    D2MU_DT2  = - MU/T/T + DMU_DT/T - R_5_2/A*DA_DT &
                - R_5_2*T*( D2A_DT2/A - (DA_DT/A)**TWO )
    D2MU_DTDY = DMU_DY/T - R_5_2*T*( D2A_DTDY/A - DA_DT*DA_DY/A/A )
    D2MU_DY2  = - R_5_2*T*( D2A_DY2/A - (DA_DY/A)**TWO )

    IF (derivatives == 2) RETURN

  END SUBROUTINE MU_SUB

END MODULE Mu_Mod

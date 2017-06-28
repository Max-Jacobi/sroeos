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
MODULE Nuclear_Mass_Mod

  USE Kind_Types_Mod, ONLY : DP, I4B
  USE Physical_Constants_Mod, ONLY : ZERO, ONE, TWO, THREE, FOUR, PI
  USE Surface_Properties_Mod, ONLY : fix_heavy_nuclei_size, A00

  IMPLICIT NONE

CONTAINS

  SUBROUTINE NUCLEAR_MASS_SUB (derivatives,N,R,DR_DU,DR_DN,DR_DT,DR_DY,        &
          D2R_DU2,D2R_DUDN,D2R_DUDT,D2R_DUDY,D2R_DN2,D2R_DNDT,D2R_DNDY,        &
          D2R_DT2,D2R_DTDY,D2R_DY2,A,DA_DU,DA_DN,DA_DT,DA_DY,D2A_DU2,D2A_DUDN, &
          D2A_DUDT,D2A_DUDY,D2A_DN2,D2A_DNDT,D2A_DNDY,D2A_DT2,D2A_DTDY,D2A_DY2)

    IMPLICIT NONE

    INTEGER(I4B), INTENT(IN) :: derivatives
    REAL(DP), INTENT(IN)  :: N, R, DR_DU, DR_DN, DR_DT, DR_DY
    REAL(DP), INTENT(IN)  :: D2R_DU2, D2R_DUDN, D2R_DUDT, D2R_DUDY
    REAL(DP), INTENT(IN)  :: D2R_DN2, D2R_DNDT, D2R_DNDY
    REAL(DP), INTENT(IN)  :: D2R_DT2, D2R_DTDY, D2R_DY2
    REAL(DP), INTENT(OUT) :: A, DA_DU, DA_DN, DA_DT, DA_DY
    REAL(DP), INTENT(OUT) :: D2A_DU2, D2A_DUDN, D2A_DUDT, D2A_DUDY
    REAL(DP), INTENT(OUT) :: D2A_DN2, D2A_DNDT, D2A_DNDY
    REAL(DP), INTENT(OUT) :: D2A_DT2, D2A_DTDY, D2A_DY2

    A = FOUR*PI/THREE*N*R**THREE

    DA_DU = ZERO
    DA_DN = ZERO
    DA_DT = ZERO
    DA_DY = ZERO

    D2A_DU2  = ZERO
    D2A_DUDN = ZERO
    D2A_DUDT = ZERO
    D2A_DUDY = ZERO

    D2A_DN2  = ZERO
    D2A_DNDT = ZERO
    D2A_DNDY = ZERO

    D2A_DT2  = ZERO
    D2A_DTDY = ZERO
    D2A_DY2  = ZERO

    IF (fix_heavy_nuclei_size) THEN
      ! A = A00
      RETURN
    ENDIF

    IF (derivatives == 0) RETURN

    DA_DU = THREE*A/R*DR_DU
    DA_DN = THREE*A/R*DR_DN + A/N
    DA_DT = THREE*A/R*DR_DT
    DA_DY = THREE*A/R*DR_DY

    IF (derivatives == 1) RETURN

    D2A_DU2  = DA_DU*( DA_DU/A + D2R_DU2 /DR_DU - DR_DU/R )
    D2A_DUDN = DA_DU*( DA_DN/A + D2R_DUDN/DR_DU - DR_DN/R )
    D2A_DUDT = DA_DU*( DA_DT/A + D2R_DUDT/DR_DU - DR_DT/R )
    D2A_DUDY = DA_DU*( DA_DY/A + D2R_DUDY/DR_DU - DR_DY/R )

    D2A_DN2  = THREE*A/R*( DA_DN*DR_DN/A + D2R_DN2  - DR_DN**TWO/R ) &
             - A/N/N + DA_DN/N
    D2A_DNDT = THREE*A/R*( DA_DT*DR_DN/A + D2R_DNDT - DR_DN*DR_DT/R ) + DA_DT/N
    D2A_DNDY = THREE*A/R*( DA_DY*DR_DN/A + D2R_DNDY - DR_DN*DR_DY/R ) + DA_DY/N

    D2A_DT2  = DA_DT*( DA_DT/A + D2R_DT2 /DR_DT - DR_DT/R )
    D2A_DTDY = DA_DT*( DA_DY/A + D2R_DTDY/DR_DT - DR_DY/R )
    D2A_DY2  = DA_DY*( DA_DY/A + D2R_DY2 /DR_DY - DR_DY/R )

    IF (derivatives == 2) RETURN

  END SUBROUTINE NUCLEAR_MASS_SUB

END MODULE Nuclear_Mass_Mod

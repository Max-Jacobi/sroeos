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
MODULE H_Mod

  USE Kind_Types_Mod, ONLY : DP, I4B
  USE Physical_Constants_Mod, ONLY : ZERO, ONE, TWO
  USE Surface_Properties_Mod, ONLY : P1 => Surface_p
  USE Critical_Temperature_Mod

  IMPLICIT NONE

CONTAINS

  ! obtain function H = H(y_i,T) and its derivatives
    SUBROUTINE H_SUB (derivatives,y,T,H,DH_DT,DH_DY,D2H_DT2,D2H_DTDY,D2H_DY2)

      IMPLICIT NONE

      INTEGER(I4B), INTENT(IN) :: derivatives
      REAL(DP), INTENT(IN) :: y, T
      REAL(DP), INTENT(OUT) :: H, DH_DT, DH_DY, D2H_DT2, D2H_DY2, D2H_DTDY
      REAL(DP) :: Tc, F, DF_DT, DTc_DY, DF_DTc, DF_DY
      REAL(DP) :: D2F_DT2, D2Tc_Dy2, D2F_DTc2, D2F_DY2, D2F_DTDY

      H        = ZERO
      DH_DT    = ZERO
      DH_Dy    = ZERO
      D2H_DT2  = ZERO
      D2H_Dy2  = ZERO
      D2H_DTDy = ZERO

      IF (T>T_crit) RETURN

      Tc = Temp_c(y)

      IF (T>Tc) RETURN

      F = ONE-(T/Tc)**TWO
      H = F**P1

      IF (derivatives == 0) RETURN

      DF_DT  = - TWO*T/Tc**TWO

      DTc_DY = DTemp_c_Dyi(y)
      DF_DTc = TWO*T*T/Tc**THREE
      DF_DY  = DF_DTc*DTc_DY

      DH_DY = P1*F**(P1-ONE)*DF_DY
      DH_DT = P1*F**(P1-ONE)*DF_DT

      IF (derivatives == 1) RETURN

      D2F_DT2  = DF_DT/T

      D2Tc_Dy2 = D2Temp_c_Dyi2(y)
      D2F_DTc2 = - THREE*DF_DTc/Tc
      D2F_DY2  = D2F_DTc2*DTc_DY*DTc_DY + DF_DTc*D2Tc_Dy2

      D2F_DTDY = -TWO*DF_DT/Tc*DTc_DY

      D2H_DT2  = P1*(P1-ONE)*F**(P1-TWO)*DF_DT**TWO  + P1*F**(P1-ONE)*D2F_DT2
      D2H_DY2  = P1*(P1-ONE)*F**(P1-TWO)*DF_DY**TWO  + P1*F**(P1-ONE)*D2F_DY2
      D2H_DTDY = P1*(P1-ONE)*F**(P1-TWO)*DF_DY*DF_DT + P1*F**(P1-ONE)*D2F_DTDY

      IF (derivatives == 2) RETURN

    END SUBROUTINE H_SUB

END MODULE H_Mod

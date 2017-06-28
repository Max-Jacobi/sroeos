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

!   TODO : MAKE THIS A SUBROUTINE LIKE BETA, SIGMA, AND H
!          THAT DEPENDS ON variable "derivatives"
!
MODULE Critical_Temperature_Mod

  USE Kind_Types_Mod, ONLY : DP, I4B
  USE Physical_Constants_Mod, ONLY : ONE, TWO, THREE, FOUR
  USE Surface_Properties_Mod, ONLY : T_crit, T_crit_coeffs

  IMPLICIT NONE

CONTAINS

  FUNCTION Temp_c (y_i) RESULT(T_c)

    IMPLICIT NONE

    REAL(DP), INTENT(IN) :: y_i
    REAL(DP) :: T_c
    REAL(DP) :: neutron_excess, delta

    neutron_excess = one-two*y_i
    delta = neutron_excess**two

    T_c = T_crit*(  T_crit_coeffs(1) + T_crit_coeffs(2)*delta &
    + T_crit_coeffs(3)*delta**two + T_crit_coeffs(4)*delta**three )

  END FUNCTION Temp_c

  FUNCTION DTemp_c_Dyi (y_i) RESULT(DT_c_Dyi)

    IMPLICIT NONE

    REAL(DP), INTENT(IN) :: y_i
    REAL(DP) :: T_c, DT_c_Dyi, Ddelta_dyi
    REAL(DP) :: neutron_excess, delta

    neutron_excess = one-two*y_i
    delta = neutron_excess**two

    Ddelta_dyi = -FOUR*neutron_excess

    DT_c_Dyi = T_crit*Ddelta_dyi*(T_crit_coeffs(2) &
              + two*T_crit_coeffs(3)*delta &
              + three*T_crit_coeffs(4)*delta**two)

  END FUNCTION DTemp_c_Dyi

  FUNCTION D2Temp_c_Dyi2 (y_i) RESULT(D2T_c_Dyi2)

    IMPLICIT NONE

    REAL(DP), INTENT(IN) :: y_i
    REAL(DP) :: T_c, DT_c_Dyi, D2T_c_Dyi2, Ddelta_dyi, D2delta_dyi2
    REAL(DP) :: neutron_excess, delta

    neutron_excess = one-two*y_i
    delta = neutron_excess**two

    Ddelta_dyi = -FOUR*neutron_excess
    D2delta_dyi2 = 8.0d0

    D2T_c_Dyi2 = T_crit*D2delta_dyi2*(T_crit_coeffs(2) &
      + two*T_crit_coeffs(3)*delta + three*T_crit_coeffs(4)*delta**two) &
      + T_crit*Ddelta_dyi*Ddelta_dyi*(two*T_crit_coeffs(3) &
      + 6.0D0*T_crit_coeffs(4)*delta)

  END FUNCTION D2Temp_c_Dyi2

  FUNCTION Critical_Temperature(y) RESULT(Temp_crit)
!   If temperature is larger than critical temperature nuclei are unstable
!   and, thus, matter is made of uniform gas of p, n, and alphas.
!   Therefore, there is no need to look for nuclei phase for T > Temp_crit
    REAL(DP), INTENT(IN) :: y

    REAL(DP) :: Temp_crit_1, Temp_crit_2, Temp_crit, x
!   Critical T for nuclei with proton fraction y
    Temp_crit_1 = Temp_c(y)
!   Critical T for nuclei with proton fraction x = Max(.16,y+0.04)
!    This is done for cases where proton fraction is too low and no
!    nuclei is possible with said proton fraction
    x = MAX(0.16D0,y+0.04D0)
    Temp_crit_2 = Temp_c(x)
!   MAX T TO TRY TO SOLVE SYSTEM INCLUDING HEAVY NUCLEI
    Temp_crit = MAX(Temp_crit_1,Temp_crit_2)


  END FUNCTION Critical_Temperature

END MODULE Critical_Temperature_Mod

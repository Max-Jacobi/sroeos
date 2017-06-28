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

MODULE Phase_Space_Input_Mod

  USE Kind_Types_Mod, ONLY : DP, I4B

  IMPLICIT NONE

  SAVE

  ! phase space lower and upper bounds
  REAL(DP) :: Yp_min, Yp_max, Yp_spacing, Yp_step
  REAL(DP) :: Log10n_min, Log10n_max, steps_per_decade_in_n, log10n_spacing
  REAL(DP) :: Log10T_min, Log10T_max, steps_per_decade_in_T, log10T_spacing
  INTEGER(I4B) :: steps_in_Yp

  ! initial and final step numbers and total number of steps
  INTEGER(I4B) :: Yp_ini, Yp_fin
  INTEGER(I4B) :: n_ini, n_fin
  INTEGER(I4B) :: T_ini, T_fin

END MODULE Phase_Space_Input_Mod

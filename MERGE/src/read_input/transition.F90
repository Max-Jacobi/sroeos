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

MODULE Transition_Input_Mod

  USE Kind_Types_Mod, ONLY : DP

  IMPLICIT NONE

  SAVE

  ! phase space lower and upper bounds
  REAL(DP) :: n_transition, n_delta, n_tolerance
  REAL(DP) :: Log10nt_max, Log10nt_min

END MODULE Transition_Input_Mod

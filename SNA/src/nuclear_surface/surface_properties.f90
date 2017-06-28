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

MODULE Surface_Properties_Mod

  USE Kind_Types_Mod, ONLY : DP, I4B, LGCL

  IMPLICIT NONE

  SAVE

! Surface properties
  LOGICAL(LGCL) :: fit_surface, fit_temperature
  REAL(DP) :: Surface_sigma, Surface_p, Surface_q, Surface_lambda
! Critical temperature for two phase coexistence
  REAL(DP) :: T_crit
  REAL(DP), DIMENSION(4) :: T_crit_coeffs
  REAL(DP), DIMENSION(100) :: T_table

! Heavy nuclei properties
  LOGICAL(LGCL) :: fix_heavy_nuclei_size
  REAL(DP) :: A0, A00

END MODULE Surface_Properties_Mod

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
MODULE ELE_MOD

  USE Kind_Types_Mod, ONLY : DP

  IMPLICIT NONE

  REAL(DP) :: ele_ener, ele_pres, ele_entr
  REAL(DP) :: ele_denerdt, ele_dpresdt, ele_dentrdt
  REAL(DP) :: ele_denerdd, ele_dpresdd, ele_dentrdd
  REAL(DP) :: ele_denerdy, ele_dpresdy, ele_dentrdy
  ! REAL(DP) :: gam1, etaele, sound

END MODULE ELE_MOD

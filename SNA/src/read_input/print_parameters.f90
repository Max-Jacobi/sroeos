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
MODULE Print_Parameters_Mod

  USE Kind_Types_Mod, ONLY : LGCL

  IMPLICIT NONE

  SAVE

  ! print parameters
  LOGICAL(LGCL) :: print_p, print_s, print_e
  LOGICAL(LGCL) :: print_muh, print_mun, print_mup
  LOGICAL(LGCL) :: print_dpdn, print_dsdn, print_dmuhdn
  LOGICAL(LGCL) :: print_dpdt, print_dsdt, print_dmuhdt
  LOGICAL(LGCL) :: print_dpdy, print_dsdy, print_dmuhdy
  LOGICAL(LGCL) :: print_abar, print_zbar, print_r, print_u
  LOGICAL(LGCL) :: print_xn, print_xp, print_xa, print_xh
  LOGICAL(LGCL) :: print_meff_n, print_meff_p

END MODULE Print_Parameters_Mod

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
MODULE MERGE_ELE_MOD

  USE Kind_Types_Mod, ONLY : DP

  IMPLICIT NONE

  REAL(DP) :: a_in, a_n, dadn, d2adn2
  REAL(DP) :: sna_free, nse_free, free
  REAL(DP) :: sna_denerdd, sna_denerdy, sna_denerdt
  REAL(DP) :: nse_denerdd, nse_denerdy, nse_denerdt

END MODULE MERGE_ELE_MOD

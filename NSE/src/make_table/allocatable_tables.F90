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
MODULE Allocatable_Tables_Mod

  USE Kind_Types_Mod, ONLY : DP, I4B

  IMPLICIT NONE

  SAVE

  INTEGER(I4B) :: pointsrho, pointstemp, pointsye
  INTEGER(I4B) :: n_log10n, n_log10T, n_Yp

  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: p_tab, e_tab, s_tab
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: dpdn_tab, dmudn_tab, dsdn_tab
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: dpdt_tab, dmudt_tab, dsdt_tab
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: dpdy_tab, dmudy_tab, dsdy_tab
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: xn_tab, xp_tab, xa_tab
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: xh_tab, xl_tab
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: ah_tab, zh_tab
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: al_tab, zl_tab
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: muh_tab,mup_tab,mun_tab
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: Meff_n_tab, Meff_p_tab
  REAL(DP), DIMENSION(:), ALLOCATABLE :: log10n_tab, log10T_tab, yp_tab
  LOGICAL, DIMENSION(:,:,:), ALLOCATABLE :: fix

END MODULE Allocatable_Tables_Mod

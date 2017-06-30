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
MODULE MERGE_NSE_MOD

  USE Kind_Types_Mod, ONLY : DP

  IMPLICIT NONE

  SAVE

  REAL(DP), dimension(:,:,:), allocatable :: &
    nse_merge_xa, nse_merge_xh, nse_merge_xp, nse_merge_xn, nse_merge_xl, &
    nse_merge_abar, nse_merge_zbar, nse_merge_albar, nse_merge_zlbar, &
    nse_merge_pres, nse_merge_entr, nse_merge_ener, &
    nse_merge_mu_p, nse_merge_mu_n, nse_merge_mu_e, &
    nse_merge_dsdn, nse_merge_dsdt, nse_merge_dsdy, &
    nse_merge_dpdn, nse_merge_dpdt, nse_merge_dpdy, &
    nse_merge_dedn, nse_merge_dedt, nse_merge_dedy, &
    nse_merge_r, nse_merge_u, nse_merge_meffn, nse_merge_meffp

END MODULE MERGE_NSE_MOD

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
MODULE MERGE_SNA_MOD

  USE Kind_Types_Mod, ONLY : DP

  IMPLICIT NONE

  SAVE

  REAL(DP), dimension(:,:,:), allocatable :: &
    sna_merge_xa, sna_merge_xp, sna_merge_xn, sna_merge_xl, sna_merge_xh, &
    sna_merge_abar, sna_merge_zbar, sna_merge_albar, sna_merge_zlbar, &
    sna_merge_pres, sna_merge_entr, sna_merge_ener, &
    sna_merge_mu_p, sna_merge_mu_n, sna_merge_mu_e, &
    sna_merge_dsdn, sna_merge_dsdt, sna_merge_dsdy, &
    sna_merge_dpdn, sna_merge_dpdt, sna_merge_dpdy, &
    sna_merge_dedn, sna_merge_dedt, sna_merge_dedy, &
    sna_merge_r,  sna_merge_u, sna_merge_meffn, sna_merge_meffp

END MODULE MERGE_SNA_MOD

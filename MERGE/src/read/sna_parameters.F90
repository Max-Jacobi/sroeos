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
MODULE SNA_MOD

  USE Kind_Types_Mod, ONLY : DP

  IMPLICIT NONE

  SAVE

  REAL(DP) :: sna_ener, sna_pres, sna_entr, sna_free
  REAL(DP) :: sna_mu_hat, sna_mu_n, sna_mu_p, sna_mu_e
  REAL(DP) :: sna_dpresdd, sna_dpresdt, sna_dpresdy
  REAL(DP) :: sna_dentrdd, sna_dentrdt, sna_dentrdy
  REAL(DP) :: sna_dmuhdd,  sna_dmuhdt,  sna_dmuhdy
  REAL(DP) :: sna_xn, sna_xp, sna_xa, sna_xh
  REAL(DP) :: sna_abar, sna_zbar, sna_xl, sna_albar, sna_zlbar
  REAL(DP) :: sna_r, sna_u, sna_meffn, sna_meffp

END MODULE SNA_MOD

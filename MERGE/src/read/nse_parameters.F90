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
MODULE NSE_MOD

  USE Kind_Types_Mod, ONLY : DP

  IMPLICIT NONE

  REAL(DP) :: nse_ener, nse_pres, nse_entr
  REAL(DP) :: nse_mu_hat, nse_mu_n, nse_mu_p, nse_mu_e
  REAL(DP) :: nse_dpresdd, nse_dpresdt, nse_dpresdy
  REAL(DP) :: nse_dentrdd, nse_dentrdt, nse_dentrdy
  REAL(DP) :: nse_dmuhdd,  nse_dmuhdt,  nse_dmuhdy
  REAL(DP) :: nse_xn, nse_xp, nse_xa, nse_xh
  REAL(DP) :: nse_abar, nse_zbar, nse_xl, nse_albar, nse_zlbar
  REAL(DP) :: nse_r, nse_u, nse_meffn, nse_meffp

END MODULE NSE_MOD

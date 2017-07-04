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
Module Save_to_Table_Mod

  USE Kind_Types_Mod,        ONLY : I4B, DP
  USE Physical_Constants_Mod, ONLY : ZERO, ONE, TWO
  USE Allocatable_Tables_Mod
  USE Main_output_Mod
  USE, INTRINSIC :: IEEE_ARITHMETIC

  IMPLICIT NONE

CONTAINS

  SUBROUTINE save_to_table (i_n, i_t, i_y, n, Yp, T)

    IMPLICIT NONE

    INTEGER(I4B), INTENT(IN) :: i_n, i_t, i_y
    REAL(DP), INTENT(IN)     :: n, T, Yp

    p_tab(i_n, i_t, i_y) = P
    e_tab(i_n, i_t, i_y) = E/n
    s_tab(i_n, i_t, i_y) = S/n

    dpdn_tab(i_n, i_t, i_y) = DP_DN
    dsdn_tab(i_n, i_t, i_y) = DS_DN/n-S/n**TWO
    dmudn_tab(i_n, i_t, i_y) = Dmu_no_Dn - Dmu_po_Dn

    dpdt_tab(i_n, i_t, i_y) = DP_DT
    dsdt_tab(i_n, i_t, i_y) = DS_DT/n
    dmudt_tab(i_n, i_t, i_y) = Dmu_no_DT - Dmu_po_DT

    dpdy_tab(i_n, i_t, i_y) = DP_DY
    dsdy_tab(i_n, i_t, i_y) = DS_DY/n
    dmudy_tab(i_n, i_t, i_y) =  Dmu_no_Dy - Dmu_po_Dy

    muh_tab(i_n, i_t, i_y)  = mu_no - mu_po
    mun_tab(i_n, i_t, i_y)  = mu_no
    mup_tab(i_n, i_t, i_y)  = mu_po

    meff_n_tab(i_n, i_t, i_y) = Meff_no
    meff_p_tab(i_n, i_t, i_y) = Meff_po

    xn_tab(i_n, i_t, i_y)   = x_no
    xp_tab(i_n, i_t, i_y)   = x_po
    xa_tab(i_n, i_t, i_y)   = x_alpha
    xh_tab(i_n, i_t, i_y)   = x_heavy

    abar_tab(i_n, i_t, i_y) = A_heavy
    zbar_tab(i_n, i_t, i_y) = Z_heavy

    if (ieee_is_nan(A_heavy) .OR. A_heavy < ONE) abar_tab(i_n, i_t, i_y) = 1.d0
    if (ieee_is_nan(Z_heavy) .OR. Z_heavy < Yp ) zbar_tab(i_n, i_t, i_y) = Yp

    u_tab(i_n, i_t, i_y) = u
    r_tab(i_n, i_t, i_y) = r

    if (ieee_is_nan(u)) u_tab(i_n,i_t,i_y) = 0.d0
    if (ieee_is_nan(r)) r_tab(i_n,i_t,i_y) = 0.d0

  end subroutine save_to_table

END MODULE Save_to_Table_Mod

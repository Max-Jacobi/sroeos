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
  USE Physical_Constants_Mod, ONLY : ZERO, ONE, TWO, Neut_Prot_Mass_Diff
  USE Allocatable_Tables_Mod
  USE Main_output_Mod
  USE Skyrme_Bulk_Mod, ONLY : SKYRME_BULK_PROPERTIES
  USE, INTRINSIC :: IEEE_ARITHMETIC

  IMPLICIT NONE

CONTAINS

  SUBROUTINE save_to_table (i_n, i_t, i_y, n, Yp, T)

    IMPLICIT NONE

    INTEGER(I4B), INTENT(IN) :: i_n, i_t, i_y
    REAL(DP), INTENT(IN)     :: n, T, Yp
    REAL(DP) :: b_nn, b_np, b_n, b_y, b_meffn, b_meffp, b_etan, b_etap
    REAL(DP) :: b_taun, b_taup, b_vn, b_vp, b_mun, b_mup, b_p, b_f

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

!   single-particle potential difference dU = U_n - U_p of the (outside)
!   nucleon gas, for beta-process mean-field corrections a la Hempel (2015).
!   v_p from SKYRME_BULK_PROPERTIES contains -(m_n-m_p) since SRO measures
!   energies relative to the neutron mass; remove it so dU is interaction only.
    CALL SKYRME_BULK_PROPERTIES(log10(max(n_no,1.d-100)), &
                                log10(max(n_po,1.d-100)), T, &
                                b_nn, b_np, b_n, b_y, b_meffn, b_meffp, &
                                b_etan, b_etap, b_taun, b_taup, &
                                b_vn, b_vp, b_mun, b_mup, b_p, b_f)
    dU_tab(i_n, i_t, i_y) = b_vn - b_vp - Neut_Prot_Mass_Diff

    xn_tab(i_n, i_t, i_y)   = x_no
    xp_tab(i_n, i_t, i_y)   = x_po
    xa_tab(i_n, i_t, i_y)   = x_alpha
    xh_tab(i_n, i_t, i_y)   = x_heavy

    abar_tab(i_n, i_t, i_y) = A_heavy
    zbar_tab(i_n, i_t, i_y) = Z_heavy

    if (ieee_is_nan(A_heavy) .OR. A_heavy < ONE) abar_tab(i_n, i_t, i_y) = 1.d0
    if (ieee_is_nan(Z_heavy) .OR. Z_heavy < Yp ) zbar_tab(i_n, i_t, i_y) = Yp

    u_tab(i_n, i_t, i_y) = u
!   r (heavy-nucleus / Wigner-Seitz radius) is only computed when a
!   non-uniform solution is chosen. u and r are THREADPRIVATE module
!   variables (Main_Output_Mod), so on the uniform branch r is never
!   reassigned and the table keeps whatever the previously processed point
!   left behind. u does not have this problem: the uniform branch sets
!   log10_u = -290 in GET_OUTPUT, so u self-clears. Measured consequence of
!   the stale r: values 0-4.73 fm at points with no nucleus, overlapping the
!   genuine 3.27-6.97 fm range, so only x_heavy tells them apart -- and the
!   table then depends on the density grid's extent (the first points of each
!   T row inherit the last, lowest-n point of the previous row).
    IF (u > 1.0d-100) THEN
      r_tab(i_n, i_t, i_y) = r
    ELSE
      r_tab(i_n, i_t, i_y) = 0.d0
    ENDIF

    if (ieee_is_nan(u)) u_tab(i_n,i_t,i_y) = 0.d0
    if (ieee_is_nan(r)) r_tab(i_n,i_t,i_y) = 0.d0

  end subroutine save_to_table

END MODULE Save_to_Table_Mod

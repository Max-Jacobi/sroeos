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
MODULE Main_output_Mod

  USE Kind_Types_Mod, ONLY : DP

  SAVE

  REAL(DP) :: F_o, F_i, F_alpha, F_TR, F_SC
  REAL(DP) :: u, r, x_no, x_po, x_alpha, x_heavy
  REAL(DP) :: n_no, n_po, n_ni, n_pi, n_alpha, n_heavy
  REAL(DP) :: A_heavy, Z_heavy, rad, F, P, S, E
  REAL(DP) :: DF_Dn, DF_Dy, DF_DT, DP_Dn, DP_Dy, DP_DT
  REAL(DP) :: DS_Dn, DS_Dy, DS_DT, DE_Dn, DE_Dy, DE_DT
  REAL(DP) :: mu_no, mu_po, Meff_no, Meff_po
  REAL(DP) :: Dmu_no_DT, Dmu_no_Dn, Dmu_no_Dy
  REAL(DP) :: Dmu_po_DT, Dmu_po_Dn, Dmu_po_Dy
  REAL(DP) :: dLog10_n_no_dn, dLog10_n_po_dn, dLog10_u_dn
  REAL(DP) :: dLog10_n_no_dT, dLog10_n_po_dT, dLog10_u_dT

!$OMP THREADPRIVATE( F_o, F_i, F_alpha, F_TR, F_SC )
!$OMP THREADPRIVATE( u, r, x_no, x_po, x_alpha, x_heavy )
!$OMP THREADPRIVATE( n_no, n_po, n_ni, n_pi, n_alpha, n_heavy )
!$OMP THREADPRIVATE( A_heavy, Z_heavy, rad, F, P, S, E )
!$OMP THREADPRIVATE( DF_Dn, DF_Dy, DF_DT, DP_Dn, DP_Dy, DP_DT )
!$OMP THREADPRIVATE( DS_Dn, DS_Dy, DS_DT, DE_Dn, DE_Dy, DE_DT )
!$OMP THREADPRIVATE( mu_no, mu_po, Meff_no, Meff_po )
!$OMP THREADPRIVATE( Dmu_no_DT, Dmu_no_Dn, Dmu_no_Dy )
!$OMP THREADPRIVATE( Dmu_po_DT, Dmu_po_Dn, Dmu_po_Dy )
!$OMP THREADPRIVATE( dLog10_n_no_dn, dLog10_n_po_dn, dLog10_u_dn )
!$OMP THREADPRIVATE( dLog10_n_no_dT, dLog10_n_po_dT, dLog10_u_dT )

END MODULE Main_output_Mod
